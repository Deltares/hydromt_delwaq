"""Implement demission model class."""

import logging
from os.path import isfile, join
from pathlib import Path
from typing import Dict, List

import geopandas as gpd
import numpy as np
import pandas as pd
import pyproj
import xarray as xr
from hydromt import hydromt_step
from hydromt.model import Model, processes
from hydromt.model.components import GeomsComponent

from hydromt_delwaq.components import (
    DelwaqHydromapsComponent,
    DelwaqStaticdataComponent,
    DemissionConfigComponent,
    DemissionForcingComponent,
    DemissionGeometryComponent,
)
from hydromt_delwaq.utils import DATADIR, write_waqgeomfile
from hydromt_delwaq.workflows import (
    config,
    emissions,
    forcing,
    geometry,
    monitoring,
    roads,
    segments,
)

__all__ = ["DemissionModel"]
__hydromt_eps__ = ["DemissionModel"]  # core entrypoints
logger = logging.getLogger(__name__)


class DemissionModel(Model):
    """Demission model class."""

    name: str = "demission"

    _MAPS = {
        "flwdir": "ldd",
        "lndslp": "slope",
        "rivmsk": "river",
        "strord": "streamorder",
        "soil_theta_s": "porosity",
    }
    _FORCING = {
        "temp": "tempair",
        "temp_dew": "temp_dew",
        "ssr": "rad",
        "wind": "vwind",
        "tcc": "cloudfrac",
    }

    def __init__(
        self,
        root: str | Path | None = None,
        mode: str = "w",
        config_filename: str | Path | None = None,
        data_libs: List[str | Path] | None = None,
    ):
        """Initialize a model.

        Parameters
        ----------
        root : str, optional
            Model root, by default None
        mode : {'r','r+','w', 'w+'}, optional
            read/append/write mode, by default "w"
        config_filename : str or Path, optional
            Model simulation configuration file, by default None.
            Note that this is not the HydroMT model setup configuration file!
        data_libs : List[str, Path], optional
            List of data catalog configuration files, by default None
        """
        if config_filename is None:
            config_filename = "emission.inp"

        components = {
            "config": DemissionConfigComponent(
                self,
                filename=str(config_filename),
            ),
            "staticdata": DelwaqStaticdataComponent(self),
            "hydromaps": DelwaqHydromapsComponent(self, region_component="staticdata"),
            "geoms": GeomsComponent(
                self, filename="geoms/{name}.geojson", region_component="staticdata"
            ),
            "geometry": DemissionGeometryComponent(
                self, filename="config/B7_geometry.inc"
            ),
            "forcing": DemissionForcingComponent(self, region_component="staticdata"),
        }

        super().__init__(
            root,
            components=components,
            mode=mode,
            region_component="staticdata",
            data_libs=data_libs,
        )

        # d-emission specific
        self.timestepsecs = 86400

    ## Properties
    # Components
    @property
    def config(self) -> DemissionConfigComponent:
        """Return the config component."""
        return self.components["config"]

    @property
    def staticdata(self) -> DelwaqStaticdataComponent:
        """Return the staticdata component."""
        return self.components["staticdata"]

    @property
    def hydromaps(self) -> DelwaqHydromapsComponent:
        """Return the hydromaps component."""
        return self.components["hydromaps"]

    @property
    def geoms(self) -> GeomsComponent:
        """Return the geoms component."""
        return self.components["geoms"]

    @property
    def geometry(self) -> DemissionGeometryComponent:
        """Return the geometry component."""
        return self.components["geometry"]

    @property
    def forcing(self) -> DemissionForcingComponent:
        """Return the forcing component."""
        return self.components["forcing"]

    @hydromt_step
    def setup_basemaps(
        self,
        region: Dict,
        maps: List[str] = ["rivmsk", "lndslp", "strord"],
    ):
        """
        Prepare demission schematization using the hydromodel region and resolution.

        Maps used and derived from the hydromodel are stored in a specific
        hydromodel attribute. Depending on the global option ``type``, build a
        one-substance D-Emission ("EM") model case.
        No specific arguments are needed for this function.

        Adds model layers:

        * **ldd** hydromap: flow direction [-]
        * **modelmap** hydromap: mask map [bool]
        * **ptid** hydromap: unique ID of Delwaq segments [-]
        * **river** map: river mask map [-]
        * **pointer** poi: delwaq pointer between segments
        * **B3_nrofseg** config: number of segments
        * **B3_attributes** config: delwaq complete attributes
        * **B4_nrofexch** config: number of exchanges
        * **B5_boundlist** config: list of boundaries and types
        * **B7_surf** config: surface data

        Parameters
        ----------
        region : dict
            Dictionary describing region of interest.
            Currently supported format is {'wflow': 'path/to/wflow_model'}
        maps: list of str
            List of variables from hydromodel to add to grid.
            By default ['rivmsk', 'lndslp', 'strord'].
        """
        # Initialise hydromodel from region
        hydromodel = processes.region.parse_region_other_model(region)

        if hydromodel.name != "wflow_sbm":
            raise NotImplementedError(
                "Demission build function only implemented for wflow_sbm base model."
            )

        logger.info("Preparing EM basemaps from hydromodel.")

        ### Select and build hydromaps from model ###
        # Initialise hydromaps
        ds_hydro = segments.hydromaps(hydromodel)
        # Add to hydromaps
        rmdict = {k: v for k, v in self._MAPS.items() if k in ds_hydro.data_vars}
        self.hydromaps.set(ds_hydro.rename(rmdict))

        # Build segment ID and segment ID down and add to hydromaps
        nrofseg, da_ptid, da_ptiddown = segments.pointer(
            self.hydromaps.data, build_pointer=False
        )
        self.hydromaps.set(da_ptid, name="ptid")
        self.hydromaps.set(da_ptiddown, name="ptiddown")

        ### Initialise grid with segment ID down, streamorder, river and slope ###
        ds_stat = segments.maps_from_hydromodel(
            hydromodel,
            maps=maps,
        )
        ds_stat["ptiddown"] = self.hydromaps.data["ptiddown"].squeeze(drop=True)

        rmdict = {k: v for k, v in self._MAPS.items() if k in ds_stat.data_vars}
        self.staticdata.set(ds_stat.rename(rmdict))
        self.staticdata.data.coords["mask"] = self.hydromaps.data["mask"]

        ### Geometry data ###
        hydromodel_grid = hydromodel.components.get(
            hydromodel._region_component_name
        ).data

        geometry_data = geometry.compute_geometry(
            ds=hydromodel_grid,
            mask=self.staticdata.data["mask"],
            fpaved_name=hydromodel._MAPS["soil_compacted_fraction"],
            fopenwater_name=hydromodel._MAPS["land_water_fraction"],
        )
        self.geometry.set(
            pd.DataFrame(
                geometry_data, columns=(["TotArea", "fPaved", "fUnpaved", "fOpenWater"])
            )
        )

        ### Config ###
        configs = config.base_config(
            nrofseg=nrofseg,
            nrofexch=0,
            layer_name="EM",
            add_surface=True,
        )
        self.config.update(data=configs)

    @hydromt_step
    def setup_monitoring(
        self,
        mon_points: str | Path = None,
        mon_areas: str = None,
    ):
        """Prepare Delwaq monitoring points and areas options.

        Adds model layers:

        * **monitoring_points** map: monitoring points segments
        * **monitoring_areas** map: monitoring areas ID
        * **B2_nrofmon** config: number of monitoring points and areas

        Parameters
        ----------
        mon_points : {'segments', 'data_source', 'path to station location'}
            Either source from DataCatalog, path to a station location dataset
            or if segments, all segments are monitored.
        mon_areas : str {'subcatch', 'riverland'}
            Either subcatch from hydromodel or 'riverland' to split river and land
            cells.
        """
        logger.info("Setting monitoring points and areas")

        # Read monitoring points source
        if mon_points is not None:
            logger.info(f"Reading monitoring points from {mon_points}")
            if mon_points == "segments":
                nb_points, monpoints = monitoring.monitoring_points_from_dataarray(
                    self.hydromaps.data["ptid"]
                )
                gdf = None
            else:
                kwargs = {}
                if isfile(mon_points):
                    kwargs["metadata"] = {"crs": self.crs}
                gdf = self.data_catalog.get_geodataframe(
                    mon_points,
                    geom=self.basins,
                    # assert_gtype="Point",
                    source_kwargs=kwargs,
                )
                (
                    nb_points,
                    monpoints,
                    gdf,
                ) = monitoring.monitoring_points_from_geodataframe(
                    gdf,
                    ds_like=self.hydromaps.data,
                )

            # Add to grid
            self.staticdata.set(monpoints.rename("monpoints"))
            # Add to geoms if mon_points is not segments
            if gdf is not None:
                self.geoms.set(gdf, name="monpoints")
        else:
            logger.info("No monitoring points set in the config file, skipping")
            nb_points = 0

        # Monitoring areas domain
        if mon_areas is not None:
            nb_areas, monareas = monitoring.monitoring_areas(
                mon_areas,
                ds=self.hydromaps.data,
            )
            self.staticdata.set(monareas, name="monareas")

            # Add to geoms
            gdf_areas = (
                self.staticdata.data["monareas"].astype(np.int32).raster.vectorize()
            )
            self.geoms.set(gdf_areas, name="monareas")
        else:
            logger.info("No monitoring areas set in the config file, skipping")
            nb_areas = 0

        # Config
        lines_ini = {
            "l1": f";{nb_points} monitoring points",
            "l2": f";{nb_areas} monitoring areas",
            "l3": f"{nb_points+nb_areas} ; nr of monitoring points/areas",
        }
        for option in lines_ini:
            self.config.set(f"B2_nrofmon.{option}", lines_ini[option])

    @hydromt_step
    def setup_emission_raster(
        self,
        emission_fn: str | Path | xr.DataArray,
        scale_method: str = "average",
        fillna_method: str = "zero",
        fillna_value: int = 0.0,
        area_division: bool = False,
    ):
        """Prepare one or several emission map from raster data.

        Adds model layer:

        * **emission_fn** map: emission data map

        Parameters
        ----------
        emission_fn : {'GHS-POP_2015'...}
            Name of raster emission map source.
        scale_method : str {'nearest', 'average', 'mode'}
            Method for resampling
        fillna_method : str {'nearest', 'zero', 'value'}
            Method to fill NaN values. Either nearest neighbour, zeros or user defined
            value.
        fillna_value : float
            If fillna_method is set to 'value', NaNs in the emission maps will be
            replaced by this value.
        area_division : boolean
            If needed do the resampling in cap/m2 (True) instead of cap (False)
        """
        logger.info(f"Preparing '{emission_fn}' map.")
        # process raster emission maps
        da = self.data_catalog.get_rasterdataset(
            emission_fn, geom=self.region, buffer=2
        )
        ds_emi = emissions.emission_raster(
            da=da,
            ds_like=self.staticdata.data,
            method=scale_method,
            fillna_method=fillna_method,
            fillna_value=fillna_value,
            area_division=area_division,
        )
        ds_emi = ds_emi.to_dataset(name=emission_fn)
        self.staticdata.set(ds_emi)

    @hydromt_step
    def setup_emission_vector(
        self,
        emission_fn: str | xr.DataArray,
        col2raster: str = "",
        rasterize_method: str = "value",
    ):
        """Prepare emission map from vector data.

        Adds model layer:

        * **emission_fn** map: emission data map

        Parameters
        ----------
        emission_fn : {'GDP_world'...}
            Name of raster emission map source.
        col2raster : str
            Name of the column from the vector file to rasterize.
            Can be left empty if the selected method is set to "fraction".
        rasterize_method : str
            Method to rasterize the vector data. Either {"value", "fraction", "area"}.
            If "value", the value from the col2raster is used directly in the raster.
            If "fraction", the fraction of the grid cell covered by the vector file is
            returned.
            If "area", the area of the grid cell covered by the vector file is returned.
        """
        logger.info(f"Preparing '{emission_fn}' map.")
        gdf_org = self.data_catalog.get_geodataframe(
            emission_fn,
            geom=self.basins,
        )
        if gdf_org.empty:
            logger.warning(
                f"No shapes of {emission_fn} found within region, "
                "setting to default value."
            )
            ds_emi = self.hydromaps.data["basins"].copy() * 0.0
            ds_emi.attrs.update(_FillValue=0.0)
        else:
            ds_emi = emissions.emission_vector(
                gdf=gdf_org,
                ds_like=self.staticdata.data,
                col_name=col2raster,
                method=rasterize_method,
                mask_name="mask",
            )
        ds_emi = ds_emi.to_dataset(name=emission_fn)
        self.staticdata.set(ds_emi)

    @hydromt_step
    def setup_emission_mapping(
        self,
        region_fn: str | Path | gpd.GeoDataFrame,
        mapping_fn: str | Path | None = None,
    ):
        """
        Derive several emission maps based on administrative boundaries.

        For several emission types administrative classes ('fid' column) are
        remapped to model parameter values based on lookup tables. The data is
        remapped at its original resolution and then resampled to the model
        resolution based using the average value, unless noted differently.

        Adds model layers:

        * **region_fn** map: emission data with classification from source_name [-]
        * **emission factor X** map: emission data from mapping file to classification

        Parameters
        ----------
        region_fn : {["gadm_level1", "gadm_level2", "gadm_level3"]}
            Name or list of names of data source in data_sources.yml file.

            * Required variables: ['ID']
        mapping_fn : str, optional
            Path to the emission mapping file corresponding to region_fn.
        """
        logger.info(f"Preparing region_fn related parameter maps for {region_fn}.")
        if mapping_fn is None:
            logger.warning("Using default mapping file.")
            mapping_fn = join(DATADIR, "admin_bound", f"{region_fn}_mapping.csv")
        # process emission factor maps
        gdf_org = self.data_catalog.get_geodataframe(
            region_fn, geom=self.basins, dst_crs=self.crs
        )
        # Rasterize the GeoDataFrame to get the areas mask of
        # administrative boundaries with their ids
        gdf_org["ID"] = gdf_org["ID"].astype(np.int32)
        # make sure index_col always has name fid in source dataset (use rename in
        # data_sources.yml or local_sources.yml to rename column used for mapping
        ds_admin = self.hydromaps.data.raster.rasterize(
            gdf_org,
            col_name="ID",
            nodata=0,
            all_touched=True,
            dtype=None,
            sindex=False,
        )

        # add admin_bound map
        ds_admin_maps = emissions.admin(
            da=ds_admin,
            ds_like=self.staticdata.data,
            source_name=region_fn,
            fn_map=mapping_fn,
        )
        rmdict = {k: v for k, v in self._MAPS.items() if k in ds_admin_maps.data_vars}
        self.staticdata.set(ds_admin_maps.rename(rmdict))

    @hydromt_step
    def setup_roads(
        self,
        roads_fn: str | Path | gpd.GeoDataFrame,
        highway_list: str | List[str],
        country_list: str | List[str],
        non_highway_list: str | List[str] | None = None,
        country_fn: str | Path | gpd.GeoDataFrame | None = None,
    ):
        """Prepare roads statistics needed for emission modelling.

        Adds model layers:

        * **km_highway_country** map: emission data with for each grid cell the total km
          of highway for the country the grid cell is in [km highway/country]
        * **km_other_country** map: emission data with for each grid cell the total km
          of non highway roads for the country the grid cell is in [km other
          road/country]
        * **km_highway_cell** map: emission data containing highway length per cell
          [km highway/cell]
        * **km_other_cell** map: emission data containing non-highway length per cell
          [km non-highway/cell]
        * **emission factor X** map: emission data from mapping file to classification

        Parameters
        ----------
        roads_fn: str
            Name of road data source in data_sources.yml file.

            * Required variables: ['road_type']

            * Optional variable: ['length', 'country_code']
        highway_list: str or list of str
            List of highway roads in the type variable of roads_fn.
        non_highway_list: str or list of str, optional.
            List of non highway roads in the type variable of roads_fn. If not provided
            takes every roads except the ones in highway_list.
        country_list: str or list of str, optional.
            List of countries for the model area in country_fn and optionally in
            country_code variable of roads_fn.
        country_fn: str, optional.
            Name of country boundaries data source in data_sources.yml file.

            * Required variables: ['country_code']

        """
        ### Country statistics ###
        # Mask the road data with countries of interest
        if country_fn is not None:
            if not isinstance(country_list, list):
                country_list = [country_list]
            gdf_country = self.data_catalog.get_geodataframe(country_fn)
            gdf_country = gdf_country.astype({"country_code": "str"})
            gdf_country = gdf_country.iloc[
                np.isin(gdf_country["country_code"], country_list)
            ]
            # Read the roads data and mask with country geom
            gdf_roads = self.data_catalog.get_geodataframe(
                roads_fn, geom=gdf_country, variables=["road_type"]
            )
            # Compute country statistics
            ds_country = roads.roads_emissions_country(
                gdf_roads=gdf_roads,
                gdf_country=gdf_country,
                ds_like=self.staticdata.data,
                highway_list=highway_list,
                non_highway_list=non_highway_list,
            )
            # Add to grid
            self.staticdata.set(ds_country)

        ### Road statistics per segment ###
        # Filter road gdf with model mask instead of whole country
        mask = self.staticdata.data["mask"].astype("int32").raster.vectorize()
        gdf_roads = self.data_catalog.get_geodataframe(
            roads_fn, geom=mask, variables=["road_type"]
        )
        # Compute segment statistics
        ds_segments = roads.roads_emissions_segments(
            gdf_roads=gdf_roads,
            ds_like=self.staticdata.data,
            highway_list=highway_list,
            non_highway_list=non_highway_list,
        )
        self.staticdata.set(ds_segments)

    @hydromt_step
    def setup_hydrology_forcing(
        self,
        hydro_forcing_fn: str | Path | xr.Dataset,
        starttime: str,
        endtime: str,
        timestepsecs: int,
        include_transport: bool = True,
    ):
        """Prepare Demission hydrological fluxes.

        Without transport, the fluxes required are rainfall (precip), runoff
        from paved (runPav) and unpaved (runUnp) areas, and infiltration
        (infilt). With transport, additional fluxes are required for exfiltration
        (exfilt), root zone soil moisture (vwcproot), and land (q_land) and subsurface
        (q_ss) runoff.

        All fluxes are in m3/s except for soil moisture which is in % (volumetric water
        content of the soil pores for the root zone ie volumetric water content of the
        root zone divided by the porosity).

        If fluxes are given in mm, they are converted to m3/s using the grid cell area.
        The unit of the fluxes can be defined in the data catalog in the `attrs`
        properties.

        Adds model layers:

        * **hydrology** dynmap: fluxes for D-Emission.
        * **B7_hydrology** config: names of fluxes included in hydrology.bin
        * **B1_timestamp** config: timestamp at the beginning of the simulation.
        * **B2_outputtimes** config: start and end timestamp for the delwaq outputs.
        * **B2_sysclock** config: model timestep.
        * **B2_timers** config: timers info (start, end, timestep...).

        Parameters
        ----------
        hydro_forcing_fn : {'name in local_sources.yml'}
            Either None, or name in a local or global data_sources.yml file.

            * Required variables without transport: ['time', 'precip', 'runPav',
              'runUnp', 'infilt']
            * Required variables with transport: ['time', 'precip', 'runPav',
              'runUnp', 'infilt', 'exfilt*', 'vwcproot', 'q_land', 'q_ss']

        startime : str
            Timestamp of the start of Delwaq simulation. Format: YYYY-mm-dd HH:MM:SS
        endtime : str
            Timestamp of the end of Delwaq simulation.Format: YYYY-mm-dd HH:MM:SS
        timestepsecs : int
            Model timestep in seconds.
        include_transport : bool, optional
            If False (default), only use the vertical fluxes for emission [precip,
            runPav, runUnp, infilt, totflw].
            If True, includes additional fluxes for land and subsurface transport
            [precip, runPav, runUnp, infilt, exfilt, vwcproot, q_land, q_ss].
        """
        logger.info(f"Setting dynamic data from hydrology source {hydro_forcing_fn}.")
        # read dynamic data
        # Normally region extent is by default exactly the same as hydromodel
        ds_in = self.data_catalog.get_rasterdataset(
            hydro_forcing_fn,
            geom=self.region,
            time_range=(starttime, endtime),
        )

        # Update model timestepsecs attribute
        self.timestepsecs = timestepsecs

        # Compute hydrology forcing
        ds = forcing.hydrology_forcing_em(
            ds=ds_in,
            ds_model=self.hydromaps.data,
            timestepsecs=timestepsecs,
            include_transport=include_transport,
        )
        # Rename xdim and ydim
        if ds.raster.x_dim != self.staticdata.data.raster.x_dim:
            ds = ds.rename({ds.raster.x_dim: self.staticdata.data.raster.x_dim})
        if ds.raster.y_dim != self.staticdata.data.raster.y_dim:
            ds = ds.rename({ds.raster.y_dim: self.staticdata.data.raster.y_dim})

        # Add to forcing
        self.forcing.set(ds)

        # Add timers info to config
        time_config = config.time_config(
            starttime=starttime,
            endtime=endtime,
            timestepsecs=timestepsecs,
            sysclock_format="days",
        )
        self.config.update(data=time_config)

        # Add hydrology variables info to config
        if include_transport:
            l2 = (
                "Rainfall RunoffPav RunoffUnp Infiltr "
                "Exfiltr vwcproot Overland Subsurface"
            )
        else:
            l2 = "Rainfall RunoffPav RunoffUnp Infiltr TotalFlow"
        lines_ini = {
            "l1": "SEG_FUNCTIONS",
            "l2": l2,
        }
        for option in lines_ini:
            self.config.set(f"B7_hydrology.{option}", lines_ini[option])

    @hydromt_step
    def setup_climate_forcing(
        self,
        climate_fn: str | Path | xr.Dataset,
        starttime: str,
        endtime: str,
        timestepsecs: int,
        climate_vars: list = ["temp", "temp_dew", "ssr", "wind10_u", "wind10_v", "tcc"],
        temp_correction: bool = False,
        dem_forcing_fn: str = None,
    ):
        """Prepare Delwaq climate fluxes.

        Adds model layers:

        * **climate** dynmap: climate fuxes for climate_vars
        * **B7_climate** config: names of climate fluxes included in climate.bin

        Parameters
        ----------
        climate_fn : str
            Climate data

            * Required variables: variables listed in climate_vars
        startime : str
            Timestamp of the start of Delwaq simulation. Format: YYYY-mm-dd HH:MM:SS
        timestepsecs: int
            Delwaq model timestep in seconds.
        endtime : str
            Timestamp of the end of Delwaq simulation.Format: YYYY-mm-dd HH:MM:SS
        climate_vars : str list, optional
            List of climate_vars to consider. By default ["temp", "temp_dew", "ssr",
            "wind_u", "wind_v", "tcc"]
        temp_correction : bool, optional
            If True temperature are corrected using elevation lapse rate,
            by default False.
        dem_forcing_fn : str, default None
            Elevation data source with coverage of entire meteorological forcing domain.
            If temp_correction is True and dem_forcing_fn is provided this is used in
            combination with elevation at model resolution to correct the temperature.
        """
        logger.info(f"Setting dynamic data from climate source {climate_fn}.")

        # read dynamic data
        ds = self.data_catalog.get_rasterdataset(
            climate_fn,
            geom=self.region,
            buffer=2,
            variables=climate_vars,
            time_range=(starttime, endtime),
            single_var_as_array=False,
        )
        dem_forcing = None
        if dem_forcing_fn is not None:
            dem_forcing = self.data_catalog.get_rasterdataset(
                dem_forcing_fn,
                geom=ds.raster.box,  # clip dem with forcing bbox for full coverage
                buffer=2,
                variables=["elevtn"],
            ).squeeze()

        ds_out = forcing.climate_forcing(
            ds=ds,
            ds_model=self.hydromaps.data,
            timestepsecs=timestepsecs,
            temp_correction=temp_correction,
            dem_forcing=dem_forcing,
        )

        # Rename and add to forcing
        rmdict = {k: v for k, v in self._FORCING.items() if k in ds_out.data_vars}
        ds_out = ds_out.rename(rmdict)
        # Rename xdim and ydim
        if ds_out.raster.x_dim != self.staticdata.data.raster.x_dim:
            ds_out = ds_out.rename(
                {ds_out.raster.x_dim: self.staticdata.data.raster.x_dim}
            )
        if ds_out.raster.y_dim != self.staticdata.data.raster.y_dim:
            ds_out = ds_out.rename(
                {ds_out.raster.y_dim: self.staticdata.data.raster.y_dim}
            )
        self.forcing.set(ds_out)

        # Update config
        lines_ini = {
            "l1": "SEG_FUNCTIONS",
            "l2": " ".join(str(x) for x in ds_out.data_vars),
        }
        for option in lines_ini:
            self.config.set(f"B7_climate.{option}", lines_ini[option])

    # I/O
    @hydromt_step
    def read(self):
        """Read the complete model schematization and configuration from file."""
        self.config.read()
        self.geoms.read()
        self.hydromaps.read()
        self.staticdata.read()
        self.geometry.read()
        self.forcing.read()

    @hydromt_step
    def write(self):
        """Write the complete model schematization and configuration to file."""
        logger.info(f"Write model data to {self.root}")
        self.write_data_catalog()
        _ = self.config.data  # try to read default if not yet set
        self.staticdata.write()
        self.geoms.write()
        self.config.write()
        self.hydromaps.write()
        self.geometry.write()
        self.forcing.write()

    @hydromt_step
    def write_waqgeom(self):
        """Write Delwaq netCDF geometry file (config/B3_waqgeom.nc)."""
        # Add waqgeom.nc file to allow Delwaq to save outputs in nc format
        fname = join(self.root.path, "config", "B3_waqgeom.nc")
        write_waqgeomfile(self.hydromaps.data, fname)

    ## D-Emission specific properties

    @property
    def basins(self):
        """Derive basins from geoms or hydromaps."""
        if "basins" in self.geoms.data:
            gdf = self.geoms.data["basins"]
        elif "basins" in self.hydromaps.data:
            gdf = self.hydromaps.data["basins"].raster.vectorize()
            gdf.crs = pyproj.CRS.from_user_input(self.crs)
            self.geoms.set(gdf, name="basins")
        return gdf

    @property
    def nrofseg(self):
        """Number of segments."""
        # from config
        nseg = self.config.get_value("B3_nrofseg.l1", fallback="0 ; nr of segments")
        nseg = int(nseg.split(";")[0])
        return nseg

    @property
    def nrofexch(self):
        """Number of exchanges."""
        nexch = self.nrofseg * 5
        return nexch

    @property
    def fluxes(self):
        """List of fluxes."""
        # from config
        fl = self.config.get_value(
            "B7_hydrology.l2",
            fallback=(
                "Rainfall RunoffPav RunoffUnp Infiltr " + "Exfiltr Overland Subsurface"
            ),
        )
        fl = fl.split(" ")
        return fl
