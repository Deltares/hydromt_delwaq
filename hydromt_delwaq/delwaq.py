"""Implement delwaq model class."""

import logging
from os.path import isfile, join
from pathlib import Path
from typing import Dict, List

import geopandas as gpd
import numpy as np
import pyproj
import xarray as xr
import xugrid as xu
from hydromt import hydromt_step
from hydromt.model import Model, processes
from hydromt.model.components import GeomsComponent

from hydromt_delwaq.components import (
    DelwaqConfigComponent,
    DelwaqForcingComponent,
    DelwaqHydromapsComponent,
    DelwaqPointerComponent,
    DelwaqStaticdataComponent,
)
from hydromt_delwaq.utils import DATADIR
from hydromt_delwaq.workflows import config, emissions, forcing, segments

__all__ = ["DelwaqModel"]
__hydromt_eps__ = ["DelwaqModel"]  # core entrypoints
logger = logging.getLogger(f"hydromt.{__name__}")


class DelwaqModel(Model):
    """Delwaq model class."""

    name: str = "delwaq"

    _MAPS = {
        # "mask": "modelmap",
        "flwdir": "ldd",
        "lndslp": "slope",
        "N": "manning",
        "rivmsk": "river",
        "strord": "streamorder",
        "thetaS": "porosity",
        "reslocs": "reservoirs",
        "lakelocs": "lakes",
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
        config_fn : str or Path, optional
            Model simulation configuration file, by default None.
            Note that this is not the HydroMT model setup configuration file!
        data_libs : List[str, Path], optional
            List of data catalog configuration files, by default None
        """
        components = {
            "config": DelwaqConfigComponent(
                self,
                filename=str(config_filename),
            ),
            "staticdata": DelwaqStaticdataComponent(self),
            "pointer": DelwaqPointerComponent(self),
            "hydromaps": DelwaqHydromapsComponent(self, region_component="staticdata"),
            "forcing": DelwaqForcingComponent(self, region_component="staticdata"),
            "geoms": GeomsComponent(
                self, filename="geoms/{name}.geojson", region_component="staticdata"
            ),
        }

        super().__init__(
            root,
            components=components,
            mode=mode,
            region_component="staticdata",
            data_libs=data_libs,
        )

        # delwaq specific
        self.timestepsecs = 86400

    ## Properties
    # Components
    @property
    def config(self) -> DelwaqConfigComponent:
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
    def pointer(self) -> DelwaqPointerComponent:
        """Return the pointer component."""
        return self.components["pointer"]

    ## SETUP METHODS
    @hydromt_step
    def setup_basemaps(
        self,
        region: Dict,
        mask: str = "basins",
        surface_water: str = "sfw",
        boundaries: List[str] = ["bd"],
        fluxes: List[str] = ["sfw>sfw", "bd>sfw"],
        maps: List[str] = ["rivmsk", "lndslp", "strord", "N"],
    ):
        """
        Prepare delwaq schematization using the hydromodel region and resolution.

        Maps used and derived from the hydromodel are stored in a specific\
        hydromodel attribute. Build a D-Water Quality ("WQ") case.\
        The pointer will be created based on the ``fluxes`` list
        between surface water and boundaries.

        Adds model layers:

        * **ldd** hydromap: flow direction [-]
        * **modelmap** hydromap: mask map [bool]
        * **ptid** hydromap: unique ID of Delwaq segments [-]
        * **river** map: river mask map [-]
        * **surface** map: surface map [m2]
        * **length** map: length map [m]
        * **width** map: width map [m]
        * **pointer** poi: delwaq pointer between segments
        * **B3_nrofseg** config: number of segments
        * **B3_attributes** config: delwaq complete attributes
        * **B4_nrofexch** config: number of exchanges
        * **B5_boundlist** config: list of boundaries and types
        * **B7_flow** config: list of flow names
        * **B7_volume** config: list of volume names

        Parameters
        ----------
        region : dict
            Dictionary describing region of interest.
            Currently supported format is {'wflow': 'path/to/wflow_model'}
        mask : str, optional
            Mask used to define Delwaq segments. Either "basins" (default) or "rivers".
        surface_water : str, optional
            Name of the surface water layer. Used to identify fluxes to and from surface
            waters in the fluxes list. By default 'sfw'.
        boundaries: list of str, optional
            List of names of boundaries to include. By default a unique boundary called
            'bd'.
        fluxes: list of str
            List of fluxes to include between surface water/boundaries. Name convention
            is '{surface_water}>{boundary_name}' for a flux from the surface water to a
            boundary, ex 'sfw>bd'. By default ['sfw>sfw', 'bd>sfw'] for runoff and
            inwater.
            Names in the fluxes list should match name in the hydrology_fn source in
            setup_hydrology_forcing.
        maps: list of str
            List of variables from hydromodel to add to grid.
            By default ['rivmsk', 'lndslp', 'strord', 'N'].
        """
        # Initialise hydromodel from region
        hydromodel = processes.region.parse_region_other_model(region)

        if hydromodel.name != "wflow_sbm":
            raise NotImplementedError(
                "Delwaq build function only implemented for wflow_sbm base model."
            )

        logger.info("Preparing WQ basemaps from hydromodel.")

        ### Select and build hydromaps from model ###
        # Initialise hydromaps
        ds_hydro = segments.hydromaps(hydromodel, mask=mask)
        # Add to hydromaps
        rmdict = {k: v for k, v in self._MAPS.items() if k in ds_hydro.data_vars}
        self.hydromaps.set(ds_hydro.rename(rmdict))

        # Build segment ID and add to hydromaps
        # Prepare delwaq pointer file for WAQ simulation (not needed with EM)
        nrofseg, da_ptid, da_ptiddown, pointer, bd_id, bd_type = segments.pointer(
            self.hydromaps.data,
            build_pointer=True,
            surface_water=surface_water,
            boundaries=boundaries,
            fluxes=fluxes,
        )
        self.hydromaps.set(da_ptid.rename("ptid"))
        self.hydromaps.set(da_ptiddown.rename("ptiddown"))
        # Initialise pointer object with schematisation attributes
        self.pointer.set("pointer", value=pointer)
        self.pointer.set("nrofseg", value=nrofseg)
        self.pointer.set("surface_water", value=surface_water)
        self.pointer.set("boundaries", value=bd_type)
        self.pointer.set("fluxes", value=fluxes)
        ### Initialise grid river and slope ###
        ds_stat = segments.maps_from_hydromodel(
            hydromodel,
            maps=maps,
        )
        ds_stat["mask"] = self.hydromaps.data["mask"]
        ds_stat = ds_stat.set_coords("mask")
        rmdict = {k: v for k, v in self._MAPS.items() if k in ds_stat.data_vars}
        self.staticdata.set(ds_stat.rename(rmdict))

        ### Add geometry ###
        ds_geom = segments.geometrymaps(
            hydromodel,
        )
        self.staticdata.set(ds_geom)

        ### Config ###
        configs = config.base_config(
            nrofseg=nrofseg,
            nrofexch=self.nrofexch,
            layer_name=surface_water,
            add_surface=False,
            boundaries=bd_id,
            boundaries_type=bd_type,
            fluxes=fluxes,
            volumes=[surface_water],
        )
        self.config.update(data=configs)

    @hydromt_step
    def setup_monitoring(
        self,
        mon_points: str = None,
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
        monpoints = None
        mv = -999
        # Read monitoring points source
        if mon_points is not None:
            if mon_points == "segments":
                monpoints = self.hydromaps.data["ptid"]
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
                gdf = gdf.to_crs(self.crs)
                if gdf.index.size == 0:
                    logger.warning(
                        f"No {mon_points} gauge locations found within domain"
                    )
                else:
                    gdf.index.name = "index"
                    monpoints = self.staticdata.data.raster.rasterize(
                        gdf, col_name="index", nodata=mv
                    )
            logger.info(f"Gauges locations read from {mon_points}")
            # Add to grid
            self.staticdata.set(monpoints.rename("monpoints"))
            # Number of monitoring points
            points = monpoints.values.flatten()
            points = points[points != mv]
            nb_points = len(points)
            # Add to geoms if mon_points is not segments
            if mon_points != "segments":
                self.geoms.set(gdf, name="monpoints")
        else:
            logger.info("No monitoring points set in the config file, skipping")
            nb_points = 0

        # Monitoring areas domain
        monareas = None
        if mon_areas is not None:
            if mon_areas == "subcatch":  # subcatch
                basins = xr.where(
                    self.hydromaps.data["basins"] > 0,
                    self.hydromaps.data["basins"],
                    mv,
                ).astype(np.int32)
                # Number or monitoring areas
                areas = basins.values.flatten()
                areas = areas[areas != mv]
                nb_areas = len(np.unique(areas))
                monareas = basins
            elif mon_areas == "riverland":  # riverland
                # seperate areas for land cells (1) and river cells (2)
                lr_areas = xr.where(
                    self.hydromaps.data["river"],
                    2,
                    xr.where(self.hydromaps.data["basins"], 1, mv),
                ).astype(np.int32)
                # Apply the current model mask
                lr_areas = lr_areas.where(self.hydromaps.data["mask"], mv)
                lr_areas.raster.set_nodata(mv)
                # Number or monitoring areas
                areas = lr_areas.values.flatten()
                areas = areas[areas != mv]
                nb_areas = len(np.unique(areas))
                monareas = lr_areas
            else:
                raise ValueError(
                    f"Unknown monitoring area type {mon_areas}. "
                    "Valid options are 'subcatch' or 'riverland'."
                )

            # Add to grid
            monareas.attrs.update(_FillValue=mv)
            monareas.attrs["mon_areas"] = mon_areas
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
    def setup_hydrology_forcing(
        self,
        hydro_forcing_fn: str | Path | xr.Dataset,
        starttime: str,
        endtime: str,
        timestepsecs: int = 86400,
        add_volume_offset: bool = True,
        min_volume: float = 0.1,
        override: List = [],
    ):
        """Prepare Delwaq hydrological fluxes.

        As the fluxes order should precisely match the pointer defined in
        setup_basemaps, the variables names in ``hydro_forcing_fn`` should match names
        defined in the ``fluxes`` argument of setup_basemaps. These names should also
        have been saved in the file config/B7_flow.inc.

        If several sub-variables in ``hydro_forcing_fn`` need to be summed up to get
        the expected flux in pointer, they can be named {flux_name_in_pointer}_{number}
        (eg "sfw>sfw_1" and "sfw>sfw_2" to get "sfw>sfw") and the function will sum them
        on the fly. To remove (- instead of +) use unit_mult attribute of the data
        catalogwith -1 for the sub-variables of interest.
        To override rather than sum fluxes, use the ``override`` argument and name of
        the concerned flux (eg "sfw>sfw"). The flux are overwritten (excluding nodata)
        in the order of the list, so the last one will be the one used.

        Unit conversions are possible from mm/area to m3/s for fluxes, volumes should be
        provided directly in m3. For conversion from mm to m3/s, it is possible to
        specify over wich surface area the mm are calculated. If 'mm' (default), the
        cellarea is assumed. Else, you can use 'mm/{surfacearea}' where {surfacearea}
        should be a map available in hydromaps (rivarea, lakearea, resarea).

        Adds model layers:

        * **B1_timestamp** config: timestamp at the beginning of the simulation.
        * **B2_outputtimes** config: start and end timestamp for the delwaq outputs.
        * **B2_sysclock** config: model timestep.
        * **B2_timers** config: timers info (start, end, timstep...).
        * **flow** dynmap: water fluxes [m3/s]
        * **volume** dynmap: water volumes [m3]

        Parameters
        ----------
        hydro_forcing_fn : {'None', 'name in local_sources.yml'}
            Either None, or name in a local or global data_sources.yml file.
            Names of fluxes should match those as set in the ``fluxes`` list of
            setup_basemaps methods (pointer creation).

            * Required variables (to be deprecated): ['time', 'run', 'vol' or 'lev',
              'inwater']

        startime : str
            Timestamp of the start of Delwaq simulation. Format: YYYY-mm-dd HH:MM:SS
        endtime : str
            Timestamp of the end of Delwaq simulation.Format: YYYY-mm-dd HH:MM:SS
        timestepsecs : int
            Model timestep in seconds. By default 86400.
        add_volume_offset : bool, optional
            Delwaq needs water volumes at the beginning of the timestep.
            In some models, like wflow, volumes are written at the end of the timestep
            and therefore an offset of one timestep needs to be added for consistency.
        min_volume : float, optional
            Add a minimum value to all the volumes in Delwaq to avoid zeros. Default is
            0.1m3.
        override : list, optional
            List of fluxes/volumes to override for non nodata values rather than sum.
            The last one (based on _number) has priority.
        """
        if not self.data_catalog.contains_source(hydro_forcing_fn):
            logger.warning(
                f"None or Invalid source {hydro_forcing_fn} ",
                "skipping setup_hydrology_forcing.",
            )
            return
        logger.info(f"Setting dynamic data from hydrology source {hydro_forcing_fn}.")
        # read dynamic data
        # Normally region extent is by default exactly the same as hydromodel
        ds = self.data_catalog.get_rasterdataset(hydro_forcing_fn, geom=self.region)

        # Prepare hydrology forcing
        ds_out = forcing.hydrology_forcing(
            ds=ds,
            ds_model=self.hydromaps.data,
            timestepsecs=timestepsecs,
            time_tuple=(starttime, endtime),
            fluxes=self.fluxes,
            surface_water=self.surface_water,
            add_volume_offset=add_volume_offset,
            min_volume=min_volume,
            override=override,
        )
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

        # Update model timestepsecs attribute
        self.timestepsecs = timestepsecs

        # Update config times
        time_config = config.time_config(
            starttime=starttime,
            endtime=endtime,
            timestepsecs=timestepsecs,
        )
        self.config.update(data=time_config)

    @hydromt_step
    def setup_sediment_forcing(
        self,
        sediment_fn: str | Path | xr.Dataset,
        starttime: str,
        endtime: str,
        timestepsecs: int = 86400,
        particle_class: List[str] = ["IM1", "IM2", "IM3", "IM4", "IM5"],
    ):
        """Prepare Delwaq sediment fluxes.

        Adds model layers:

        * **sediment** dynmap: sediment particles from land erosion (fErodIM*)
          [g/timestep]
        * **B7_sediment** config: names of sediment fluxes included in sediment.bin

        Parameters
        ----------
        sediment_fn : {'None', 'name in local_sources.yml'}
            Either None, or name in a local or global data_sources.yml file. Can contai
            several particule sizes (IM1, IM2 etc)

            * Required variables: ['ErodIM*']
        startime : str
            Timestamp of the start of Delwaq simulation. Format: YYYY-mm-dd HH:MM:SS
        endtime : str
            Timestamp of the end of Delwaq simulation.Format: YYYY-mm-dd HH:MM:SS
        timestepsecs: int
            Delwaq model timestep in seconds. By default 86400 for daily.
        particle_class : str list, optional
            List of particle classes to consider. By default 5 classes: ['IM1', 'IM2',
            'IM3', 'IM4', 'IM5']
        """
        logger.info(f"Setting dynamic data from sediment source {sediment_fn}.")
        # select particle classes
        sed_vars = [f"Erod{x}" for x in particle_class]
        # read data
        ds = self.data_catalog.get_rasterdataset(
            sediment_fn,
            geom=self.region,
            time_tuple=(starttime, endtime),
            variables=sed_vars,
        )
        # Prepare sediment forcing
        ds_out = forcing.sediment_forcing(
            ds=ds,
            ds_model=self.hydromaps.data,
            timestepsecs=timestepsecs,
        )
        # Rename xdim and ydim
        if ds_out.raster.x_dim != self.staticdata.data.raster.x_dim:
            ds_out = ds_out.rename(
                {ds_out.raster.x_dim: self.staticdata.data.raster.x_dim}
            )
        if ds_out.raster.y_dim != self.staticdata.data.raster.y_dim:
            ds_out = ds_out.rename(
                {ds_out.raster.y_dim: self.staticdata.data.raster.y_dim}
            )
        # Add to forcing
        self.forcing.set(ds_out)

        # Update config
        lines_ini = {
            "l1": "SEG_FUNCTIONS",
            "l2": " ".join(str(x) for x in sed_vars),
        }
        for option in lines_ini:
            self.config.set(f"B7_sediment.{option}", lines_ini[option])

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
            time_tuple=(starttime, endtime),
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
            emission_fn, geom=self.basins, dst_crs=self.crs
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

    # I/O
    @hydromt_step
    def read(self):
        """Read the complete model schematization and configuration from file."""
        self.config.read()
        self.staticdata.read()
        self.hydromaps.read()
        self.geoms.read()
        self.pointer.read()
        # self.forcing.read()

    @hydromt_step
    def write(self):
        """Write the complete model schematization and configuration to file."""
        logger.info(f"Write model data to {self.root.path}")
        # if in r, r+ mode, only write updated components
        if not self.root.is_writing_mode():
            logger.warning("Cannot write in read-only mode")
            return
        self.write_data_catalog()
        _ = self.config.data  # try to read default if not yet set

        self.staticdata.write()
        self.hydromaps.write()
        self.geoms.write()
        self.pointer.write()
        self.forcing.write()
        self.config.write()

    ## DELWAQ specific data and methods

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
        """Fast accessor to nrofseg property of pointer."""
        if "nrofseg" in self.pointer.data:
            nseg = self.pointer.data["nrofseg"]
        else:
            # from config
            nseg = self.config.get_value("B3_nrofseg.l1", fallback="0 ; nr of segments")
            nseg = int(nseg.split(";")[0])
            self.pointer.set("nrofseg", value=nseg)
        return nseg

    @property
    def nrofexch(self):
        """Fast accessor to nrofexch property of pointer."""
        if "nrofexch" in self.pointer.data:
            nexch = self.pointer.data["nrofexch"]
        elif "pointer" in self.pointer.data:
            nexch = self.pointer.data["pointer"].shape[0]
            self.pointer.set("nrofexch", value=nexch)
        else:
            # from config
            nexch = self.config.get_value(
                "B4_nrofexch.l1",
                fallback="0 0 0 ; x, y, z direction",
            )
            nexch = int(nexch.split(" ")[0])
            self.pointer.set("nrofexch", value=nexch)
        return nexch

    @property
    def surface_water(self):
        """Fast accessor to surface_water name property of pointer."""
        if "surface_water" in self.pointer.data:
            sfw = self.pointer.data["surface_water"]
        else:
            # from config
            nl = 7
            sfw = self.config.get_value(
                f"B3_attributes.l{nl}",
                fallback="     1*01 ; sfw",
            )
            sfw = sfw.split(";")[1][1:]
            self.pointer.set("surface_water", value=sfw)
        return sfw

    @property
    def fluxes(self):
        """Fast accessor to fluxes property of pointer."""
        if "fluxes" in self.pointer.data:
            fl = self.pointer.data["fluxes"]
        else:
            # from config
            fl = self.config.get_value(
                "B7_flow.l2",
                fallback="sfw>sfw inw>sfw",
            )
            fl = fl.split(" ")
            self.pointer.set("fluxes", value=fl)
        return fl

    @hydromt_step
    def write_waqgeom(self):
        """Write Delwaq netCDF geometry file (config/B3_waqgeom.nc)."""
        # Add waqgeom.nc file to allow Delwaq to save outputs in nc format
        # TODO: Update for several layers
        logger.info("Writting waqgeom.nc file")
        ptid = self.hydromaps.data["ptid"].copy()
        # For now only 1 comp is supported
        ptid = ptid.squeeze(drop=True)
        if len(ptid.dims) != 2:
            raise ValueError("Only 2D (1 comp) supported for waqgeom.nc")
        ptid_mv = ptid.raster.nodata
        # PCR cell id's start at 1, we need it zero based, and NaN set to -1
        ptid = xr.where(ptid == ptid_mv, -1, ptid - 1)
        np_ptid = ptid.values
        # Wflow map dimensions
        m, n = np_ptid.shape
        # Number of segments in horizontal dimension
        int(np.max(np_ptid)) + 1  # because of the zero based
        # Get LDD map
        np_ldd = self.hydromaps.data["ldd"].squeeze(drop=True).values
        np_ldd[np_ldd == self.hydromaps.data["ldd"].raster.nodata] = 0

        # print("Input DataArray: ")
        # print(ptid.dims, ptid)
        ptid = xr.where(ptid == -1, np.nan, ptid)
        if ptid.raster.y_dim != "y":
            ptid = ptid.rename({ptid.raster.y_dim: "y"})
        if ptid.raster.x_dim != "x":
            ptid = ptid.rename({ptid.raster.x_dim: "x"})
        # ptid = ptid.rename({"lat": "y", "lon": "x"})
        da_ptid = xu.UgridDataArray.from_structured2d(ptid)
        da_ptid = da_ptid.dropna(dim=da_ptid.ugrid.grid.face_dimension)
        da_ptid.ugrid.set_crs(crs=self.crs)  # "EPSG:4326"
        uda_waqgeom = da_ptid.ugrid.to_dataset(
            optional_attributes=True
        )  # .to_netcdf("updated_ugrid.nc")
        uda_waqgeom.coords["projected_coordinate_system"] = -2147483647

        epsg_nb = int(self.crs.to_epsg())
        uda_waqgeom["projected_coordinate_system"].attrs.update(
            dict(
                epsg=epsg_nb,
                grid_mapping_name="Unknown projected",
                longitude_of_prime_meridian=0.0,
                inverse_flattening=298.257223563,
                epsg_code=f"{self.crs}",
                value="value is equal to EPSG code",
            )
        )

        grid = da_ptid.ugrid.grid

        bounds = grid.face_node_coordinates
        x_bounds = bounds[..., 0]
        y_bounds = bounds[..., 1]

        name_x = "mesh2d_face_x_bnd"
        name_y = "mesh2d_face_y_bnd"
        uda_waqgeom["mesh2d_face_x"].attrs["bounds"] = name_x
        uda_waqgeom["mesh2d_face_y"].attrs["bounds"] = name_y
        uda_waqgeom["mesh2d_face_x"].attrs["units"] = "degrees_east"
        uda_waqgeom["mesh2d_face_y"].attrs["units"] = "degrees_north"
        uda_waqgeom[name_x] = xr.DataArray(
            x_bounds, dims=(grid.face_dimension, "mesh2d_nMax_face_nodes")
        )
        uda_waqgeom[name_y] = xr.DataArray(
            y_bounds, dims=(grid.face_dimension, "mesh2d_nMax_face_nodes")
        )

        # Write the waqgeom.nc file
        fname = join(self.root.path, "config", "B3_waqgeom.nc")
        uda_waqgeom.to_netcdf(path=fname, mode="w")
