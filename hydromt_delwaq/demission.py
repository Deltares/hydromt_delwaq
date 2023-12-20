"""Implement demission model class"""

import os
from os.path import join, isfile, basename
import glob
import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
import pyproj
import logging
import struct
from datetime import datetime
import time as t
from typing import List

import hydromt
from hydromt import workflows


from hydromt_wflow.wflow import WflowModel
from .delwaq import DelwaqModel

from .workflows import emissions, segments, roads
from . import DATADIR

__all__ = ["DemissionModel"]

logger = logging.getLogger(__name__)

# specify pcraster map types for non scalar (float) data types
PCR_VS_MAP = {
    "modelmap": "bool",
    "basins": "nominal",
    "ldd": "ldd",
    "ptid": "ordinal",
}


class DemissionModel(DelwaqModel):
    """This is the demission model class"""

    _NAME = "demission"
    _CONF = "demission.inp"
    _DATADIR = DATADIR
    _CF = dict()  # configreader kwargs
    _GEOMS = {}
    _MAPS = {
        "basmsk": "modelmap",
        "flwdir": "ldd",
        "lndslp": "slope",
        "N": "manning",
        "rivmsk": "river",
        "strord": "streamorder",
        "thetaS": "porosity",
    }
    _FOLDERS = [
        "hydromodel",
        "staticdata",
        "geoms",
        "config",
        "dynamicdata",
        "fews",
    ]

    def __init__(
        self,
        root=None,
        mode="w",
        config_fn=None,
        hydromodel_name="wflow",
        hydromodel_root=None,
        data_libs=None,
        deltares_data=False,
        logger=logger,
    ):
        super().__init__(
            root=root,
            mode=mode,
            config_fn=config_fn,
            hydromodel_name=hydromodel_name,
            hydromodel_root=hydromodel_root,
            data_libs=data_libs,
            deltares_data=deltares_data,
            logger=logger,
        )

    def setup_basemaps(
        self,
        region,
        maps=["rivmsk", "lndslp", "strord"],
    ):
        """
        Setup the demission model schematization using the hydromodel region and resolution.

        Maps used and derived from the hydromodel are stored in a specific
        hydromodel attribute. Depending on the global option ``type``, build a
        one-substance D-Emission ("EM") model case.
        No specific arguments are needed for this function.

        Adds model layers:

        * **ldd** hydromap: flow direction [-]
        * **modelmap** hydromap: mask map [bool]
        * **ptid** hydromap: unique ID of Delwaq segments [-]
        * **streamorder** map: Strahler stream order map. [-]
        * **slope** map: slope map [-]
        * **pointer** poi: delwaq pointer between segments

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
        kind, region = workflows.parse_region(region, logger=self.logger)
        if kind != "model":
            raise ValueError("Delwaq model can only be built from 'model' region.")
        hydromodel = region["mod"]
        if hydromodel._NAME != "wflow":
            raise NotImplementedError(
                "Demission build function only implemented for wflow base model."
            )
        else:
            self.hydromodel_name = hydromodel._NAME
            self.hydromodel_root = hydromodel.root

        self.logger.info(f"Preparing EM basemaps from hydromodel.")

        ### Select and build hydromaps from model ###
        # Initialise hydromaps
        ds_hydro = segments.hydromaps(hydromodel)
        # Add mask
        da_mask = ds_hydro["basmsk"]
        ds_hydro = ds_hydro.drop_vars(["rivmsk", "basmsk"])
        ds_hydro.coords["mask"] = da_mask
        ds_hydro["modelmap"] = da_mask.astype(np.int32)
        ds_hydro["modelmap"].raster.set_nodata(0)
        # Add to hydromaps
        rmdict = {k: v for k, v in self._MAPS.items() if k in ds_hydro.data_vars}
        self.set_hydromaps(ds_hydro.rename(rmdict))

        # Build segment ID and segment ID down and add to hydromaps
        nrofseg, da_ptid, da_ptiddown = segments.pointer(
            self.hydromaps, build_pointer=False
        )
        self.set_hydromaps(da_ptid, name="ptid")
        self.set_hydromaps(da_ptiddown, name="ptiddown")

        ### Initialise grid with segment ID down, streamorder, river and slope ###
        ds_stat = segments.maps_from_hydromodel(
            hydromodel, compartments=self.compartments, maps=maps
        )
        ds_stat["ptiddown"] = self.hydromaps["ptiddown"].squeeze(drop=True)

        rmdict = {k: v for k, v in self._MAPS.items() if k in ds_stat.data_vars}
        self.set_grid(ds_stat.rename(rmdict))
        self.grid.coords["mask"] = self.hydromaps["mask"]

        ### Geometry data ###
        surface = emissions.gridarea(hydromodel.grid)
        surface.raster.set_nodata(np.nan)
        geom = surface.rename("surface").to_dataset()
        # For EM type build a pointer like object and add to self.geometry
        geom["fPaved"] = hydromodel.grid["PathFrac"]
        geom["fOpenWater"] = hydromodel.grid["WaterFrac"]
        geom["fUnpaved"] = (
            (geom["fPaved"] * 0.0 + 1.0) - geom["fPaved"] - geom["fOpenWater"]
        )
        geom["fUnpaved"] = xr.where(geom["fUnpaved"] < 0.0, 0.0, geom["fUnpaved"])

        mask = self.grid["mask"].values.flatten()
        for dvar in ["surface", "fPaved", "fUnpaved", "fOpenWater"]:
            data = geom[dvar].values.flatten()
            data = data[mask].reshape(nrofseg, 1)
            if dvar == "surface":
                geometry = data
            else:
                geometry = np.hstack((geometry, data))
        self._geometry = pd.DataFrame(
            geometry, columns=(["TotArea", "fPaved", "fUnpaved", "fOpenWater"])
        )

        ### Config ###
        # For now, only one compartment in EM and WQ
        nrofcomp = 1
        # B3_nrofseg
        lines_ini = {
            "l1": f"{nrofseg} ; nr of segments",
        }
        for option in lines_ini:
            self.set_config("B3_nrofseg", option, lines_ini[option])
        # B3_attributes
        lines_ini = {
            "l1": "      ; DELWAQ_COMPLETE_ATTRIBUTES",
            "l2": " 1    ; one block with input",
            "l3": " 2    ; number of attributes, they are :",
            "l4": "     1     2",
            "l5": " 1    ; file option in this file",
            "l6": " 1    ; option without defaults",
            "l7": f"     {int(nrofseg/nrofcomp)}*01 ; EM",
            "l8": " 0    ; no time dependent attributes",
        }
        for option in lines_ini:
            self.set_config("B3_attributes", option, lines_ini[option])
        # B4_nrofexch
        nrexch = 0
        lines_ini = {
            "l1": f"{nrexch} 0 0 ; x, y, z direction",
        }
        for option in lines_ini:
            self.set_config("B4_nrofexch", option, lines_ini[option])
        # B5_boundlist
        lines_ini = {
            "l1": f";'NodeID' 'Number' 'Type'",
        }
        for option in lines_ini:
            self.set_config("B5_boundlist", option, lines_ini[option])
        # B7_Surf
        lines_ini = {
            "l1": f"PARAMETERS Surf ALL DATA {nrofseg}*1.0",
        }
        for option in lines_ini:
            self.set_config("B7_surf", option, lines_ini[option])

    def setup_emission_raster(
        self,
        emission_fn: str,
        scale_method: str = "average",
        fillna_method: str = "zero",
        fillna_value: int = 0.0,
        area_division: bool = False,
    ):
        """Setup one or several emission map from raster data.

        Adds model layer:

        * **emission_fn** map: emission data map

        Parameters
        ----------
        emission_fn : {'GHS-POP_2015'...}
            Name of raster emission map source.
        scale_method : str {'nearest', 'average', 'mode'}
            Method for resampling
        fillna_method : str {'nearest', 'zero', 'value'}
            Method to fill NaN values. Either nearest neighbour, zeros or user defined value.
        fillna_value : float
            If fillna_method is set to 'value', NaNs in the emission maps will be replaced by this value.
        area_division : boolean
            If needed do the resampling in cap/m2 (True) instead of cap (False)
        comp_emi: str
            Name of the model compartment recaiving the emission data (by default surface water 'sfw').
        """
        self.logger.info(f"Preparing '{emission_fn}' map.")
        # process raster emission maps
        da = self.data_catalog.get_rasterdataset(
            emission_fn, geom=self.region, buffer=2
        )
        ds_emi = emissions.emission_raster(
            da=da,
            ds_like=self.grid,
            method=scale_method,
            fillna_method=fillna_method,
            fillna_value=fillna_value,
            area_division=area_division,
            logger=self.logger,
        )
        self.set_grid(ds_emi.rename(emission_fn))

    def setup_emission_vector(
        self,
        emission_fn: str,
        col2raster: str = "",
        rasterize_method: str = "value",
    ):
        """Setup emission map from vector data.

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
            Method to rasterize the vector data. Either {"value", "fraction"}.
            If "value", the value from the col2raster is used directly in the raster.
            If "fraction", the fraction of the grid cell covered by the vector file is returned.
        comp_emi: str
            Name of the model compartment recaiving the emission data (by default surface water 'sfw').
        """
        self.logger.info(f"Preparing '{emission_fn}' map.")
        gdf_org = self.data_catalog.get_geodataframe(
            emission_fn, geom=self.basins, dst_crs=self.crs
        )
        if gdf_org.empty:
            self.logger.warning(
                f"No shapes of {emission_fn} found within region, setting to default value."
            )
            ds_emi = self.hydromaps["basins"].copy() * 0.0
            ds_emi.attrs.update(_FillValue=0.0)
        else:
            ds_emi = emissions.emission_vector(
                gdf=gdf_org,
                ds_like=self.grid,
                col_name=col2raster,
                method=rasterize_method,
                mask_name="mask",
                logger=self.logger,
            )
        self.set_grid(ds_emi.rename(emission_fn))

    def setup_emission_mapping(
        self,
        region_fn,
        mapping_fn=None,
    ):
        """This component derives several emission maps based on administrative
        boundaries.

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
        comp_emi: str
            Name of the model compartment recaiving the emission data (by default surface water 'sfw').
        """
        self.logger.info(
            f"Preparing administrative boundaries related parameter maps for {region_fn}."
        )
        if mapping_fn is None:
            self.logger.warning(f"Using default mapping file.")
            mapping_fn = join(DATADIR, "admin_bound", f"{region_fn}_mapping.csv")
        # process emission factor maps
        gdf_org = self.data_catalog.get_geodataframe(
            region_fn, geom=self.basins, dst_crs=self.crs
        )
        # Rasterize the GeoDataFrame to get the areas mask of administrative boundaries with their ids
        gdf_org["ID"] = gdf_org["ID"].astype(np.int32)
        # make sure index_col always has name fid in source dataset (use rename in data_sources.yml or
        # local_sources.yml to rename column used for mapping (INDEXCOL), if INDEXCOL name is not fid:
        # rename:
        #   INDEXCOL: fid)\
        ds_admin = self.hydromaps.raster.rasterize(
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
            ds_like=self.grid,
            source_name=region_fn,
            fn_map=mapping_fn,
            logger=self.logger,
        )
        rmdict = {k: v for k, v in self._MAPS.items() if k in ds_admin_maps.data_vars}
        self.set_grid(ds_admin_maps.rename(rmdict))

    def setup_roads(
        self,
        roads_fn,
        highway_list,
        country_list,
        non_highway_list=None,
        country_fn=None,
    ):
        """Setup roads statistics needed for emission modelling.

        Adds model layers:

        * **km_highway_country** map: emission data with for each grid cell the total km of highway for the country the grid cell is in [km highway/country]
        * **km_other_country** map: emission data with for each grid cell the total km of non highway roads for the country the grid cell is in [km other road/country]
        * **km_highway_cell** map: emission data containing highway length per cell [km highway/cell]
        * **km_other_cell** map:
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
            List of non highway roads in the type variable of roads_fn. If not provided takes every roads except the ones in highway_list.
        country_list: str or list of str, optional.
            List of countries for the model area in country_fn and optionnally in country_code variable of roads_fn.
        country_fn: str, optional.
            Name of country boundaries data source in data_sources.yml file.

            * Required variables: ['country_code']

        """
        if roads_fn is None:
            return
        if roads_fn not in self.data_catalog:
            self.logger.warning(f"Invalid source '{roads_fn}', skipping setup_roads.")
            return
        # Convert string to lists
        if not isinstance(highway_list, list):
            highway_list = [highway_list]
        if not isinstance(country_list, list):
            country_list = [country_list]

        # Mask the road data with countries of interest
        gdf_country = self.data_catalog.get_geodataframe(country_fn, dst_crs=self.crs)
        gdf_country = gdf_country.astype({"country_code": "str"})
        gdf_country = gdf_country.iloc[
            np.isin(gdf_country["country_code"], country_list)
        ]
        # bbox = gdf_country.total_bounds

        # Read the roads data and mask with country geom
        gdf_roads = self.data_catalog.get_geodataframe(
            roads_fn, dst_crs=self.crs, geom=gdf_country, variables=["road_type"]
        )
        # Make sure road type is str format
        gdf_roads["road_type"] = gdf_roads["road_type"].astype(str)
        # Get non_highway_list
        if non_highway_list is None:
            road_types = np.unique(gdf_roads["road_type"].values)
            non_highway_list = road_types[~np.isin(road_types, highway_list)].tolist()
        elif not isinstance(non_highway_list, list):
            non_highway_list = [non_highway_list]

        # Feature filter for highway and non-highway
        feature_filter = {
            "hwy": {"road_type": highway_list},
            "nnhwy": {"road_type": non_highway_list},
        }

        ### Country statistics ###
        # Loop over feature_filter
        for name, colfilter in feature_filter.items():
            self.logger.info(
                f"Computing {name} roads statistics per country of interest"
            )
            # Filter gdf
            colname = [k for k in colfilter.keys()][0]
            subset_roads = gdf_roads.iloc[
                np.isin(gdf_roads[colname], colfilter.get(colname))
            ]
            gdf_country = roads.zonal_stats(
                gdf=subset_roads,
                zones=gdf_country,
                variables=["length"],
                stats=["sum"],
                method="sjoin",
            )
            gdf_country = gdf_country.rename(
                columns={"length_sum": f"{name}_length_sum"}
            )
            # Convert from m to km
            gdf_country[f"{name}_length_sum"] = gdf_country[f"{name}_length_sum"] / 1000

            # Rasterize statistics
            da_emi = emissions.emission_vector(
                gdf=gdf_country,
                ds_like=self.grid,
                col_name=f"{name}_length_sum",
                method="value",
                mask_name="mask",
                logger=self.logger,
            )
            self.set_grid(da_emi.rename(f"{name}_length_sum_country"))

        ### Road statistics per segment ###
        # Filter road gdf with model mask
        mask = self.grid["mask"].astype("int32").raster.vectorize()
        gdf_roads = self.data_catalog.get_geodataframe(
            roads_fn, dst_crs=self.crs, geom=mask, variables=["road_type"]
        )
        # Make sure road type is str format
        gdf_roads["road_type"] = gdf_roads["road_type"].astype(str)

        # Loop over feature_filter
        for name, colfilter in feature_filter.items():
            self.logger.info(f"Computing {name} roads statistics per segment")
            # Filter gdf
            colname = [k for k in colfilter.keys()][0]
            subset_roads = gdf_roads.iloc[
                np.isin(gdf_roads[colname], colfilter.get(colname))
            ]
            ds_roads = roads.zonal_stats_grid(
                gdf=subset_roads,
                ds_like=self.grid,
                variables=["length"],
                stats=["sum"],
                mask_name="mask",
                method="overlay",
            )
            # Convert from m to km and rename
            ds_roads[f"{name}_length"] = ds_roads["length_sum"] / 1000
            ds_roads[f"{name}_length"].attrs.update(_FillValue=0)
            self.set_grid(ds_roads)

    def setup_hydrology_forcing(
        self,
        hydro_forcing_fn: str,
        starttime: str,
        endtime: str,
        timestepsecs: int,
        include_transport: bool = False,
        **kwargs,
    ):
        """Setup Demission hydrological fluxes.

        Adds model layers:

        * **B1_timestamp.inc** config: timestamp at the beginning of the simulation.
        * **B2_outputtimes.inc** config: start and end timestamp for the delwaq outputs.
        * **B2_sysclock.inc** config: model timestep.
        * **B2_timers.inc** config: timers info (start, end, timstep...).

        In EM mode, adds:

        * **hydrology.bin** dynmap: fluxes for EM (Rainfall RunoffPav  RunoffUnp Infiltr TotalFlow) [mm]
        * **B7_hydrology.inc** config: names of fluxes included in hydrology.bin

        Parameters
        ----------
        hydro_forcing_fn : {'None', 'name in local_sources.yml'}
            Either None, or name in a local or global data_sources.yml file.

            * Required variables for EM: ['time', 'precip', 'runPav', 'runUnp', 'infilt', 'exfilt*', 'q_land', 'q_ss']

        startime : str
            Timestamp of the start of Delwaq simulation. Format: YYYY-mm-dd HH:MM:SS
        endtime : str
            Timestamp of the end of Delwaq simulation.Format: YYYY-mm-dd HH:MM:SS
        timestepsecs : int
            Model timestep in seconds.
        include_transport : bool, optional
            If False (default), only use the vertical fluxes for emission [precip, runPav, runUnp, infilt, totflw].
            If True, includes additional fluxes for land and subsurface trasnport [precip, runPav, runUnp, infilt, exfilt, q_land, q_ss].
        """
        if hydro_forcing_fn not in self.data_catalog:
            self.logger.warning(
                f"None or Invalid source '{hydro_forcing_fn}', skipping setup_hydrology_forcing."
            )
            return
        self.logger.info(
            f"Setting dynamic data from hydrology source {hydro_forcing_fn}."
        )
        # read dynamic data
        # Normally region extent is by default exactly the same as hydromodel
        ds_in = self.data_catalog.get_rasterdataset(
            hydro_forcing_fn,
            geom=self.region,
            time_tuple=(starttime, endtime),
        )

        # Select variables based on model type
        vars = ["precip", "runPav", "runUnp", "infilt"]
        if include_transport:
            vars.extend(["q_land", "q_ss"])
            ex_vars = [v for v in ds_in.data_vars if v.startswith("exfilt")]
            vars.extend(ex_vars)
        ds = ds_in[vars]

        # Update model timestepsecs attribute
        self.timestepsecs = timestepsecs

        # Unit conversion (from mm to m3/s)
        for dvar in ds.data_vars.keys():
            if ds[dvar].attrs.get("unit") == "mm":
                attrs = ds[dvar].attrs.copy()
                surface = emissions.gridarea(ds)
                ds[dvar] = ds[dvar] * surface / (1000 * self.timestepsecs)
                ds[dvar].attrs.update(attrs)  # set original attributes
                ds[dvar].attrs.update(unit="m3/s")

        # Sum up exfiltrattion or add totflw depending on include_transport
        if include_transport:
            # Add exfiltration (can be split into several variables)
            # Check if the flux is split into several variables
            ex_vars = [v for v in ds.data_vars if v.startswith("exfilt")]
            if len(ex_vars) > 1:  # need to sum
                attrs = ds[ex_vars[0]].attrs.copy()
                nodata = ds[ex_vars[0]].raster.nodata
                # Cover exfilt with zeros (some negative zeros in wflow outputs?)
                ds["exfilt"] = ds[ex_vars].fillna(0).to_array().sum("variable")
                ds["exfilt"] = ds["exfilt"].where(ds["exfilt"] > 0.0, 0.0)
                ds["exfilt"].attrs.update(attrs)
                ds["exfilt"].raster.set_nodata(nodata)
                ds = ds.drop_vars(ex_vars)
            elif ex_vars[0] != "exfilt":
                ds = ds.rename({ex_vars[0]: "exfilt"})
        else:
            # Add total flow
            ds["totflw"] = ds["precip"].copy()

        # align forcing file with hydromaps
        # as hydro forcing comes from hydro model, it should be aligned with hydromaps
        if not ds.raster.identical_grid(self.hydromaps):
            self.logger.warning(
                "hydro_forcing_fn and model grid are not identical. Reprojecting."
            )
            ds = ds.raster.reproject_like(self.hydromaps)

        # Add _FillValue to the data attributes
        for dvar in ds.data_vars.keys():
            nodata = ds[dvar].raster.nodata
            if nodata is not None:
                ds[dvar].attrs.update(_FillValue=nodata)
            else:
                ds[dvar].attrs.update(_FillValue=-9999.0)

        # Add mask
        ds.coords["mask"] = xr.DataArray(
            dims=ds.raster.dims,
            coords=ds.raster.coords,
            data=self.hydromaps["modelmap"].values,
            attrs=dict(_FillValue=self.hydromaps["mask"].raster.nodata),
        )

        # Add to forcing
        self.set_forcing(ds)

        # Add time info to config
        ST = datetime.strptime(starttime, "%Y-%m-%d %H:%M:%S")
        ET = datetime.strptime(endtime, "%Y-%m-%d %H:%M:%S")
        # B1_timestamp
        lines_ini = {
            "l1": f"'T0: {ST.strftime('%Y.%m.%d %H:%M:%S')}  (scu=       1s)'",
        }
        for option in lines_ini:
            self.set_config("B1_timestamp", option, lines_ini[option])

        # B2_outputtimes
        STstr = ST.strftime("%Y/%m/%d-%H:%M:%S")
        ETstr = ET.strftime("%Y/%m/%d-%H:%M:%S")
        timestep = pd.Timedelta(ds.time.values[1] - ds.time.values[0])
        hours = int(timestep.seconds / 3600)
        minutes = int(timestep.seconds / 60)
        seconds = int(timestep.seconds - minutes * 60)
        minutes -= hours * 60
        timestepstring = "%03d%02d%02d%02d" % (timestep.days, hours, minutes, seconds)
        lines_ini = {
            "l1": f"  {STstr}  {ETstr}  {timestepstring} ; mon/bal",
            "l2": f"  {STstr}  {ETstr}  {timestepstring} ; map",
            "l3": f"  {STstr}  {ETstr}  {timestepstring} ; his",
        }
        for option in lines_ini:
            self.set_config("B2_outputtimes", option, lines_ini[option])

        # B2_sysclock
        timestepsec = timestep.days * 86400 + timestep.seconds
        lines_ini = {
            "l1": f"  1 'DDHHMMSS' 'DDHHMMSS'  ; system clock",
        }
        for option in lines_ini:
            self.set_config("B2_sysclock", option, lines_ini[option])

        # B2_timers
        lines_ini = {
            "l1": f"  {STstr} ; start time",
            "l2": f"  {ETstr} ; stop time",
            "l3": "  0 ; timestep constant",
            "l4": "; dddhhmmss format for timestep",
            "l5": f"{timestepstring} ; timestep",
        }
        for option in lines_ini:
            self.set_config("B2_timers", option, lines_ini[option])

        # B2_timers_only
        lines_ini = {
            "l1": f"  {STstr} ; start time",
            "l2": f"  {ETstr} ; stop time",
            "l3": "  0 ; timestep constant",
        }
        for option in lines_ini:
            self.set_config("B2_timers_only", option, lines_ini[option])

        # Add var info to config
        lines_ini = {
            "l1": "SEG_FUNCTIONS",
            "l2": "Rainfall RunoffPav RunoffUnp Infiltr Exfiltr Overland Subsurface ",
        }
        for option in lines_ini:
            self.set_config("B7_hydrology", option, lines_ini[option])

    # I/O
    def read(self):
        """Method to read the complete model schematization and configuration from file."""
        # self.read_config()
        self.read_geoms()
        self.read_hydromaps()
        # self.read_pointer()
        self.read_grid()
        # self.read_fewsadapter()
        # self.read_forcing()
        self.logger.info("Model read")

    def write(self):
        """Method to write the complete model schematization and configuration to file."""
        self.logger.info(f"Write model data to {self.root}")
        # if in r, r+ mode, only write updated components
        self.write_grid()
        self.write_geoms()
        self.write_config()
        self.write_hydromaps()
        self.write_geometry()
        self.write_forcing()

    def read_geometry(self):
        """Read Delwaq EM geometry file"""
        raise NotImplementedError()

    def write_geometry(self):
        """Write geometry at <root/staticdata> in ASCII and binary format."""
        if not self._write:
            raise IOError("Model opened in read-only mode")
        if self._geometry is not None:
            self.logger.info("Writting geometry file in root/staticdata")
            fname = join(self.root, "config", "B7_geometry")

            # Write ASCII file
            exfile = open(fname + ".inc", "w")
            print(";Geometry of the EM compartment", file=exfile)
            print("PARAMETERS TotArea fPaved fUnpaved fOpenWater ALL", file=exfile)
            print("DATA", file=exfile)
            np.savetxt(exfile, self._geometry.values, fmt="%10.4f")
            exfile.close()

            # Write binary file
            # Flatten the geometry data and repeat them for each compartment
            geometry_data = np.tile(self._geometry.values.flatten(), 1)
            artow = np.array(geometry_data, dtype=np.float32).copy()
            # Define dummy time
            timear = np.array(0, dtype=np.int32)
            # Open and write the data
            fp = open(fname + ".bin", "wb")
            tstr = timear.tobytes() + artow.tobytes()
            fp.write(tstr)
            fp.close()

            # Write corresponding def file of the geometry
            fpa = open(fname + "-parameters.inc", "w")
            print("PARAMETERS TotArea fPaved fUnpaved fOpenWater", file=fpa)
            fpa.close()

    def write_forcing(self, write_nc=False):
        """Write forcing at <root/dynamicdata> in binary format and NetCDF (if write_nc is True)."""
        if not self._write:
            raise IOError("Model opened in read-only mode")
        if not self.forcing:
            self.logger.warning(
                "Warning: no forcing available, skipping write_forcing."
            )
            return

        # Go from dictionnary to xr.DataSet
        ds_out = xr.Dataset()
        for name, da in self.forcing.items():
            ds_out[name] = da

        # To avoid appending data to existing file, first delete all the .dat files
        dynDir = join(self.root, "dynamicdata")
        if os.path.exists(dynDir):
            filelist = os.listdir(dynDir)
            for f in filelist:
                os.remove(os.path.join(dynDir, f))

        # Filter data with mask
        for dvar in ds_out.data_vars:
            # nodata = ds_out[dvar].raster.nodata
            # Change the novalue outside of mask for julia compatibilty
            ds_out[dvar] = xr.where(ds_out["mask"], ds_out[dvar], -9999.0)
            ds_out[dvar].attrs.update(_FillValue=-9999.0)

        self.logger.info("Writing dynamicmap files.")
        # Netcdf format
        if write_nc:
            fname = join(self.root, "dynamicdata", "dynamicmaps.nc")
            ds_out = ds_out.drop_vars(["mask", "spatial_ref"], errors="ignore")
            ds_out.to_netcdf(path=fname)

        # Binary format
        timesteps = np.arange(0, len(ds_out.time.values))
        timestepstamp = np.arange(
            0, (len(ds_out.time.values) + 1) * int(self.timestepsecs), self.timestepsecs
        )
        for i in timesteps:
            self.logger.info(
                f"Writting dynamic data for timestep {i+1}/{timesteps[-1]+1}"
            )
            fname = join(self.root, "dynamicdata", "hydrology.bin")
            datablock = []
            dvars = ["precip", "runPav", "runUnp", "infilt"]
            if "totflw" in ds_out.data_vars:
                dvars.extend(["totflw"])
            else:
                dvars.extend(["exfilt", "q_land", "q_ss"])
            for dvar in dvars:
                nodata = ds_out[dvar].raster.nodata
                data = ds_out[dvar].isel(time=i).values.flatten()
                data = data[data != nodata]
                datablock = np.append(datablock, data)
            self.dw_WriteSegmentOrExchangeData(
                timestepstamp[i], fname, datablock, 1, WriteAscii=False
            )
            # sediment
            if "B7_sediment" in self.config:
                sedname = join(self.root, "dynamicdata", "sediment.dat")
                sed_vars = self.get_config("B7_sediment.l2").split(" ")
                sedblock = []
                for dvar in sed_vars:
                    nodata = ds_out[dvar].raster.nodata
                    data = (
                        ds_out[dvar]
                        .isel(time=i)
                        .transpose("comp", ...)
                        .values.flatten()
                    )
                    data = data[data != nodata]
                    sedblock = np.append(sedblock, data)
                self.dw_WriteSegmentOrExchangeData(
                    timestepstamp[i], sedname, sedblock, 1, WriteAscii=False
                )
            # climate
            if "B7_climate" in self.config:
                climname = join(self.root, "dynamicdata", "climate.dat")
                clim_vars = self.get_config("B7_climate.l2").split(" ")
                climblock = []
                for dvar in clim_vars:
                    nodata = ds_out[dvar].raster.nodata
                    data = (
                        ds_out[dvar]
                        .isel(time=i)
                        .transpose("comp", ...)
                        .values.flatten()
                    )
                    data = data[data != nodata]
                    climblock = np.append(climblock, data)
                self.dw_WriteSegmentOrExchangeData(
                    timestepstamp[i], climname, climblock, 1, WriteAscii=False
                )

    @property
    def nrofseg(self):
        """Fast accessor to nrofseg property of pointer"""
        if "nrofseg" in self.pointer:
            nseg = self.pointer["nrofseg"]
        else:
            # from config
            nseg = self.get_config("B3_nrofseg.l1", "0 ; nr of segments")
            nseg = int(nseg.split(";")[0])
            self.set_pointer(nseg, "nrofseg")
        return nseg

    @property
    def nrofexch(self):
        """Fast accessor to nrofexch property of pointer"""
        if "nrofexch" in self.pointer:
            nexch = self.pointer["nrofexch"]
        else:
            nexch = self.nrofseg * 5
            self.set_pointer(nexch, "nrofexch")
        return nexch

    @property
    def nrofcomp(self):
        """Fast accessor to nrofcomp property of pointer"""
        if "nrofcomp" in self.pointer:
            ncomp = self.pointer["nrofcomp"]
        else:
            ncomp = 1
            self.set_pointer(ncomp, "nrofexch")
        return ncomp

    @property
    def compartments(self):
        return ["EM"]
