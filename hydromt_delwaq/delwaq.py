"""Implement delwaq model class."""

import glob
import logging
import os
import struct
from datetime import datetime
from os.path import isfile, join
from pathlib import Path
from typing import Dict, List, Optional, Union

import geopandas as gpd
import numpy as np
import pandas as pd
import pyproj
import xarray as xr
import xugrid as xu
from hydromt import io, workflows
from hydromt.models.model_grid import GridModel
from tqdm import tqdm

from . import DATADIR
from .workflows import emissions, forcing, hydrology, segments

__all__ = ["DelwaqModel"]

logger = logging.getLogger(__name__)

# specify pcraster map types for non scalar (float) data types
PCR_VS_MAP = {
    "modelmap": "bool",
    "basins": "nominal",
    "ldd": "ldd",
    "ptid": "ordinal",
}


class DelwaqModel(GridModel):
    """Delwaq model class."""

    _NAME = "delwaq"
    _CLI_ARGS = {"region": "setup_basemaps"}
    _CONF = "delwaq.inp"
    _DATADIR = DATADIR
    _CF = dict()  # configreader kwargs
    _GEOMS = {}
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
        logger=logger,
    ):
        super().__init__(
            root=root,
            mode=mode,
            config_fn=config_fn,
            data_libs=data_libs,
            logger=logger,
        )

        # delwaq specific
        self.hydromodel_name = hydromodel_name
        self.hydromodel_root = hydromodel_root

        self._hydromaps = None  # extract of hydromodel grid
        self._pointer = (
            None  # dictionnary of pointer values and related model attributes
        )
        self._geometry = None
        self._fewsadapter = None

        self.timestepsecs = 86400

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
        * **streamorder** map: Strahler stream order map. [-]
        * **slope** map: slope map [-]
        * **river** map: river mask map [-]
        * **surface** map: surface map [m2]
        * **length** map: length map [m]
        * **width** map: width map [m]
        * **pointer** poi: delwaq pointer between segments

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
            is '{surface_water}>{boundary_name}' for a flux from athe surface water to a
            boundary, ex 'sfw>bd'. By default ['sfw>sfw', 'bd>sfw'] for runoff and
            inwater.
            Names in the fluxes list should match name in the hydrology_fn source in
            setup_hydrology_forcing.
        maps: list of str
            List of variables from hydromodel to add to grid.
            By default ['rivmsk', 'lndslp', 'strord', 'N'].
        """
        # Initialise hydromodel from region
        kind, region = workflows.parse_region(region)  # , logger=self.logger)
        if kind != "model":
            raise ValueError("Delwaq model can only be built from 'model' region.")
        hydromodel = region["mod"]
        if hydromodel._NAME != "wflow":
            raise NotImplementedError(
                "Delwaq build function only implemented for wflow base model."
            )
        else:
            self.hydromodel_name = hydromodel._NAME
            self.hydromodel_root = hydromodel.root

        self.logger.info("Preparing WQ basemaps from hydromodel.")

        ### Select and build hydromaps from model ###
        # Initialise hydromaps
        ds_hydro = segments.hydromaps(hydromodel)
        # Add mask
        if mask == "rivers":
            da_mask = ds_hydro["rivmsk"]
        else:
            da_mask = ds_hydro["basmsk"]
        ds_hydro = ds_hydro.drop_vars(["rivmsk", "basmsk"])
        ds_hydro.coords["mask"] = da_mask
        ds_hydro["modelmap"] = da_mask.astype(np.int32)
        ds_hydro["modelmap"].raster.set_nodata(0)
        # Add to hydromaps
        rmdict = {k: v for k, v in self._MAPS.items() if k in ds_hydro.data_vars}
        self.set_hydromaps(ds_hydro.rename(rmdict))

        # Build segment ID and add to hydromaps
        # Prepare delwaq pointer file for WAQ simulation (not needed with EM)
        nrofseg, da_ptid, da_ptiddown, pointer, bd_id, bd_type = segments.pointer(
            self.hydromaps,
            build_pointer=True,
            surface_water=surface_water,
            boundaries=boundaries,
            fluxes=fluxes,
        )
        self.set_hydromaps(da_ptid, name="ptid")
        self.set_hydromaps(da_ptiddown, name="ptiddown")
        # Initialise pointer object with schematisation attributes
        self.set_pointer(pointer, name="pointer")
        self.set_pointer(nrofseg, name="nrofseg")
        self.set_pointer(surface_water, name="surface_water")
        self.set_pointer(bd_type, name="boundaries")
        self.set_pointer(fluxes, name="fluxes")

        ### Initialise grid river and slope ###
        ds_stat = segments.maps_from_hydromodel(
            hydromodel, maps=maps, logger=self.logger
        )
        ds_stat["mask"] = da_mask
        ds_stat = ds_stat.set_coords("mask")
        rmdict = {k: v for k, v in self._MAPS.items() if k in ds_stat.data_vars}
        self.set_grid(ds_stat.rename(rmdict))

        ### Add geometry ###
        ds_geom = segments.geometrymaps(
            hydromodel,
        )
        self.set_grid(ds_geom)

        ### Config ###
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
            "l7": f"     {nrofseg}*01 ; {surface_water}",
            "l8": " 0    ; no time dependent attributes",
        }
        for option in lines_ini:
            self.set_config("B3_attributes", option, lines_ini[option])
        # B4_nrofexch
        lines_ini = {
            "l1": f"{self.nrofexch} 0 0 ; x, y, z direction",
        }
        for option in lines_ini:
            self.set_config("B4_nrofexch", option, lines_ini[option])
        # B5_boundlist
        lines_ini = {
            "l1": ";'NodeID' 'Number' 'Type'",
        }
        for i in range(len(bd_id)):
            lstr = "l" + str(i + 2)
            lines_ini.update(
                {lstr: f"'BD_{int(bd_id[i])}' '{int(bd_id[i])}' '{bd_type[i]}'"}
            )
        for option in lines_ini:
            self.set_config("B5_boundlist", option, lines_ini[option])
        # B7_fluxes
        lines_ini = {
            "l1": "SEG_FUNCTIONS",
            "l2": " ".join(fluxes),
        }
        for option in lines_ini:
            self.set_config("B7_fluxes", option, lines_ini[option])

    def setup_monitoring(self, mon_points: str = None, mon_areas: str = None, **kwargs):
        """Prepare Delwaq monitoring points and areas options.

        Adds model layers:

        * **monitoring_points** map: monitoring points segments
        * **monitoring_areas** map: monitoring areas ID
        * **B2_nrofmon.inc** config: number of monitoring points and areas

        Parameters
        ----------
        mon_points : {'segments', 'data_source', 'path to station location'}
            Either source from DataCatalog, path to a station location dataset
            or if segments, all segments are monitored.
        mon_areas : str {'subcatch', 'riverland'}
            Either subcatch from hydromodel or 'riverland' to split river and land
            cells.
        """
        self.logger.info("Setting monitoring points and areas")
        monpoints = None
        mv = -999
        # Read monitoring points source
        if mon_points is not None:
            if mon_points == "segments":
                monpoints = self.hydromaps["ptid"]
            else:
                kwargs = {}
                if isfile(mon_points):
                    kwargs.update(crs=self.crs)
                gdf = self.data_catalog.get_geodataframe(
                    mon_points, geom=self.basins, assert_gtype="Point", **kwargs
                )
                gdf = gdf.to_crs(self.crs)
                if gdf.index.size == 0:
                    self.logger.warning(
                        f"No {mon_points} gauge locations found within domain"
                    )
                else:
                    gdf.index.name = "index"
                    monpoints = self.grid.raster.rasterize(
                        gdf, col_name="index", nodata=mv
                    )
            self.logger.info(f"Gauges locations read from {mon_points}")
            # Add to grid
            self.set_grid(monpoints.rename("monpoints"))
            # Number of monitoring points
            points = monpoints.values.flatten()
            points = points[points != mv]
            nb_points = len(points)
            # Add to geoms if mon_points is not segments
            if mon_points != "segments":
                self.set_geoms(gdf, name="monpoints")
        else:
            self.logger.info("No monitoring points set in the config file, skipping")
            nb_points = 0

        # Monitoring areas domain
        monareas = None
        if mon_areas is not None:
            if mon_areas == "subcatch":  # subcatch
                basins = self.hydromaps["basins"].copy()
                # Number or monitoring areas
                areas = basins.values.flatten()
                areas = areas[areas != mv]
                nb_areas = len(np.unique(areas))
                monareas = basins
            else:  # riverland
                # seperate areas for land cells (1) and river cells (2)
                lr_areas = xr.where(
                    self.hydromaps["river"],
                    2,
                    xr.where(self.hydromaps["basins"], 1, mv),
                ).astype(np.int32)
                # Number or monitoring areas
                areas = lr_areas.values.flatten()
                areas = areas[areas != mv]
                nb_areas = len(np.unique(areas))
                monareas = lr_areas

            # Add to grid
            self.set_grid(monareas.rename("monareas"))
            self.grid["monareas"].attrs.update(_FillValue=mv)
            self.grid["monareas"].attrs["mon_areas"] = mon_areas

            # Add to geoms
            gdf_areas = self.grid["monareas"].raster.vectorize()
            self.set_geoms(gdf_areas, name="monareas")
        else:
            self.logger.info("No monitoring areas set in the config file, skipping")
            nb_areas = 0

        # Config
        lines_ini = {
            "l1": f";{nb_points} monitoring points",
            "l2": f";{nb_areas} monitoring areas",
            "l3": f"{nb_points+nb_areas} ; nr of monitoring points/areas",
        }
        for option in lines_ini:
            self.set_config("B2_nrofmon", option, lines_ini[option])

    def setup_flooding(
        self,
        hydro_fn: str,
        from_volume: bool = True,
        bankfull_rp: int = 2,
        **kwargs,
    ):
        """
        Derive bankfull volume from either discharge or volume timeseries.

        Addionnally, if floodplain processes need to be modelled, the bankfull volume
        to characterise flooding thresholds can be derived from hydrology timeseries
        data.
        The bankfull volume is defined as the volume corresponding to a return period of
        2 years (default).
        xclim.indices.stats.frequency_analysis is used to get bankfull volume based on
        annual maxima method.

        If volume timeseries are not available, the bankfull volume can be derived from
        discharge timeseries data using the flag from_volume=False.

        Parameters
        ----------
        hydro_fn : str
            Path to the discharge forcing file.

            * Required variable: ["volume"] or ["discharge"]

        from_volume : bool, optional
            If True, bankfull volume is derived from volume timeseries data.
        **kwargs
            Keywords arguments for either bankfull_volume or
            bankfull_volume_from_discharge.
        """
        # Read data
        # Normally region extent is by default exactly the same as hydromodel
        if from_volume:
            da_v = self.data_catalog.get_rasterdataset(
                hydro_fn, geom=self.region, variables=["volume"]
            )

            # Derive bankfull discharge / depth / volume
            self.logger.info("Deriving flooding characteristics from volume.")
            ds_bankfull = hydrology.bankfull_volume(
                da_v=da_v,
                bankfull_rp=bankfull_rp,
                **kwargs,
            )
        else:
            da_q = self.data_catalog.get_rasterdataset(
                hydro_fn, geom=self.region, variables=["discharge"]
            )

            # Derive bankfull discharge / depth / volume
            self.logger.info("Deriving flooding characteristics from discharge.")
            ds_bankfull = hydrology.bankfull_volume_from_discharge(
                da_q=da_q,
                ds_model=self.grid,
                bankfull_rp=bankfull_rp,
                **kwargs,
            )

        # Add to grid
        self.set_grid(ds_bankfull)

    def setup_hydrology_forcing(
        self,
        hydro_forcing_fn: str,
        starttime: str,
        endtime: str,
        timestepsecs: int = 86400,
        add_volume_offset: bool = True,
        min_volume: float = 0.1,
        override: List = [],
    ):
        """Prepare Delwaq hydrological fluxes.

        As the fluxes order should precisely macth the pointer defined in
        setup_basemaps, the variables names in ``hydro_forcing_fn`` should match names
        defined in the ``fluxes`` argument of setup_basemaps. These names should also
        have been saved in the file config/B7_fluxes.inc.

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

        * **B1_timestamp.inc** config: timestamp at the beginning of the simulation.
        * **B2_outputtimes.inc** config: start and end timestamp for the delwaq outputs.
        * **B2_sysclock.inc** config: model timestep.
        * **B2_timers.inc** config: timers info (start, end, timstep...).
        * **flow.dat** dynmap: water fluxes [m3/s]
        * **volume.dat** dynmap: water volumes [m3]
        * **area.dat** dynmap: water cross sections [m2]
        * **velocity.dat** dynmap: flow velocity [m/s]

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
        if hydro_forcing_fn not in self.data_catalog:
            self.logger.warning(
                f"None or Invalid source {hydro_forcing_fn} ",
                "skipping setup_hydrology_forcing.",
            )
            return
        self.logger.info(
            f"Setting dynamic data from hydrology source {hydro_forcing_fn}."
        )
        # read dynamic data
        # Normally region extent is by default exactly the same as hydromodel
        ds = self.data_catalog.get_rasterdataset(hydro_forcing_fn, geom=self.region)

        # Prepare hydrology forcing
        ds_out = forcing.hydrology_forcing(
            ds=ds,
            ds_model=self.hydromaps,
            timestepsecs=timestepsecs,
            time_tuple=(starttime, endtime),
            fluxes=self.fluxes,
            surface_water=self.surface_water,
            add_volume_offset=add_volume_offset,
            min_volume=min_volume,
            override=override,
            logger=self.logger,
        )
        self.set_forcing(ds_out)

        # Update model timestepsecs attribute
        self.timestepsecs = timestepsecs

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
        timestep = pd.Timedelta(ds_out.time.values[1] - ds_out.time.values[0])
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
            "l1": f"{timestepsec:7d} 'DDHHMMSS' 'DDHHMMSS'  ; system clock",
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

    def setup_sediment_forcing(
        self,
        sediment_fn: Union[str, xr.Dataset],
        starttime: str,
        endtime: str,
        timestepsecs: int = 86400,
        particle_class: List[str] = ["IM1", "IM2", "IM3", "IM4", "IM5"],
    ):
        """Prepare Delwaq sediment fluxes.

        Adds model layers:

        * **sediment.bin** dynmap: sediment particles from land erosion (fErodIM*)
          [g/timestep]
        * **B7_sediment.inc** config: names of sediment fluxes included in sediment.bin

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
        self.logger.info(f"Setting dynamic data from sediment source {sediment_fn}.")
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
            ds_model=self.hydromaps,
            timestepsecs=timestepsecs,
            logger=self.logger,
        )
        # Add to forcing
        self.set_forcing(ds_out)

        # Update config
        lines_ini = {
            "l1": "SEG_FUNCTIONS",
            "l2": " ".join(str(x) for x in sed_vars),
        }
        for option in lines_ini:
            self.set_config("B7_sediment", option, lines_ini[option])

    def setup_climate_forcing(
        self,
        climate_fn: Union[str, xr.Dataset],
        starttime: str,
        endtime: str,
        timestepsecs: int,
        climate_vars: list = ["temp", "temp_dew", "ssr", "wind10_u", "wind10_v", "tcc"],
        temp_correction: bool = False,
        dem_forcing_fn: str = None,
    ):
        """Prepare Delwaq climate fluxes.

        Adds model layers:

        * **climate.bin** dynmap: climate fuxes for climate_vars
        * **B7_climate.inc** config: names of climate fluxes included in climate.bin

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
        self.logger.info(f"Setting dynamic data from climate source {climate_fn}.")

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
            ds_model=self.hydromaps,
            timestepsecs=timestepsecs,
            temp_correction=temp_correction,
            dem_forcing=dem_forcing,
        )

        # Rename and add to forcing
        rmdict = {k: v for k, v in self._FORCING.items() if k in ds_out.data_vars}
        ds_out = ds_out.rename(rmdict)
        self.set_forcing(ds_out)

        # Update config
        lines_ini = {
            "l1": "SEG_FUNCTIONS",
            "l2": " ".join(str(x) for x in ds_out.data_vars),
        }
        for option in lines_ini:
            self.set_config("B7_climate", option, lines_ini[option])

    def setup_emission_raster(
        self,
        emission_fn: Union[str, xr.DataArray],
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
        ds_emi = ds_emi.to_dataset(name=emission_fn)
        self.set_grid(ds_emi)

    def setup_emission_vector(
        self,
        emission_fn: Union[str, xr.DataArray],
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
        self.logger.info(f"Preparing '{emission_fn}' map.")
        gdf_org = self.data_catalog.get_geodataframe(
            emission_fn, geom=self.basins, dst_crs=self.crs
        )
        if gdf_org.empty:
            self.logger.warning(
                f"No shapes of {emission_fn} found within region, "
                "setting to default value."
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
        ds_emi = ds_emi.to_dataset(name=emission_fn)
        self.set_grid(ds_emi)

    def setup_emission_mapping(
        self,
        region_fn: Union[str, gpd.GeoDataFrame],
        mapping_fn: Union[str, Path] = None,
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
        self.logger.info(f"Preparing region_fn related parameter maps for {region_fn}.")
        if mapping_fn is None:
            self.logger.warning("Using default mapping file.")
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

    # I/O
    def read(self):
        """Read the complete model schematization and configuration from file."""
        self.read_config()
        self.read_geoms()
        self.read_hydromaps()
        # self.read_pointer()
        self.read_grid()
        # self.read_fewsadapter()
        # self.read_forcing()
        self.logger.info("Model read")

    def write(self):
        """Write the complete model schematization and configuration to file."""
        self.logger.info(f"Write model data to {self.root}")
        # if in r, r+ mode, only write updated components
        self.write_grid()
        self.write_geoms()
        self.write_config()
        self.write_hydromaps()
        self.write_pointer()
        self.write_forcing()

    def read_grid(self, **kwargs):
        """Read grid at <root/staticdata> and parse to xarray."""
        fn = join(self.root, "staticdata", "staticdata.nc")
        super().read_grid(fn, **kwargs)

    def write_grid(self):
        """Write grid at <root/staticdata> in NetCDF and binary format."""
        if len(self.grid) == 0:
            self.logger.debug("No grid data found, skip writing.")
            return

        self._assert_write_mode()
        ds_out = self.grid

        # Filter data with mask
        for dvar in ds_out.data_vars:
            ds_out[dvar] = ds_out[dvar].raster.mask(mask=ds_out["mask"])

        self.logger.info("Writing staticmap files.")
        # Netcdf format
        fname = join(self.root, "staticdata", "staticdata.nc")
        # Update attributes for gdal compliance
        # ds_out = ds_out.raster.gdal_compliant(rename_dims=False)
        ds_out.to_netcdf(path=fname)

        # Binary format
        mask = ds_out["mask"].values.flatten()
        for dvar in ds_out.data_vars:
            if dvar != "monpoints" and dvar != "monareas":
                fname = join(self.root, "staticdata", dvar + ".dat")
                data = ds_out[dvar].values.flatten()
                data = data[mask]
                self.dw_WriteSegmentOrExchangeData(
                    0, fname, data, 1, WriteAscii=False, mode="w"
                )

        # Monitoring files format
        monpoints = None
        monareas = None
        if "monpoints" in ds_out.data_vars:
            monpoints = ds_out["monpoints"]
        if "monareas" in ds_out.data_vars:
            monareas = ds_out["monareas"]
        self.write_monitoring(monpoints, monareas)

    def read_geoms(self):
        """Read and geoms at <root/geoms> and parse to geopandas."""
        super().read_geoms(fn="geoms/*.geojson")

    def write_geoms(self):
        """Write grid at <root/geoms> in model ready format."""
        # to write use self.geoms[var].to_file()
        super().write_geoms(fn="geoms/{name}.geojson", driver="GeoJSON")

    def read_config(
        self,
        skip: List[str] = [
            "B4_pointer",
            "B2_stations",
            "B2_stations-balance",
            "B2_monareas",
        ],
    ):
        """Read config files in ASCII format at <root/config>."""
        # Because of template config can be read also in write mode.
        # Use fonction from hydromt core to read the main config file
        super().read_config()

        # Add the other files in the config folder
        config_fns = glob.glob(join(self.root, "config", "*.inc"))
        for fn in config_fns:
            name = os.path.splitext(os.path.basename(fn))[0]
            if name in skip:
                # Skip pointer file (should be read with read_pointer())
                # Skip monitoring files (should be read with read_monitoring())
                continue
            self.config[name] = dict()
            with open(fn) as f:
                for i, line in enumerate(f):
                    # Remove line breaks
                    self.config[name][f"l{i+1}"] = line.replace("\n", "")

    def write_config(self):
        """Write config files in ASCII format at <root/config>."""
        self._assert_write_mode()
        if self.config:
            self.logger.info("Writing model config to file.")
            for name, lines in self.config.items():
                fn_out = join(self.root, "config", f"{name}.inc")
                exfile = open(fn_out, "w")
                for lnb, line in lines.items():
                    print(line, file=exfile)
                exfile.close()

    def read_hydromaps(self, crs=None, **kwargs):
        """Read hydromaps at <root/hydromodel> and parse to xarray."""
        self._assert_read_mode()
        self._initialize_hydromaps(skip_read=True)
        if "chunks" not in kwargs:
            kwargs.update(chunks={"y": -1, "x": -1})
        # Load grid data in r+ mode to allow overwritting netcdf files
        if self._read and self._write:
            kwargs["load"] = True
        fns = glob.glob(join(self.root, "hydromodel", "*.tif"))
        if len(fns) > 0:
            ds_hydromaps = io.open_mfraster(fns, **kwargs)
        if ds_hydromaps.raster.crs is None and crs is not None:
            ds_hydromaps.raster.set_crs(crs)
        ds_hydromaps.coords["mask"] = ds_hydromaps["modelmap"].astype(bool)
        self.set_hydromaps(ds_hydromaps)

    def write_hydromaps(self):
        """Write hydromaps at <root/hydromodel> in PCRaster maps format."""
        self._assert_write_mode()
        if len(self.hydromaps) == 0:
            self.logger.debug("No grid data found, skip writing.")
            return

        ds_out = self.hydromaps
        self.logger.info("Writing hydromap files.")
        # Convert bool dtype before writting
        for var in ds_out.raster.vars:
            if ds_out[var].dtype == "bool":
                ds_out[var] = ds_out[var].astype(np.int32)
        ds_out.raster.to_mapstack(
            root=join(self.root, "hydromodel"),
            mask=True,
        )

    def read_pointer(self):
        """Read Delwaq pointer file."""
        raise NotImplementedError()

    def write_pointer(self):
        """Write pointer at <root/dynamicdata> in ASCII and binary format."""
        self._assert_write_mode()
        if self._pointer is not None and "pointer" in self._pointer:
            pointer = self.pointer["pointer"]
            self.logger.info("Writting pointer file in root/config")
            fname = join(self.root, "config", "B4_pointer")
            # Write ASCII file
            exfile = open((fname + ".inc"), "w")
            print(";Pointer for WAQ simulation in Surface Water", file=exfile)
            print(";nr of pointers is: ", str(pointer.shape[0]), file=exfile)
            np.savetxt(exfile, pointer, fmt="%10.0f")
            exfile.close()

            # Write binary file
            f = open((fname + ".poi"), "wb")
            for i in range(pointer.shape[0]):
                f.write(struct.pack("4i", *np.int_(pointer[i, :])))
            f.close()

    def read_forcing(self):
        """Read and forcing at <root/?/> and parse to dict of xr.DataArray."""
        self._assert_read_mode()
        self._forcing = dict()
        # raise NotImplementedError()

    def write_forcing(self, write_nc=False):
        """
        Write grid at <root/staticdata> in binary format.

        Can also write a NetCDF copy if write_nc is True.
        """
        self._assert_write_mode()
        if len(self.forcing) == 0:
            self.logger.debug("No forcing data found, skip writing.")
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
            ds_out[dvar] = ds_out[dvar].where(ds_out["mask"], -9999.0)
            ds_out[dvar].attrs.update(_FillValue=-9999.0)

        self.logger.info("Writing dynamicmap files.")
        # Netcdf format
        if write_nc:
            fname = join(self.root, "dynamicdata", "dynamicmaps.nc")
            ds_out = ds_out.drop_vars(["mask", "spatial_ref"], errors="ignore")
            ds_out.to_netcdf(path=fname)

        # Binary format
        # timesteps = np.arange(0, len(ds_out.time.values))
        timestepstamp = np.arange(
            0, (len(ds_out.time.values) + 1) * int(self.timestepsecs), self.timestepsecs
        )

        for i in tqdm(
            np.arange(0, len(ds_out.time.values)), desc="Writing dynamic data"
        ):
            # self.logger.info(
            #    f"Writting dynamic data for timestep {i+1}/{len(timesteps)}"
            # )
            # Flow
            flname = join(self.root, "dynamicdata", "flow.dat")
            flow_vars = self.fluxes
            flowblock = []
            for dvar in flow_vars:
                nodata = ds_out[dvar].raster.nodata
                data = ds_out[dvar].isel(time=i).values.flatten()
                data = data[data != nodata]
                flowblock = np.append(flowblock, data)
            self.dw_WriteSegmentOrExchangeData(
                timestepstamp[i], flname, flowblock, 1, WriteAscii=False
            )
            # volume
            voname = join(self.root, "dynamicdata", "volume.dat")
            vol_vars = [self.surface_water]
            volblock = []
            for dvar in vol_vars:
                nodata = ds_out[dvar].raster.nodata
                data = ds_out[dvar].isel(time=i).values.flatten()
                data = data[data != nodata]
                volblock = np.append(volblock, data)
            self.dw_WriteSegmentOrExchangeData(
                timestepstamp[i], voname, volblock, 1, WriteAscii=False
            )
            # sediment
            if "B7_sediment" in self.config:
                sedname = join(self.root, "dynamicdata", "sediment.dat")
                sed_vars = self.get_config("B7_sediment.l2").split(" ")
                sedblock = []
                for dvar in sed_vars:
                    nodata = ds_out[dvar].raster.nodata
                    data = ds_out[dvar].isel(time=i).values.flatten()
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
                    data = ds_out[dvar].isel(time=i).values.flatten()
                    data = data[data != nodata]
                    climblock = np.append(climblock, data)
                self.dw_WriteSegmentOrExchangeData(
                    timestepstamp[i], climname, climblock, 1, WriteAscii=False
                )

    def read_states(self):
        """Read states at <root/?/> and parse to dict of xr.DataArray."""
        self._assert_read_mode()
        self._states = dict()
        # raise NotImplementedError()

    def write_states(self):
        """Write states at <root/?/> in model ready format."""
        self._assert_write_mode()
        raise NotImplementedError()

    def read_results(self):
        """Read results at <root/?/> and parse to dict of xr.DataArray."""
        self._assert_read_mode()
        self._results = dict()
        # raise NotImplementedError()

    def write_results(self):
        """Write results at <root/?/> in model ready format."""
        self._assert_write_mode()
        raise NotImplementedError()

    ## DELWAQ specific data and methods

    @property
    def basins(self):
        """Derive basins from geoms or hydromaps."""
        if "basins" in self.geoms:
            gdf = self.geoms["basins"]
        elif "basins" in self.hydromaps:
            gdf = self.hydromaps["basins"].raster.vectorize()
            gdf.crs = pyproj.CRS.from_user_input(self.crs)
            self.set_geoms(gdf, name="basins")
        return gdf

    @property
    def hydromaps(self):
        """xarray.dataset representation of all hydrology maps."""
        if self._hydromaps is None:
            self._initialize_hydromaps()
        return self._hydromaps

    def _initialize_hydromaps(self, skip_read=False) -> None:
        """Initialize grid object."""
        if self._hydromaps is None:
            self._hydromaps = xr.Dataset()
            if self._read and not skip_read:
                self.read_hydromaps()

    def set_hydromaps(
        self,
        data: Union[xr.DataArray, xr.Dataset, np.ndarray],
        name: Optional[str] = None,
    ):
        """Add data to hydromaps re-using the set_grid method."""
        self._initialize_hydromaps()
        name_required = isinstance(data, np.ndarray) or (
            isinstance(data, xr.DataArray) and data.name is None
        )
        if name is None and name_required:
            raise ValueError(f"Unable to set {type(data).__name__} data without a name")
        if isinstance(data, np.ndarray):
            if data.shape != self.hydromaps.raster.shape:
                raise ValueError("Shape of data and hydromaps do not match")
            data = xr.DataArray(dims=self.hydromaps.raster.dims, data=data, name=name)
        if isinstance(data, xr.DataArray):
            if name is not None:  # rename
                data.name = name
            data = data.to_dataset()
        elif not isinstance(data, xr.Dataset):
            raise ValueError(f"cannot set data of type {type(data).__name__}")
        # force read in r+ mode
        if len(self._hydromaps) == 0:  # empty hydromaps
            self._hydromaps = data
        else:
            for dvar in data.data_vars:
                if dvar in self.hydromaps:
                    if self._read:
                        self.logger.warning(f"Replacing hydromap: {dvar}")
                self._hydromaps[dvar] = data[dvar]

    @property
    def pointer(self):
        """
        Dictionnary of schematisation attributes of a Delwaq model.

        Contains
        --------
        pointer: np.array
            Model pointer defining exchanges between segments
        surface_water: list of str
            Name of the surface water layer
        boundaries: list of str
            List of model boundaries names
        fluxes: list of str
            List of model fluxes names
        nrofseg: int
            number of segments
        nrofexch: int
            number of exchanges
        """
        if not self._pointer:
            # not implemented yet, fix later
            self._pointer = dict()
            # if self._read:
            #    self.read_pointer
        return self._pointer

    def set_pointer(self, attr, name):
        """Add model attribute property to pointer."""
        # Check that pointer attr is a four column np.array
        if name == "pointer":
            if not isinstance(attr, np.ndarray) and attr.shape[1] == 4:
                self.logger.warning(
                    "pointer values in self.pointer should be a"
                    "np.ndarray with four columns."
                )
                return
        elif np.isin(name, ["csurface_water", "boundaries", "fluxes"]):
            if not isinstance(attr, list):
                self.logger.warning(
                    f"{name} object in self.pointer should be a list of names."
                )
                return
        elif np.isin(name, ["nrofseg", "nrofexch"]):
            if not isinstance(attr, int):
                self.logger.warning(
                    f"{name} object in self.pointer should be an integer."
                )
                return
        if self._pointer is None:
            self._pointer = {name: attr}
        else:
            self._pointer[name] = attr

    @property
    def nrofseg(self):
        """Fast accessor to nrofseg property of pointer."""
        if "nrofseg" in self.pointer:
            nseg = self.pointer["nrofseg"]
        else:
            # from config
            nseg = self.get_config("B3_nrofseg.l1", fallback="0 ; nr of segments")
            nseg = int(nseg.split(";")[0])
            self.set_pointer(nseg, "nrofseg")
        return nseg

    @property
    def nrofexch(self):
        """Fast accessor to nrofexch property of pointer."""
        if "nrofexch" in self.pointer:
            nexch = self.pointer["nrofexch"]
        elif "pointer" in self.pointer:
            nexch = self.pointer["pointer"].shape[0]
            self.set_pointer(nexch, "nrofexch")
        else:
            # from config
            nexch = self.get_config(
                "B4_nrofexch.l1",
                fallback="0 0 0 ; x, y, z direction",
            )
            nexch = int(nexch.split(" ")[0])
            self.set_pointer(nexch, "nrofexch")
        return nexch

    @property
    def surface_water(self):
        """Fast accessor to surface_water name property of pointer."""
        if "surface_water" in self.pointer:
            sfw = self.pointer["surface_water"]
        else:
            # from config
            nl = 7
            sfw = self.get_config(
                f"B3_attributes.l{nl}",
                fallback="     1*01 ; sfw",
            )
            sfw = sfw.split(";")[1][1:]
            self.set_pointer(sfw, "surface_water")
        return sfw

    @property
    def fluxes(self):
        """Fast accessor to fluxes property of pointer."""
        if "fluxes" in self.pointer:
            fl = self.pointer["fluxes"]
        else:
            # from config
            fl = self.get_config(
                "B7_fluxes.l2",
                fallback="sfw>sfw inw>sfw",
            )
            fl = fl.split(" ")
            self.set_pointer(fl, "fluxes")
        return fl

    def dw_WriteSegmentOrExchangeData(
        self,
        ttime,
        fname,
        datablock,
        boundids,
        WriteAscii=True,
        mode="a",
    ):
        """
        Write a timestep to a segment/exchange data file.

        Either appends to an existing file or creates a new one.

        Input:
            - time - timestep number for this timestep
            - fname - File path of the segment/exchange data file</param>
            - datablock - array with data
            - boundids to write more than 1 block
            - WriteAscii - if True to make a copy in an ascii checkfile
            - mode - {"a", "w"} Force the writting mode, append or overwrite existing
              files.

        """
        # Supress potential NaN values to avoid error (replaced by -1.0)
        datablock[np.isnan(datablock)] = -1.0
        # Convert the array to a 32 bit float
        totareas = datablock
        for i in range(boundids - 1):
            totareas = np.vstack((totareas, datablock))

        artow = np.array(totareas, dtype=np.float32).copy()
        timear = np.array(ttime, dtype=np.int32)

        if os.path.isfile(fname) and mode == "a":  # append to existing file
            fp = open(fname, "ab")
            tstr = timear.tobytes() + artow.tobytes()
            fp.write(tstr)
            if WriteAscii:
                fpa = open(fname + ".asc", "a")
                timear.tofile(fpa, format="%d\t", sep=":")
                artow.tofile(fpa, format="%10.8f", sep="\t")
                fpa.write("\n")
        else:
            fp = open(fname, "wb")
            tstr = timear.tobytes() + artow.tobytes()
            fp.write(tstr)
            if WriteAscii:
                fpa = open(fname + ".asc", "w")
                timear.tofile(fpa, format="%d\t", sep=":")
                artow.tofile(fpa, format="%10.8f", sep="\t")
                fpa.write("\n")

        fp.close()
        if WriteAscii:
            fpa.close()

    def write_monitoring(self, monpoints, monareas):
        """
        Write monitoring files and config in ASCII format.

        Input:
            - monpoints - xr.DataArray of monitoring points location
            - monareas - xr.DataArray of monitoring areas location
        """
        ptid = self.hydromaps["ptid"].values.flatten()
        # Monitoring points
        if monpoints is not None:
            mv = monpoints.raster.nodata
            points = monpoints.values.flatten()
            id_points = ptid[points != mv]
            points = points[points != mv]
            nb_points = len(points)
            names_points = np.array(
                ["'Point" + x1 + "_Sfw'" for x1 in points.astype(str)]
            ).reshape(nb_points, 1)
            onecol = np.repeat(1, nb_points).reshape(nb_points, 1)
            balcol = np.repeat("NO_BALANCE", nb_points).reshape(nb_points, 1)
            stations = np.hstack(
                (names_points, balcol, onecol, id_points.reshape(nb_points, 1))
            )
            stations_balance = np.hstack(
                (names_points, onecol, id_points.reshape(nb_points, 1))
            )
            # Write to file
            for name in ["stations", "stations-balance"]:
                fname = join(self.root, "config", "B2_" + name + ".inc")
                exfile = open(fname, "w")
                print(";Written by hydroMT", file=exfile)
                if name == "stations":
                    np.savetxt(exfile, stations, fmt="%.20s")
                else:
                    np.savetxt(exfile, stations_balance, fmt="%.20s")
                exfile.close()
        else:
            fname = join(self.root, "config", "B2_stations.inc")
            exfile = open(fname, "w")
            print(";Written by hydroMT: no monitoring points were set.", file=exfile)
            exfile.close()

        # Monitoring areas
        if monareas is not None:
            mv = monareas.raster.nodata
            areas = monareas.values.flatten()
            id_areas = ptid[areas != mv]
            areas = areas[areas != mv]
            # Write to file
            fname = join(self.root, "config", "B2_monareas.inc")
            exfile = open(fname, "w")
            print(";Written by hydroMT", file=exfile)
            for i in np.unique(areas):
                id_areasi = id_areas[areas == i]
                # Reshape id_areasi as max characters /line in ASCII file is 1000
                # Max allowed number ID of cells has 20 characters -> 50 cells / row
                NTOT = len(id_areasi)
                # Number of complete rows
                NCOMP = int(len(id_areasi) / 50)
                areai_1 = id_areasi[0 : NCOMP * 50].reshape(NCOMP, 50)
                areai_2 = id_areasi[NCOMP * 50 : NTOT]
                areai_2 = areai_2.reshape(1, len(areai_2))
                if monareas.attrs["mon_areas"] == "riverland":
                    if i == 1:
                        print(f"'{'land'}'        {NTOT}", file=exfile)
                    else:  # i = 2
                        print(f"'{'river'}'        {NTOT}", file=exfile)
                else:  # 'subcatch' or 'compartments'
                    print(f"'{i}'        {NTOT}", file=exfile)
                np.savetxt(exfile, areai_1, fmt="%10.20s")
                np.savetxt(exfile, areai_2, fmt="%10.20s")
            exfile.close()
        else:
            fname = join(self.root, "config", "B2_monareas.inc")
            exfile = open(fname, "w")
            print(";Written by hydroMT: no monitoring areas were set.", file=exfile)
            exfile.close()

    def write_waqgeom(self):
        """Write Delwaq netCDF geometry file (config/B3_waqgeom.nc)."""
        # Add waqgeom.nc file to allow Delwaq to save outputs in nc format
        # TODO: Update for several layers
        self.logger.info("Writting waqgeom.nc file")
        ptid = self.hydromaps["ptid"].copy()
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
        np_ldd = self.hydromaps["ldd"].squeeze(drop=True).values
        np_ldd[np_ldd == self.hydromaps["ldd"].raster.nodata] = 0

        # print("Input DataArray: ")
        # print(ptid.dims, ptid)
        ptid = xr.where(ptid == -1, np.nan, ptid)
        if ptid.raster.y_dim != "y":
            ptid = ptid.rename({ptid.raster.y_dim: "y"})
        if ptid.raster.x_dim != "x":
            ptid = ptid.rename({ptid.raster.x_dim: "x"})
        # ptid = ptid.rename({"lat": "y", "lon": "x"})
        da_ptid = xu.UgridDataArray.from_structured(ptid)
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

        # Write the waqgeom.nc file
        fname = join(self.root, "config", "B3_waqgeom.nc")
        uda_waqgeom.to_netcdf(path=fname, mode="w")
