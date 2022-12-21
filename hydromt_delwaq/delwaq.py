"""Implement delwaq model class"""

import os
from os.path import join, isfile, basename
import glob
import numpy as np
import pandas as pd
import xarray as xr
import pyproj
import logging
import struct
from datetime import datetime
import time as t

import hydromt
from hydromt.models.model_api import Model
from hydromt import workflows, flw, io


from hydromt_wflow.wflow import WflowModel

from .workflows import emissions, segments, roads
from . import DATADIR

__all__ = ["DelwaqModel"]

logger = logging.getLogger(__name__)

# specify pcraster map types for non scalar (float) data types
PCR_VS_MAP = {
    "modelmap": "bool",
    "basins": "nominal",
    "ldd": "ldd",
    "ptid": "ordinal",
}


class DelwaqModel(Model):
    """This is the delwaq model class"""

    _NAME = "delwaq"
    _CONF = "delwaq.inp"
    _DATADIR = DATADIR
    _CF = dict()  # configreader kwargs
    _GEOMS = {}
    _MAPS = {
        "flwdir": "ldd",
        "basmsk": "modelmap",
    }
    _FOLDERS = [
        "hydromodel",
        "staticdata",
        "staticgeoms",
        "config",
        "dynamicdata",
        "fews",
    ]

    def __init__(
        self,
        root=None,
        mode="w",
        config_fn=None,
        mtype="EM",
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
            data_libs=data_libs,
            deltares_data=deltares_data,
            logger=logger,
        )

        # delwaq specific
        self.hydromodel_name = hydromodel_name
        self.hydromodel_root = hydromodel_root

        self._hydromaps = xr.Dataset()  # extract of hydromodel staticmaps
        self._pointer = None
        self._geometry = None
        self._fewsadapter = None

        self.timestepsecs = 86400

        self.type = mtype

    def setup_basemaps(
        self,
        region,
        include_soil=False,
    ):
        """Setup the delwaq model schematization using the hydromodel region and
        resolution. 
        
        Maps used and derived from the hydromodel are stored in a specific\ 
        hydromodel attribute. Depending on the global option ``type``, build a\ 
        one-substance D-Emission ("EM") or D-Water Quality ("WQ") case.\
        For a WQ case, the pointer will be created. No arguments are needed for this\
        function (``include_soil`` will be supported later for subsurface transport of\
        substances).
        
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
        include_soil : boolean
            add a soil compartment for transport (default: False)
        """
        # Initialise hydromodel from region
        kind, region = workflows.parse_region(region, logger=self.logger)
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

        self.logger.info(f"Preparing {self.type} basemaps from hydromodel.")

        ### Select and build hydromaps from model ###
        # Initialise hydromaps
        ds_hydro = segments.hydromaps(hydromodel)
        rmdict = {k: v for k, v in self._MAPS.items() if k in ds_hydro.data_vars}
        self.set_hydromaps(ds_hydro.rename(rmdict))
        self._hydromaps.coords["mask"] = ds_hydro["basmsk"]

        # Build segment ID and add to hydromaps
        # Prepare delwaq pointer file for WAQ simulation (not needed with EM)
        if self.type == "WQ":
            nrofseg, da_ptid, da_ptiddown, pointer, bd_id, bd_type = segments.pointer(
                self.hydromaps, build_pointer=True
            )
            self.set_hydromaps(da_ptid, name="ptid")
            self.set_hydromaps(da_ptiddown, name="ptiddown")
            self._pointer = pointer
        else:
            nrofseg, da_ptid = segments.pointer(self.hydromaps, build_pointer=False)
            self.set_hydromaps(da_ptid, name="ptid")

        ### Initialise staticmaps with streamorder, river and slope ###
        ds_stat = (
            hydromodel.staticmaps[hydromodel._MAPS["strord"]]
            .rename("streamorder")
            .to_dataset()
        )
        ds_stat["slope"] = hydromodel.staticmaps[hydromodel._MAPS["lndslp"]]
        ds_stat["river"] = hydromodel.staticmaps[hydromodel._MAPS["rivmsk"]]
        self.set_staticmaps(ds_stat)
        self.staticmaps.coords["mask"] = self.hydromaps["modelmap"]

        ### Geometry data ###
        surface = emissions.gridarea(hydromodel.staticmaps)
        surface.raster.set_nodata(np.nan)
        geom = surface.rename("surface").to_dataset()
        # For WQ type, add surface and manning to staticmaps
        if self.type == "WQ":
            # For WQ surface in river cells should be river surface
            rivlen = hydromodel.staticmaps[hydromodel._MAPS["rivlen"]]
            rivwth = hydromodel.staticmaps[hydromodel._MAPS["rivwth"]]
            rivmsk = hydromodel.staticmaps[hydromodel._MAPS["rivmsk"]]
            geom["surface"] = xr.where(rivmsk, rivlen * rivwth, geom["surface"])
            geom["surface"].raster.set_nodata(np.nan)
            geom["manning"] = hydromodel.staticmaps["N"]
            self.set_staticmaps(geom)
        # For EM type build a pointer like object and add to self.geometry
        if self.type == "EM":
            geom["fPaved"] = hydromodel.staticmaps["PathFrac"]
            geom["fOpenWater"] = hydromodel.staticmaps["WaterFrac"]
            geom["fUnpaved"] = (
                (geom["fPaved"] * 0.0 + 1.0) - geom["fPaved"] - geom["fOpenWater"]
            )
            geom["fUnpaved"] = xr.where(geom["fUnpaved"] < 0.0, 0.0, geom["fUnpaved"])

            mask = self.staticmaps["mask"].values.flatten()
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
        if self.type == "EM":
            l7 = f"     {int(nrofseg/nrofcomp)}*01 ; EM"
        else:
            l7 = f"     {int(nrofseg/nrofcomp)}*01 ; Sfw"
        lines_ini = {
            "l1": "      ; DELWAQ_COMPLETE_ATTRIBUTES",
            "l2": " 1    ; one block with input",
            "l3": " 2    ; number of attributes, they are :",
            "l4": "     1     2",
            "l5": " 1    ; file option in this file",
            "l6": " 1    ; option without defaults",
            "l7": l7,
            "l8": " 0    ; no time dependent attributes",
        }
        for option in lines_ini:
            self.set_config("B3_attributes", option, lines_ini[option])
        # B4_nrofexch
        if self.type == "WQ":
            nrexch = self._pointer.shape[0]
        else:
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
        if self.type == "WQ":
            for i in range(len(bd_id)):
                lstr = "l" + str(i + 2)
                lines_ini.update(
                    {lstr: f"'BD_{int(bd_id[i])}' '{int(bd_id[i])}' '{bd_type[i]}'"}
                )
        for option in lines_ini:
            self.set_config("B5_boundlist", option, lines_ini[option])
        # B7_Surf
        if self.type == "EM":
            lines_ini = {
                "l1": f"PARAMETERS Surf ALL DATA {nrofseg}*1.0",
            }
            for option in lines_ini:
                self.set_config("B7_surf", option, lines_ini[option])

    def setup_monitoring(self, mon_points, mon_areas, **kwargs):
        """Setup Delwaq monitoring points and areas options.

        Adds model layers:

        * **monitoring_points** map: monitoring points segments
        * **monitoring_areas** map: monitoring areas ID
        * **B2_nrofmon.inc** config: number of monitoring points and areas

        Parameters
        ----------
        mon_points : {'None', 'segments', 'data_source', 'path to station location'}
            Either None, source from DataCatalog, path to a station location dataset
            or if segments, all segments are monitored.
        mon_areas : str {'None', 'compartments', 'subcatch'}
            Either None, subcatch from hydromodel or by compartments.
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
                if isfile(mon_points) and str(mon_points).endswith(".csv"):
                    kwargs.update(crs=self.crs)
                gdf = self.data_catalog.get_geodataframe(
                    str(mon_points), geom=self.basins, assert_gtype="Point", **kwargs
                )
                gdf = gdf.to_crs(self.crs)
                if gdf.index.size == 0:
                    self.logger.warning(
                        f"No {mon_points} gauge locations found within domain"
                    )
                else:
                    gdf.index.name = "index"
                    monpoints = self.staticmaps.raster.rasterize(
                        gdf, col_name="index", nodata=mv
                    )
            self.logger.info(f"Gauges locations read from {mon_points}")
            # Add to staticmaps
            self.set_staticmaps(monpoints.rename("monpoints"))
            # Number of monitoring points
            points = monpoints.values.flatten()
            points = points[points != mv]
            nb_points = len(points)
            # Add to staticgeoms if mon_points is not segments
            if mon_points != "segments":
                self.set_staticgeoms(gdf, name="monpoints")
        else:
            self.logger.info("No monitoring points set in the config file, skipping")
            nb_points = 0

        # Monitoring areas domain
        monareas = None
        if mon_areas == "subcatch" or mon_areas == "compartments":
            monareas = self.hydromaps["basins"]
            if mon_areas == "compartments":
                monareas = xr.where(self.hydromaps["modelmap"], 1, mv).astype(np.int32)
            # Add to staticmaps
            self.set_staticmaps(monareas.rename("monareas"))
            self.staticmaps["monareas"].attrs.update(_FillValue=mv)
            # Number or monitoring areas
            areas = monareas.values.flatten()
            areas = areas[areas != mv]
            nb_areas = len(np.unique(areas))
            # Add to staticgeoms
            gdf_areas = self.staticmaps["monareas"].raster.vectorize()
            self.set_staticgeoms(gdf_areas, name="monareas")
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

    def setup_hydrology_forcing(
        self,
        hydro_forcing_fn,
        starttime,
        endtime,
        timestepsecs,
        add_volume_offset=True,
        min_volume=0.1,
        **kwargs,
    ):
        """Setup Delwaq hydrological fluxes.

        Adds model layers:

        * **B1_timestamp.inc** config: timestamp at the beginning of the simulation.
        * **B2_outputtimes.inc** config: start and end timestamp for the delwaq outputs.
        * **B2_sysclock.inc** config: model timestep.
        * **B2_timers.inc** config: timers info (start, end, timstep...).

        In EM mode, adds:

        * **hydrology.bin** dynmap: fluxes for EM (Rainfall RunoffPav  RunoffUnp Infiltr TotalFlow) [mm]
        * **B7_hydrology.inc** config: names of fluxes included in hydrology.bin

        In WQ mode, adds:

        * **flow.dat** dynmap: fluxes for WQ (SurfaceRunoff Inwater) [m3/s]
        * **volume.dat** dynmap: water volumes [m3]
        * **area.dat** dynmap: water cross sections [m2]
        * **velocity.dat** dynmap: flow velocity [m/s]

        Parameters
        ----------
        hydro_forcing_fn : {'None', 'name in local_sources.yml'}
            Either None, or name in a local or global data_sources.yml file.

            * Required variables for EM: ['time', 'precip', 'infilt', 'runPav', 'runUnp']

            * Required variables for WQ (option1): ['time', 'run', 'vol' or 'lev', 'inwater']

            * Required variables for WQ (option2): ['time', 'runRiv' and 'runLand', 'volRiv' and 'volLand' or 'levRiv' and 'levLand', 'inwaterRiv' and 'inwaterLand']

            * Optional variable for WQ: ['inwaterInternal']
        startime : str
            Timestamp of the start of Delwaq simulation. Format: YYYY-mm-dd HH:MM:SS
        endtime : str
            Timestamp of the end of Delwaq simulation.Format: YYYY-mm-dd HH:MM:SS
        timestepsecs : int
            Model timestep in seconds.
        add_volume_offset : bool, optional
            Delwaq needs water volumes at the beginning of the timestep.
            In some models, like wflow, volumes are written at the end of the timestep and therefore
            an offset of one timestep needs to be added for consistency.
        min_volume : float, optional
            Add a minimum value to all the volumes in Delwaq to avoid zeros. Default is 0.1m3.
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
        # TODO when forcing workflow is ready nice clipping/resampling can be added
        # Normally region extent is by default exactly the same as hydromodel
        ds = self.data_catalog.get_rasterdataset(hydro_forcing_fn, geom=self.region)
        #        # Check if latitude is from south to north
        #        ys = ds.raster.ycoords
        #        resy = np.mean(np.diff(ys.values))
        #        if resy >= 0:
        #            ds = ds.reindex({ds.raster.y_dim: ds[ds.raster.y_dim][::-1]})

        # align forcing file with hydromaps
        ds = ds.raster.reproject_like(self.hydromaps)

        # Add _FillValue to the data attributes
        for dvar in ds.data_vars.keys():
            nodata = ds[dvar].raster.nodata
            if nodata is not None:
                ds[dvar].attrs.update(_FillValue=nodata)
            else:
                ds[dvar].attrs.update(_FillValue=-9999.0)

        # Update model timestepsecs attribute
        self.timestepsecs = timestepsecs
        # Select variables based on model type
        if self.type == "EM":
            ds = ds[["precip", "runPav", "runUnp", "infilt"]]
            # Add total flow
            ds["totflw"] = ds["precip"].copy()
        else:
            # If needed compute total variables from land + river
            if "runLand" in ds.data_vars and "runRiv" in ds.data_vars:
                attrs = ds["runLand"].attrs.copy()
                ds["run"] = ds["runRiv"].fillna(0) + ds["runLand"]
                ds["run"].attrs.update(attrs)
            if "inwaterLand" in ds.data_vars and "inwaterRiv" in ds.data_vars:
                attrs = ds["inwaterLand"].attrs.copy()
                ds["inwater"] = ds["inwaterRiv"].fillna(0) + ds["inwaterLand"]
                ds["inwater"].attrs.update(attrs)
            # In wflow the inwater flux contains internal water fluxes between the land and river surface waters components
            # This needs to be substracted from the inwater
            if "inwaterInternal" in ds.data_vars:
                attrs = ds["inwater"].attrs.copy()
                ds["inwater"] = ds["inwater"] - ds["inwaterInternal"]
                ds["inwater"].attrs.update(attrs)

            # If needed, compute volume from level using surface
            if "vol" not in ds.data_vars and "volRiv" not in ds.data_vars:
                if "levLand" in ds.data_vars and "levRiv" in ds.data_vars:
                    surface = emissions.gridarea(ds)
                    volL = ds["levLand"] * surface
                    rivlen = self.hydromaps["rivlen"]
                    rivwth = self.hydromaps["rivwth"]
                    volR = ds["levRiv"] * rivlen * rivwth
                    ds["vol"] = volL + volR.fillna(0) + min_volume
                    ds["vol"].attrs.update(units="m3")
                    ds["vol"].attrs.update(_FillValue=ds["levLand"].raster.nodata)
                else:
                    ds["vol"] = ds["lev"] * self.staticmaps["surface"] + min_volume
                    ds["vol"].attrs.update(units="m3")
                    ds["vol"].attrs.update(_FillValue=ds["lev"].raster.nodata)
            else:
                if "volLand" in ds.data_vars and "volRiv" in ds.data_vars:
                    attrs = ds["volLand"].attrs.copy()
                    ds["vol"] = ds["volRiv"].fillna(0) + ds["volLand"] + min_volume
                    ds["vol"].attrs.update(attrs)
                else:
                    # In order to avoid zero volumes, a basic minimum value of 0.0001 m3 is added to all volumes
                    attrs = ds["vol"].attrs.copy()
                    ds["vol"] = ds["vol"] + min_volume
                    ds["vol"].attrs.update(attrs)

            # Select variables
            ds = ds[["run", "inwater", "vol"]]

        # Select times
        # Add offset for the volumes if needed
        if add_volume_offset and self.type == "WQ":
            # Get the freq of ds and add + 1 offset
            times = pd.to_datetime(ds["time"].values)
            times.freq = pd.infer_freq(times)
            times_offset = times + times.freq
            vol = ds["vol"].copy()
            ds = ds.drop_vars("vol")
            vol["time"] = times_offset
            ds = ds.merge(vol)
        # Sell times to starttime and endtime
        ds = ds.sel(time=slice(starttime, endtime))

        # Unit conversion (from mm to m3/s)
        for dvar in ds.data_vars.keys():
            if ds[dvar].attrs.get("units") == "mm":
                attrs = ds[dvar].attrs.copy()
                surface = emissions.gridarea(ds)
                ds[dvar] = ds[dvar] * surface / (1000 * self.timestepsecs)
                ds[dvar].attrs.update(attrs)  # set original attributes
                ds[dvar].attrs.update(unit="m3/s")

        ds.coords["mask"] = xr.DataArray(
            dims=ds.raster.dims,
            coords=ds.raster.coords,
            data=self.hydromaps["modelmap"].values,
            attrs=dict(_FillValue=self.hydromaps["mask"].raster.nodata),
        )

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
        if self.type == "EM":
            lines_ini = {
                "l1": f"  1 'DDHHMMSS' 'DDHHMMSS'  ; system clock",
            }
        else:
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

        # Add var info to config
        if self.type == "EM":
            lines_ini = {
                "l1": "SEG_FUNCTIONS",
                "l2": "Rainfall RunoffPav RunoffUnp Infiltr TotalFlow ",
            }
            for option in lines_ini:
                self.set_config("B7_hydrology", option, lines_ini[option])

    def setup_sediment_forcing(
        self,
        sediment_fn,
        starttime,
        endtime,
        particle_class=["IM1", "IM2", "IM3", "IM4", "IM5"],
        **kwargs,
    ):
        """Setup Delwaq hydrological fluxes.

        Adds model layers:

        * **sediment.bin** dynmap: sediment particles from land erosion (fErodIM*) [g/timestep]
        * **B7_sediment.inc** config: names of sediment fluxes included in sediment.bin

        Parameters
        ----------
        sediment_fn : {'None', 'name in local_sources.yml'}
            Either None, or name in a local or global data_sources.yml file. Can contain several
            particule sizes (IM1, IM2 etc)

            * Required variables: ['ErodIM*']
        startime : str
            Timestamp of the start of Delwaq simulation. Format: YYYY-mm-dd HH:MM:SS
        endtime : str
            Timestamp of the end of Delwaq simulation.Format: YYYY-mm-dd HH:MM:SS
        particle_class : str list, optional
            List of particle classes to consider. By default 5 classes: ['IM1', 'IM2', 'IM3', 'IM4', 'IM5']
        """
        if sediment_fn not in self.data_catalog:
            self.logger.warning(
                f"None or Invalid source '{sediment_fn}', skipping setup_sediment_forcing."
            )
            return
        self.logger.info(f"Setting dynamic data from sediment source {sediment_fn}.")

        # read dynamic data
        ds = self.data_catalog.get_rasterdataset(sediment_fn, geom=self.region)
        # align forcing file with hydromaps
        ds = ds.raster.reproject_like(self.hydromaps)

        # Select time and particle classes
        # Sell times to starttime and endtime
        ds = ds.sel(time=slice(starttime, endtime))
        sed_vars = [f"Erod{x}" for x in particle_class]
        ds = ds[sed_vars]

        # Add basin mask
        ds.coords["mask"] = xr.DataArray(
            dims=ds.raster.dims,
            coords=ds.raster.coords,
            data=self.hydromaps["modelmap"].values,
            attrs=dict(_FillValue=self.hydromaps["mask"].raster.nodata),
        )
        # Add _FillValue to the data attributes
        for dvar in ds.data_vars.keys():
            nodata = ds[dvar].raster.nodata
            if nodata is not None:
                ds[dvar].attrs.update(_FillValue=nodata)
            else:
                ds[dvar].attrs.update(_FillValue=-9999.0)
            # Fill with zeros inside mask and keep NaN outside
            ds[dvar] = ds[dvar].fillna(0).where(ds["mask"], ds[dvar].raster.nodata)

        self.set_forcing(ds)

        lines_ini = {
            "l1": "SEG_FUNCTIONS",
            "l2": " ".join(str(x) for x in sed_vars),
        }
        for option in lines_ini:
            self.set_config("B7_sediment", option, lines_ini[option])

    def setup_emission_raster(
        self,
        emission_fn,
        scale_method="average",
        fillna_method="zero",
        fillna_value=0.0,
        area_division=False,
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
        """
        if emission_fn is None:
            self.logger.warning(
                "Source name set to None, skipping setup_emission_raster."
            )
            return
        if emission_fn not in self.data_catalog:
            self.logger.warning(
                f"Invalid source '{emission_fn}', skipping setup_emission_raster."
            )
            return

        self.logger.info(f"Preparing '{emission_fn}' map.")
        # process raster emission maps
        da = self.data_catalog.get_rasterdataset(
            emission_fn, geom=self.region, buffer=2
        )
        ds_emi = emissions.emission_raster(
            da=da,
            ds_like=self.staticmaps,
            method=scale_method,
            fillna_method=fillna_method,
            fillna_value=fillna_value,
            area_division=area_division,
            logger=self.logger,
        )
        self.set_staticmaps(ds_emi.rename(emission_fn))

    def setup_emission_vector(
        self,
        emission_fn,
        col2raster="",
        rasterize_method="value",
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
            Method to rasterize the vector data. Either {"value", "fraction", "area"}.
            If "value", the value from the col2raster is used directly in the raster.
            If "fraction", the fraction of the grid cell covered by the vector file is returned.
            If "area", the area of the grid cell covered by the vector file is returned.
        """
        if emission_fn is None:
            self.logger.warning(
                "Source name set to None, skipping setup_emission_vector."
            )
            return
        if emission_fn not in self.data_catalog:
            self.logger.warning(
                f"Invalid source '{emission_fn}', skipping setup_emission_vector."
            )
            return

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
                ds_like=self.staticmaps,
                col_name=col2raster,
                method=rasterize_method,
                mask_name="mask",
                logger=self.logger,
            )
        self.set_staticmaps(ds_emi.rename(emission_fn))

    def setup_emission_mapping(self, region_fn, mapping_fn=None):
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
        """
        if region_fn is None:
            self.logger.warning(
                "Source name set to None, skipping setup_emission_mapping."
            )
            return
        if region_fn not in self.data_catalog:
            self.logger.warning(
                f"Invalid source '{region_fn}', skipping setup_admin_bound."
                "\nCheck if {source_i} exists in data_sources.yml or local_sources.yml"
            )
        else:
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
                ds_like=self.staticmaps,
                source_name=region_fn,
                fn_map=mapping_fn,
                logger=self.logger,
            )
            rmdict = {
                k: v for k, v in self._MAPS.items() if k in ds_admin_maps.data_vars
            }
            self.set_staticmaps(ds_admin_maps.rename(rmdict))

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
                ds_like=self.staticmaps,
                col_name=f"{name}_length_sum",
                method="value",
                mask_name="mask",
                logger=self.logger,
            )
            self.set_staticmaps(da_emi.rename(f"{name}_length_sum_country"))

        ### Road statistics per segment ###
        # Filter road gdf with model mask
        mask = self.staticmaps["mask"].astype("int32").raster.vectorize()
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
                ds_like=self.staticmaps,
                variables=["length"],
                stats=["sum"],
                mask_name="mask",
                method="overlay",
            )
            # Convert from m to km and rename
            ds_roads[f"{name}_length"] = ds_roads["length_sum"] / 1000
            ds_roads[f"{name}_length"].attrs.update(_FillValue=0)
            self.set_staticmaps(ds_roads)

    # I/O
    def read(self):
        """Method to read the complete model schematization and configuration from file."""
        # self.read_config()
        self.read_staticgeoms()
        self.read_hydromaps()
        # self.read_pointer()
        self.read_staticmaps()
        # self.read_fewsadapter()
        # self.read_forcing()
        self.logger.info("Model read")

    def write(self):
        """Method to write the complete model schematization and configuration to file."""
        self.logger.info(f"Write model data to {self.root}")
        # if in r, r+ mode, only write updated components
        if self._staticmaps or not self._read:
            self.write_staticmaps()
        if self._staticgeoms or not self._read:
            self.write_staticgeoms()
        if self._config or not self._read:
            self.write_config()
        if self._hydromaps or not self._read:
            self.write_hydromaps()
        if self._pointer is not None or not self._read:
            self.write_pointer()
        if self._geometry is not None or not self._read:
            self.write_geometry()
        #        if self._fewsadapter or not self._read:
        #            self.write_fewsadapter()
        if self._forcing or not self._read:
            self.write_forcing()

    def read_staticmaps(self, crs=None, **kwargs):
        """Read staticmaps at <root/staticdata> and parse to xarray"""
        fn = join(self.root, "staticdata", "staticmaps.nc")
        if not self._write:
            # start fresh in read-only mode
            self._staticmaps = xr.Dataset()
        if fn is not None and isfile(fn):
            self.logger.info(f"Read staticmaps from {fn}")
            # FIXME: we need a smarter (lazy) solution for big models which also
            # works when overwriting / apending data in thet same source!
            ds = xr.open_dataset(
                fn, mask_and_scale=False, decode_coords="all", **kwargs
            ).load()
            ds.close()
            self.set_staticmaps(ds)

    def write_staticmaps(self):
        """Write staticmaps at <root/staticdata> in NetCDF and binary format."""
        if not self._write:
            raise IOError("Model opened in read-only mode")
        ds_out = self.staticmaps

        # Filter data with mask
        for dvar in ds_out.data_vars:
            ds_out[dvar] = ds_out[dvar].raster.mask(mask=ds_out["mask"])

        self.logger.info("Writing staticmap files.")
        # Netcdf format
        fname = join(self.root, "staticdata", "staticmaps.nc")
        # Update attributes for gdal compliance
        # ds_out = ds_out.raster.gdal_compliant(rename_dims=False)
        if os.path.isfile(fname):
            ds_out.to_netcdf(path=fname, mode="a")
        else:
            ds_out.to_netcdf(path=fname)

        # Binary format
        for dvar in ds_out.data_vars:
            if dvar != "monpoints" and dvar != "monareas":
                fname = join(self.root, "staticdata", dvar + ".dat")
                data = ds_out[dvar].values.flatten()
                mask = ds_out["mask"].values.flatten()
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

    def read_staticgeoms(self):
        """Read and staticgeoms at <root/staticgeoms> and parse to geopandas"""
        if not self._write:
            self._staticgeoms = dict()  # fresh start in read-only mode
        fns = glob.glob(join(self.root, "staticgeoms", "*.geojson"))
        if len(fns) > 1:
            self.logger.info("Reading model staticgeom files.")
        for fn in fns:
            name = basename(fn).split(".")[0]
            self.set_staticgeoms(io.open_vector(fn), name=name)

    def write_staticgeoms(self):
        """Write staticmaps at <root/staticgeoms> in model ready format"""
        # to write use self.staticgeoms[var].to_file()
        if not self._write:
            raise IOError("Model opened in read-only mode")
        if self.staticgeoms:
            self.logger.info("Writing model staticgeom to file.")
            for name, gdf in self.staticgeoms.items():
                fn_out = join(self.root, "staticgeoms", f"{name}.geojson")
                gdf.to_file(fn_out, driver="GeoJSON")

    def write_config(self):
        """Write config files in ASCII format at <root/config>."""
        if not self._write:
            raise IOError("Model opened in read-only mode")
        if self.config:
            self.logger.info("Writing model config to file.")
            for name, lines in self.config.items():
                fn_out = join(self.root, "config", f"{name}.inc")
                exfile = open(fn_out, "w")
                for lnb, line in lines.items():
                    print(line, file=exfile)
                exfile.close()

    def read_hydromaps(self, crs=None, **kwargs):
        """Read hydromaps at <root/hydromodel> and parse to xarray"""
        if self._read and "chunks" not in kwargs:
            kwargs.update(chunks={"y": -1, "x": -1})
        fns = glob.glob(join(self.root, "hydromodel", f"*.tif"))
        if len(fns) > 0:
            self._hydromaps = io.open_mfraster(fns, **kwargs)
        if self._hydromaps.raster.crs is None and crs is not None:
            self.set_crs(crs)
        self._hydromaps["modelmap"] = self._hydromaps["modelmap"].astype(bool)
        self._hydromaps.coords["mask"] = self._hydromaps["modelmap"]

    def write_hydromaps(self):
        """Write hydromaps at <root/hydromodel> in PCRaster maps format."""
        if not self._write:
            raise IOError("Model opened in read-only mode")
        ds_out = self.hydromaps
        if self._read:
            # only write loaded maps in 'r+' mode
            dvars = [
                dvar
                for dvar in self.hydromaps.data_vars.keys()
                if isinstance(self.hydromaps[dvar].data, np.ndarray)
            ]
            if len(dvars) > 0:
                ds_out = ds_out[dvars]
                self.logger.debug(f"Updated maps: {dvars}")
            else:
                self.logger.warning(f"No updated maps. Skipping writing to file.")
                return
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
        """Read Delwaq pointer file"""
        raise NotImplementedError()

    def write_pointer(self):
        """Write pointer at <root/dynamicdata> in ASCII and binary format."""
        if not self._write:
            raise IOError("Model opened in read-only mode")
        if self._pointer is not None:
            self.logger.info("Writting pointer file in root/config")
            fname = join(self.root, "config", "B4_pointer")
            # Write ASCII file
            exfile = open((fname + ".inc"), "w")
            print(";Pointer for WAQ simulation in Surface Water", file=exfile)
            print(";nr of pointers is: ", str(self._pointer.shape[0]), file=exfile)
            np.savetxt(exfile, self._pointer, fmt="%10.0f")
            exfile.close()

            # Write binary file
            f = open((fname + ".poi"), "wb")
            for i in range(self._pointer.shape[0]):
                f.write(struct.pack("4i", *np.int_(self._pointer[i, :])))
            f.close()

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

    def read_forcing(self):
        """Read and forcing at <root/?/> and parse to dict of xr.DataArray"""
        if not self._write:
            # start fresh in read-only mode
            self._forcing = dict()
        # raise NotImplementedError()

    def write_forcing(self, write_nc=False):
        """Write staticmaps at <root/staticdata> in binary format and NetCDF (if write_nc is True)."""
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
            if os.path.isfile(fname):
                ds_out.to_netcdf(path=fname, mode="a")
            else:
                ds_out.to_netcdf(path=fname)

        # Binary format
        timesteps = np.arange(0, len(ds_out.time.values))
        timestepstamp = np.arange(
            0, (len(ds_out.time.values) + 1) * int(self.timestepsecs), self.timestepsecs
        )
        if self.type == "EM":
            for i in timesteps:
                self.logger.info(
                    f"Writting dynamic data for timestep {i+1}/{timesteps[-1]+1}"
                )
                fname = join(self.root, "dynamicdata", "hydrology.bin")
                datablock = []
                for dvar in ["precip", "runPav", "runUnp", "infilt", "totflw"]:
                    nodata = ds_out[dvar].raster.nodata
                    data = ds_out[dvar].isel(time=i).values.flatten()
                    data = data[data != nodata]
                    datablock = np.append(datablock, data)
                self.dw_WriteSegmentOrExchangeData(
                    timestepstamp[i], fname, datablock, 1, WriteAscii=False
                )
        else:
            for i in timesteps:
                self.logger.info(
                    f"Writting dynamic data for timestep {i+1}/{len(timesteps)}"
                )
                # Flow
                flname = join(self.root, "dynamicdata", "flow.dat")
                flowblock = []
                for dvar in ["run", "inwater"]:
                    nodata = ds_out[dvar].raster.nodata
                    data = ds_out[dvar].isel(time=i).values.flatten()
                    data = data[data != nodata]
                    flowblock = np.append(flowblock, data)
                self.dw_WriteSegmentOrExchangeData(
                    timestepstamp[i], flname, flowblock, 1, WriteAscii=False
                )
                # volume
                voname = join(self.root, "dynamicdata", "volume.dat")
                nodata = ds_out["vol"].raster.nodata
                data = ds_out["vol"].isel(time=i).values.flatten()
                data = data[data != nodata]
                self.dw_WriteSegmentOrExchangeData(
                    timestepstamp[i], voname, data, 1, WriteAscii=False
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

        # Add waqgeom.nc file to allow Delwaq to save outputs in nc format
        self.logger.info("Writting waqgeom.nc file")
        self.dw_WriteWaqGeom()

    def read_states(self):
        """Read states at <root/?/> and parse to dict of xr.DataArray"""
        if not self._write:
            # start fresh in read-only mode
            self._states = dict()
        # raise NotImplementedError()

    def write_states(self):
        """write states at <root/?/> in model ready format"""
        if not self._write:
            raise IOError("Model opened in read-only mode")
        raise NotImplementedError()

    def read_results(self):
        """Read results at <root/?/> and parse to dict of xr.DataArray"""
        if not self._write:
            # start fresh in read-only mode
            self._results = dict()
        # raise NotImplementedError()

    def write_results(self):
        """write results at <root/?/> in model ready format"""
        if not self._write:
            raise IOError("Model opened in read-only mode")
        raise NotImplementedError()

    ## DELWAQ specific data and methods

    @property
    def basins(self):
        if "basins" in self.staticgeoms:
            gdf = self.staticgeoms["basins"]
        elif "basins" in self.hydromaps:
            gdf = self.hydromaps["basins"].raster.vectorize()
            gdf.crs = pyproj.CRS.from_user_input(self.crs)
            self.set_staticgeoms(gdf, name="basins")
        return gdf

    @property
    def hydromaps(self):
        """xarray.dataset representation of all hydrology maps"""
        if len(self._hydromaps) == 0:
            if self._read:
                self.read_hydromaps()
            else:
                raise ValueError("No hydromaps defined")
        return self._hydromaps

    def set_hydromaps(self, data, name=None):
        """Add data to hydromaps re-using the set_staticmaps method"""
        if name is None:
            if isinstance(data, xr.DataArray) and data.name is not None:
                name = data.name
            elif not isinstance(data, xr.Dataset):
                raise ValueError("Setting a map requires a name")
        elif name is not None and isinstance(data, xr.Dataset):
            data_vars = list(data.data_vars)
            if len(data_vars) == 1 and name not in data_vars:
                data = data.rename_vars({data_vars[0]: name})
            elif name not in data_vars:
                raise ValueError("Name not found in DataSet")
            else:
                data = data[[name]]
        if isinstance(data, xr.DataArray):
            data.name = name
            data = data.to_dataset()
        if len(self._hydromaps) == 0:  # new data
            if not isinstance(data, xr.Dataset):
                raise ValueError("First parameter map(s) should xarray.Dataset")
            self._hydromaps = data
        else:
            if isinstance(data, np.ndarray):
                if data.shape != self.shape:
                    raise ValueError("Shape of data and staticmaps do not match")
                data = xr.DataArray(dims=self.dims, data=data, name=name).to_dataset()
            for dvar in data.data_vars.keys():
                if dvar in self._hydromaps:
                    if not self._write:
                        raise IOError(
                            f"Cannot overwrite staticmap {dvar} in read-only mode"
                        )
                    elif self._read:
                        self.logger.warning(f"Overwriting staticmap: {dvar}")
                self._hydromaps[dvar] = data[dvar]

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
        Writes a timestep to a segment/exchange data file (appends to an existing
        file or creates a new one).

        Input:
            - time - timestep number for this timestep
            - fname - File path of the segment/exchange data file</param>
            - datablock - array with data
            - boundids to write more than 1 block
            - WriteAscii - if True to make a copy in an ascii checkfile
            - mode - {"a", "w"} Force the writting mode, append or overwrite existing files.

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
        Writes monitoring files and config in ASCII format.

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

                print(f"'{i}'        {NTOT}", file=exfile)
                np.savetxt(exfile, areai_1, fmt="%10.20s")
                np.savetxt(exfile, areai_2, fmt="%10.20s")
            exfile.close()
        else:
            fname = join(self.root, "config", "B2_monareas.inc")
            exfile = open(fname, "w")
            print(";Written by hydroMT: no monitoring areas were set.", file=exfile)
            exfile.close()

    def dw_WriteWaqGeom(self):
        """Writes Delwaq netCDF geometry file (config/B3_waqgeom.nc)."""
        ptid = self.hydromaps["ptid"].copy()
        ptid_mv = ptid.raster.nodata
        # PCR cell id's start at 1, we need it zero based, and NaN set to -1
        ptid = xr.where(ptid == ptid_mv, -1, ptid - 1)
        np_ptid = ptid.values
        # Wflow map dimensions
        m, n = np_ptid.shape
        # Number of segments in horizontal dimension
        nosegh = int(np.max(np_ptid)) + 1  # because of the zero based
        # Get LDD map
        np_ldd = self.hydromaps["ldd"].values
        np_ldd[np_ldd == self.hydromaps["ldd"].raster.nodata] = 0

        # Get  grid coordinates
        res_x, res_y = self.res
        xcoords = ptid.raster.xcoords.values
        ycoords = ptid.raster.ycoords.values
        left = xcoords - res_x / 2.0
        right = xcoords + res_x / 2.0
        top = ycoords - res_y / 2.0
        bottom = ycoords + res_y / 2.0

        xxul = np.tile(left, len(ycoords)).reshape((len(ycoords), len(xcoords)))
        xxlr = np.tile(right, len(ycoords)).reshape((len(ycoords), len(xcoords)))
        yyul = np.repeat(top, len(xcoords)).reshape((len(ycoords), len(xcoords)))
        yylr = np.repeat(bottom, len(xcoords)).reshape((len(ycoords), len(xcoords)))

        # Waqgeom dimensions
        n_net_node = 0
        n_net_link = 0
        n_net_link_pts = 2
        n_net_elem = nosegh
        n_net_elem_max_node = 4  # all elements are rectangles
        n_flow_link = nosegh - 1  # one per element, except for outlet
        n_flow_link_pts = 2

        # Prepare waqgeom data structures
        nodes_x = []
        nodes_y = []
        nodes_z = []
        net_links = []
        elem_nodes = np.zeros((n_net_elem, n_net_elem_max_node), dtype=np.int32)
        face_x_bnd = np.zeros((n_net_elem, n_net_elem_max_node), dtype=np.double)
        face_y_bnd = np.zeros((n_net_elem, n_net_elem_max_node), dtype=np.double)
        flow_links = np.zeros((n_flow_link, n_flow_link_pts), dtype=np.int32)
        flow_link_type = np.repeat(2, n_flow_link)
        flow_link_x = np.zeros((n_flow_link), dtype=np.float32)
        flow_link_y = np.zeros((n_flow_link), dtype=np.float32)

        # prepare all coordinates for grid csv
        nodes_x_all = []
        nodes_y_all = []
        # Keep track of nodes and links as dataset grows
        i_node = 0  # index of last node
        i_flink = 0  # index of last flow link

        # Helper function
        def add_node(i, j, corner):
            # Get coordinates
            if corner == UL:
                x = xxul[i, j]
                y = yyul[i, j]
            elif corner == LR:
                x = xxlr[i, j]
                y = yylr[i, j]
            elif corner == UR:
                x = xxlr[i, j]
                y = yyul[i, j]
            elif corner == LL:
                x = xxul[i, j]
                y = yylr[i, j]
            else:
                assert 0
            # Add node coordinates
            nodes_x.append(x)
            nodes_y.append(y)
            nodes_z.append(0)

        def add_all_nodes(i, j, corner):
            # Get coordinates
            if corner == UL:
                x = xxul[i, j]
                y = yyul[i, j]
            elif corner == LR:
                x = xxlr[i, j]
                y = yylr[i, j]
            elif corner == UR:
                x = xxlr[i, j]
                y = yyul[i, j]
            elif corner == LL:
                x = xxul[i, j]
                y = yylr[i, j]
            else:
                assert 0
            # Add node coordinates
            nodes_x_all.append(x)
            nodes_y_all.append(y)

        # Cell corners
        UL, UR, LR, LL = 0, 1, 2, 3

        # Process all cells from upper-left to lower-right
        for i in range(m):
            for j in range(n):
                # Current element index
                i_elem = int(np_ptid[i, j])
                if i_elem < 0:
                    # Skip inactive segment
                    continue

                # Get index of neighbouring elements that could have been processed before
                if i == 0:
                    i_elem_up_left = -1
                    i_elem_up = -1
                    i_elem_up_right = -1
                elif j == 0:
                    i_elem_up_left = -1
                    i_elem_up = int(np_ptid[i - 1, j])
                    i_elem_up_right = int(np_ptid[i - 1, j + 1])
                elif j == n - 1:
                    i_elem_up_left = int(np_ptid[i - 1, j - 1])
                    i_elem_up = int(np_ptid[i - 1, j])
                    i_elem_up_right = -1
                else:
                    i_elem_up_left = int(np_ptid[i - 1, j - 1])
                    i_elem_up = int(np_ptid[i - 1, j])
                    i_elem_up_right = int(np_ptid[i - 1, j + 1])

                if j == 0:
                    i_elem_left = -1
                else:
                    i_elem_left = int(np_ptid[i, j - 1])

                # Update nodes:
                # If left or upper neighbours are active, some nodes of current cell
                # have been added already.

                # UL node
                if i_elem_left < 0 and i_elem_up_left < 0 and i_elem_up < 0:
                    add_node(i, j, UL)
                    elem_nodes[i_elem, UL] = i_node
                    i_node += 1
                elif i_elem_left >= 0:
                    elem_nodes[i_elem, UL] = elem_nodes[i_elem_left, UR]
                elif i_elem_up_left >= 0:
                    elem_nodes[i_elem, UL] = elem_nodes[i_elem_up_left, LR]
                elif i_elem_up >= 0:
                    elem_nodes[i_elem, UL] = elem_nodes[i_elem_up, LL]

                # UR node
                if i_elem_up < 0 and i_elem_up_right < 0:
                    add_node(i, j, UR)
                    elem_nodes[i_elem, UR] = i_node
                    i_node += 1
                elif i_elem_up >= 0:
                    elem_nodes[i_elem, UR] = elem_nodes[i_elem_up, LR]
                elif i_elem_up_right >= 0:
                    elem_nodes[i_elem, UR] = elem_nodes[i_elem_up_right, LL]
                if i_elem_up < 0:
                    # add UL-UR link
                    net_links.append((elem_nodes[i_elem, UL], elem_nodes[i_elem, UR]))

                # LL node
                if i_elem_left < 0:
                    add_node(i, j, LL)
                    elem_nodes[i_elem, LL] = i_node
                    i_node += 1
                    # add UL-LL link
                    net_links.append((elem_nodes[i_elem, UL], elem_nodes[i_elem, LL]))
                else:
                    elem_nodes[i_elem, LL] = elem_nodes[i_elem_left, LR]

                # LR node
                add_node(i, j, LR)
                add_all_nodes(i, j, LR)
                elem_nodes[i_elem, LR] = i_node
                i_node += 1
                # add LL-LR link
                net_links.append((elem_nodes[i_elem, LL], elem_nodes[i_elem, LR]))
                # add UR-LR link
                net_links.append((elem_nodes[i_elem, UR], elem_nodes[i_elem, LR]))

                # Update flow links based on local drain direction
                # TODO: diagonal flow links between cells that have only one node in common?

                direction = np_ldd[i, j]
                i_other = -1
                if direction == 1:
                    i_other = np_ptid[i + 1, j - 1]  # to lower left
                elif direction == 2:
                    i_other = np_ptid[i + 1, j]  # to lower
                elif direction == 3:
                    i_other = np_ptid[i + 1, j + 1]  # to lower right
                elif direction == 4:
                    i_other = np_ptid[i, j - 1]  # to left
                elif direction == 6:
                    i_other = np_ptid[i, j + 1]  # to right
                elif direction == 7:
                    i_other = np_ptid[i - 1, j - 1]  # to upper right
                elif direction == 8:
                    i_other = np_ptid[i - 1, j]  # to upper
                elif direction == 9:
                    i_other = np_ptid[i - 1, j + 1]  # to upper left
                if i_other >= 0:
                    flow_links[i_flink, :] = i_elem, i_other
                    i_flink += 1

        # Convert data to numpy arrays
        nodes_x_all = np.array(nodes_x_all)
        nodes_y_all = np.array(nodes_y_all)

        nodes_x = np.array(nodes_x)
        nodes_y = np.array(nodes_y)
        nodes_z = np.array(nodes_z)
        net_links = np.array(net_links)

        # Proces all cells to derive mesh_face_x_bnd and mesh_face_y_bnd
        for bnd in range(0, n_net_elem):
            # UL, UR, LR, LL = 0, 1, 2, 3
            face_x_bnd[bnd, UL] = nodes_x[elem_nodes[bnd][UL]]
            face_x_bnd[bnd, UR] = nodes_x[elem_nodes[bnd][UR]]
            face_x_bnd[bnd, LR] = nodes_x[elem_nodes[bnd][LR]]
            face_x_bnd[bnd, LL] = nodes_x[elem_nodes[bnd][LL]]
            face_y_bnd[bnd, UL] = nodes_y[elem_nodes[bnd][UL]]
            face_y_bnd[bnd, UR] = nodes_y[elem_nodes[bnd][UR]]
            face_y_bnd[bnd, LR] = nodes_y[elem_nodes[bnd][LR]]
            face_y_bnd[bnd, LL] = nodes_y[elem_nodes[bnd][LL]]

        # Update dimensions
        n_net_node = nodes_x.shape[0]
        n_net_link = net_links.shape[0]

        # Info for the global xarray attributes
        time_string = t.strftime("%b %d %Y, %H:%M:%S")
        offset_s = -t.altzone
        offset_m = int((offset_s % 3600) / 60)
        offset_h = int((offset_s / 60 - offset_m) / 60)
        time_string2 = t.strftime("%Y-%m-%dT%H:%M:%S") + "+%02i%02i" % (
            offset_h,
            offset_m,
        )

        # Create xr dataset
        # TODO: mesh_face_x and mesh_face_y now based on LR node, correct?
        # TODO: add mesh_edge* to dataset
        ds_out = xr.Dataset(
            data_vars=dict(
                mesh=(["dim"], [-2147483647]),
                # projected_coordinate_system=(["dim"], [-2147483647]),
                # NetNode_x=(["nNetNode"], nodes_x),
                # NetNode_y=(["nNetNode"], nodes_y),
                # NetNode_z=(["nNetNode"], nodes_z),
                mesh_node_x=(["nNetNode"], nodes_x),
                mesh_node_y=(["nNetNode"], nodes_y),
                mesh_node_z=(["nNetNode"], nodes_z),
                mesh_face_x=(["nmesh_face"], nodes_x_all),
                mesh_face_y=(["nmesh_face"], nodes_y_all),
                NetLink=(["nNetLink", "nNetLinkPts"], (net_links + 1)),
                mesh_face_nodes=(
                    ["nmesh_face", "max_mesh_face_nodes"],
                    (elem_nodes + 1),
                ),
                mesh_face_x_bnd=(["nmesh_face", "max_mesh_face_nodes"], (face_x_bnd)),
                mesh_face_y_bnd=(["nmesh_face", "max_mesh_face_nodes"], (face_y_bnd)),
                FlowLink=(["nFlowLink", "nFlowLinkPts"], (flow_links + 1)),
                FlowLinkType=(["nFlowLink"], flow_link_type),
                FlowLink_xu=(["nFlowLink"], flow_link_x),
                FlowLink_yu=(["nFlowLink"], flow_link_y),
            ),
            coords=dict(
                # dim = [1],
                projected_coordinate_system=[-2147483647],
                # nNetNode = [n_net_node],
                # nNetLink = [n_net_link],
                # nNetLinkPts = [n_net_link_pts],
                # nNetElem = [n_net_elem],
                # nNetElemMaxNode = [n_net_elem_max_node],
                # nFlowLink= [n_flow_link],
                # nFlowLinkPts = [n_flow_link_pts],
            ),
            attrs=dict(
                institution="Deltares",
                references="http://www.deltares.nl",
                source=f"Wflow, Deltares, {time_string}.",
                history=f"Created on {time_string2}, wflow_delwaq.py",
                Conventions="CF-1.6 UGRID-0.9",
            ),
        )
        # Update variables attributes
        ds_out["mesh"].attrs.update(
            dict(
                long_name="Delft3D FM aggregated mesh",
                cf_role="mesh_topology",
                topology_dimension=2,
                node_coordinates="NetNode_x NetNode_y",
                face_node_connectivity="mesh_face_nodes",
                edge_node_connectivity="NetLink",
                edge_face_connectivity="FlowLink",
                face_dimension="nmesh_face",
                edge_dimension="nNetLink",
                node_dimension="nNetNode",
                face_coordinates="mesh_face_x mesh_face_y",
                edge_coordinates="FlowLink_xu FlowLink_yu",
            )
        )
        epsg_nb = int(self.crs.to_epsg())
        ds_out["projected_coordinate_system"].attrs.update(
            dict(
                epsg=epsg_nb,
                grid_mapping_name="Unknown projected",
                longitude_of_prime_meridian=0.0,
                inverse_flattening=298.257223563,
                epsg_code=f"{self.crs}",
                value="value is equal to EPSG code",
            )
        )
        # ds_out["NetNode_x"].attrs.update(
        #     dict(
        #         units="degrees_east",
        #         standard_name="longitude",
        #         long_name="longitude",
        #         mesh="mesh",
        #         location="node",
        #     )
        # )
        # ds_out["NetNode_y"].attrs.update(
        #     dict(
        #         units="degrees_north",
        #         standard_name="latitude",
        #         long_name="latitude",
        #         mesh="mesh",
        #         location="node",
        #     )
        # )
        # ds_out["NetNode_z"].attrs.update(
        #     dict(
        #         units="m",
        #         positive="up",
        #         standard_name="sea_floor_depth",
        #         long_name="Bottom level at net nodes (flow element's corners)",
        #         coordinates="NetNode_x NetNode_y",
        #     )
        # )
        ds_out["mesh_node_x"].attrs.update(
            dict(
                units="degrees_east",
                standard_name="longitude",
                long_name="longitude",
                mesh="mesh",
                location="node",
            )
        )
        ds_out["mesh_node_y"].attrs.update(
            dict(
                units="degrees_north",
                standard_name="latitude",
                long_name="latitude",
                mesh="mesh",
                location="node",
            )
        )
        ds_out["mesh_face_x"].attrs.update(
            dict(
                units="degrees_east",
                standard_name="longitude",
                long_name="longitude",
                mesh="mesh",
                location="face",
            )
        )
        ds_out["mesh_face_y"].attrs.update(
            dict(
                units="degrees_north",
                standard_name="latitude",
                long_name="latitude",
                mesh="mesh",
                location="face",
            )
        )
        ds_out["NetLink"].attrs.update(
            dict(
                long_name="link between two netnodes",
                start_index=1,
            )
        )
        ds_out["mesh_face_nodes"].attrs.update(
            dict(
                long_name="Mapping from every face to its corner nodes (counterclockwise)",
                cf_role="face_node_connectivity",
                mesh="mesh",
                location="face",
                start_index=1,
                _FillValue=0,
            )
        )
        ds_out["mesh_face_x_bnd"].attrs.update(
            dict(
                units="m",
                standard_name="projection_x_coordinate",
                long_name="x-coordinate bounds of 2D mesh face (i.e. corner coordinates)",
                mesh="mesh",
                location="face",
                start_index=1,
            )
        )
        ds_out["mesh_face_y_bnd"].attrs.update(
            dict(
                units="m",
                standard_name="projection_y_coordinate",
                long_name="y-coordinate bounds of 2D mesh face (i.e. corner coordinates)",
                mesh="mesh",
                location="face",
                start_index=1,
            )
        )
        ds_out["FlowLink"].attrs.update(
            dict(
                long_name="link/interface between two flow elements",
                start_index=1,
            )
        )
        ds_out["FlowLinkType"].attrs.update(
            dict(
                long_name="type of flowlink",
                valid_range=[1, 2],
                flag_values=[1, 2],
                flag_meanings="link_between_1D_flow_elements link_between_2D_flow_elements",
            )
        )
        ds_out["FlowLink_xu"].attrs.update(
            dict(
                units="degrees_east",
                standard_name="longitude",
                long_name="x-Coordinate of velocity point on flow link.",
            )
        )
        ds_out["FlowLink_yu"].attrs.update(
            dict(
                units="degrees_north",
                standard_name="latitude",
                long_name="y-Coordinate of velocity point on flow link.",
            )
        )

        # Write the waqgeom.nc file
        fname = join(self.root, "config", "B3_waqgeom.nc")
        ds_out.to_netcdf(path=fname, mode="w", format="NETCDF3_CLASSIC")
