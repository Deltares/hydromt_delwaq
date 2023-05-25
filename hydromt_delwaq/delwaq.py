"""Implement delwaq model class"""

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
import matplotlib.pyplot as plt
import xugrid as xu

import hydromt
from hydromt.models.model_api import Model
from hydromt import workflows, flw, io


from hydromt_wflow.wflow import WflowModel

from .workflows import emissions, segments
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
        # "mask": "modelmap",
        "flwdir": "ldd",
        "lndslp": "slope",
        "N": "manning",
        "rivmsk": "river",
        "strord": "streamorder",
        "thetaS": "porosity",
    }
    _FORCING = {
        "temp": "tempair",
        "temp_dew": "temp_dew",
        "ssr": "radsw",
        "wind": "vwind",
        "tcc": "cloudfrac",
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
        self._pointer = (
            None  # dictionnary of pointer values and related model attributes
        )
        self._geometry = None
        self._fewsadapter = None

        self.timestepsecs = 86400

    def setup_basemaps(
        self,
        region,
        mask="basins",
        compartments=["sfw"],
        boundaries=["bd"],
        fluxes=["sfw>sfw", "bd>sfw"],
        comp_attributes=["0"],
        maps=["rivmsk", "lndslp", "strord", "N"],
    ):
        """Setup the delwaq model schematization using the hydromodel region and
        resolution. 
        
        Maps used and derived from the hydromodel are stored in a specific\ 
        hydromodel attribute. Build a D-Water Quality ("WQ") case.\
        The pointer will be created based on the ``fluxes`` list
        between model compartments and boundaries.
        
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
            Mask to use to define Delwaq segments. Either "basins" (default) or "rivers".        
        compartments : list of str, optional
            List of names of compartments to include. By default one for surface waters called 'sfw'.
        boundaries: list of str, optional
            List of names of boundaries to include. By default a unique boundary called 'bd'.
        fluxes: list of str
            List of fluxes to include between compartments/boundaries. Name convention is '{compartment_name}>{boundary_name}'
            for a flux from a compartment to a boundary, ex 'sfw>bd'. By default ['sfw>sfw', 'bd>sfw'] for runoff and inwater.
            Names in the fluxes list should match name in the hydrology_fn source in setup_hydrology_forcing.
        comp_attributes: list of int
            Attribute 1 value of the B3_attributes config file. 1 or 0 for surface water. Also used to compute surface variable.
        maps: list of str
            List of variables from hydromodel to add to staticmaps. 
            By default ['rivmsk', 'lndslp', 'strord', 'N'].
        """
        # Only list of str allowed in config
        comp_attributes = [int(i) for i in comp_attributes]
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

        self.logger.info(f"Preparing WQ basemaps from hydromodel.")

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
            compartments=compartments,
            boundaries=boundaries,
            fluxes=fluxes,
        )
        self.set_hydromaps(da_ptid, name="ptid")
        self.set_hydromaps(da_ptiddown, name="ptiddown")
        # Initialise pointer object with schematisation attributes
        self.set_pointer(pointer, name="pointer")
        self.set_pointer(nrofseg, name="nrofseg")
        self.set_pointer(compartments, name="compartments")
        self.set_pointer(bd_type, name="boundaries")
        self.set_pointer(fluxes, name="fluxes")

        ### Initialise staticmaps with streamorder, river and slope ###
        ds_stat = segments.maps_from_hydromodel(hydromodel, compartments, maps=maps)
        mask = segments.extend_comp_with_duplicates(da_mask.to_dataset(), compartments)
        ds_stat["mask"] = mask[da_mask.name]
        ds_stat = ds_stat.set_coords("mask")
        rmdict = {k: v for k, v in self._MAPS.items() if k in ds_stat.data_vars}
        self.set_staticmaps(ds_stat.rename(rmdict))

        ### Add geometry ###
        ds_geom = segments.geometrymaps(
            hydromodel,
            compartments=compartments,
            comp_attributes=comp_attributes,
        )
        self.set_staticmaps(ds_geom)

        ### Config ###
        nrofcomp = len(compartments)
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
        }
        nl = 7
        for i in range(nrofcomp):
            lines_ini[
                f"l{nl}"
            ] = f"     {int(nrofseg/nrofcomp)}*{comp_attributes[i]}1 ; {compartments[i]}"
            nl += 1
        lines_ini[f"l{nl}"] = " 0    ; no time dependent attributes"
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
            "l1": f";'NodeID' 'Number' 'Type'",
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
            else:  # TODO update for several compartements
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
        if mon_areas != None:
            monareas = self.hydromaps["basins"].copy()
            monareas_tot = []
            if mon_areas == "compartments":
                nb_areas = self.nrofcomp
                for i in np.arange(1, self.nrofcomp + 1):
                    monareas_tot.append(
                        xr.where(self.hydromaps["mask"], i, mv).astype(np.int32)
                    )
            elif mon_areas == "subcatch":  # subcatch
                # Number or monitoring areas
                areas = monareas.values.flatten()
                areas = areas[areas != mv]
                nb_sub = len(np.unique(areas))
                nb_areas = nb_sub * self.nrofcomp
                for i in range(self.nrofcomp):
                    monareas_tot.append(monareas + (i * nb_sub))
            else:  # riverland
                # seperate areas for land cells (1) and river cells (2)
                lr_areas = xr.where(
                    self.hydromaps["river"],
                    2,
                    xr.where(self.hydromaps["basins"], 1, mv),
                ).astype(np.int32)
                # Number or monitoring areas
                ##plt.imshow(lr_areas)
                ##plt.show()
                ##plt.savefig('filename.png')
                areas = lr_areas.values.flatten()
                areas = areas[areas != mv]
                nb_sub = len(np.unique(areas))
                nb_areas = nb_sub * self.nrofcomp
                for i in np.arange(1, self.nrofcomp + 1):
                    monareas_tot.append(lr_areas + ((i - 1) * nb_sub))
            monareas = xr.concat(
                monareas_tot,
                pd.Index(np.arange(1, self.nrofcomp + 1, dtype=int), name="comp"),
            ).transpose("comp", ...)
            # Add to staticmaps
            self.set_staticmaps(monareas.rename("monareas"))
            self.staticmaps["monareas"].attrs.update(_FillValue=mv)
            self.staticmaps["monareas"].attrs["mon_areas"] = mon_areas

            # Add to staticgeoms (only first compartments)
            gdf_areas = self.staticmaps["monareas"].sel(comp=1).raster.vectorize()
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

        As the fluxes order should precisely macth the pointer defined in setup_basemaps, the variables names
        in ``hydro_forcing_fn`` should match names defined in the ``fluxes`` argument of setup_basemaps.
        These names should also have been saved in the file config/B7_fluxes.inc.

        If several sub-variables in ``hydro_forcing_fn`` need to be summed up to get the expected flux in pointer,
        they can be named {flux_name_in_pointer}_{number} (eg "sfw>sfw_1" and "sfw>sfw_2" to get "sfw>sfw") and the
        function will sum them on the fly. To remove (- instead of +) use unit_mult attribute of tha data catalog
        with -1 for the sub-variables of interest.

        Unit conversions are possible (after the sum!) from mm/day to m3/s for fluxes, volumes should be provided directly in m3.

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
            Names of fluxes should match those as set in the ``fluxes`` list of setup_basemaps methods (pointer creation).

            * Required variables (to be deprecated): ['time', 'run', 'vol' or 'lev', 'inwater']

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
        # Normally region extent is by default exactly the same as hydromodel
        ds = self.data_catalog.get_rasterdataset(hydro_forcing_fn, geom=self.region)
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

        # Copy of ds to be filled
        dsvar = [v for v in ds.data_vars if v.startswith(self.fluxes[0])][0]
        ds_out = hydromt.raster.full_like(ds[dsvar], lazy=True).to_dataset()
        ds_out = ds_out.sel(time=slice(starttime, endtime))
        # Array of zeros
        da_zeros = ds[dsvar] * 0.0
        da_zeros.attrs.update(units="m3/s")
        da_zeros.raster.set_nodata(-999.0)

        ### Fluxes ###
        for flux in self.fluxes:
            # Check if the flux is split into several variables
            fl_vars = [v for v in ds.data_vars if v.startswith(flux)]
            if len(fl_vars) > 0:  # flux not in ds
                attrs = ds[fl_vars[0]].attrs.copy()
            else:
                self.logger.warning(
                    f"Flux {flux} not found in hydro_forcing_fn. Using zeros."
                )
                ds[flux] = da_zeros
                attrs = da_zeros.attrs.copy()
            if len(fl_vars) > 1:  # need to sum
                ds[flux] = ds[fl_vars].fillna(0).to_array().sum("variable")
            # Unit conversion (from mm to m3/s)
            if ds[flux].attrs.get("units") == "mm":
                surface = emissions.gridarea(ds)
                ds[flux] = ds[flux] * surface / (1000 * self.timestepsecs)
            ds_out[flux] = ds[flux].sel(time=slice(starttime, endtime))
            ds_out[flux].attrs.update(attrs)
            ds_out[flux].attrs.update(unit="m3/s")

        ### Volumes ###
        # Add offset for the volumes if needed
        if add_volume_offset:
            # Get the freq of ds and add + 1 offset
            times = pd.to_datetime(ds["time"].values)
            times.freq = pd.infer_freq(times)
            times_offset = times + times.freq
        for vol in self.compartments:
            # Check if the flux is split into several variables
            vol_vars = [v for v in ds.data_vars if v.startswith(vol)]
            attrs = ds[vol_vars[0]].attrs.copy()
            if len(vol_vars) > 1:  # need to sum
                ds[vol] = ds[vol_vars].fillna(0).to_array().sum("variable")
            # Unit conversion (from mm to m3)
            if ds[vol].attrs.get("units") == "mm":
                surface = emissions.gridarea(ds)
                ds[vol] = ds[vol] * surface / (1000 * self.timestepsecs)
            # In order to avoid zero volumes, a basic minimum value of 0.0001 m3 is added to all volumes
            ds[vol] = ds[vol] + min_volume
            # Add offset for the volumes if needed
            if add_volume_offset:
                da_vol = ds[vol].copy()
                ds = ds.drop_vars(vol)
                da_vol["time"] = times_offset
                # ds = ds.merge(da_vol)
            else:
                da_vol = ds[vol]
            ds_out[vol] = da_vol.sel(time=slice(starttime, endtime))
            ds_out[vol].attrs.update(attrs)
            ds_out[vol].attrs.update(unit="m3")

        # Select variables, needed??
        variables = self.fluxes.copy()
        variables.extend(self.compartments)
        ds_out = ds_out[variables]

        # Sell times to starttime and endtime
        # ds = ds.sel(time=slice(starttime, endtime))

        ds_out.coords["mask"] = xr.DataArray(
            dims=ds_out.raster.dims,
            coords=ds_out.raster.coords,
            data=self.hydromaps["mask"].values,
            attrs=dict(_FillValue=self.hydromaps["mask"].raster.nodata),
        )

        self.set_forcing(ds_out)

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
        # Select time and particle classes
        # Sell times to starttime and endtime
        ds = ds.sel(time=slice(starttime, endtime))
        sed_vars = [f"Erod{x}" for x in particle_class]
        ds = ds[sed_vars]
        # align forcing file with hydromaps
        ds = ds.raster.reproject_like(self.hydromaps)

        # Add basin mask
        ds.coords["mask"] = xr.DataArray(
            dims=ds.raster.dims,
            coords=ds.raster.coords,
            data=self.hydromaps["mask"].values,
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

        # Add zeros to the non surface water compartments
        comp_sfw = segments.sfwcomp(self.compartments, self.config)
        ds = segments.extend_comp_with_zeros(
            ds1c=ds, comp_ds1c=comp_sfw, compartments=self.compartments
        )

        self.set_forcing(ds)

        lines_ini = {
            "l1": "SEG_FUNCTIONS",
            "l2": " ".join(str(x) for x in sed_vars),
        }
        for option in lines_ini:
            self.set_config("B7_sediment", option, lines_ini[option])

    def setup_climate_forcing(
        self,
        climate_fn: str,
        starttime: str,
        endtime: str,
        timestepsecs: int,
        climate_vars: list = ["temp", "temp_dew", "ssr", "wind10_u", "wind10_v", "tcc"],
        temp_correction: bool = False,
        dem_forcing_fn: str = None,
        **kwargs,
    ):
        """Setup Delwaq climate fluxes.

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
            List of climate_vars to consider. By default ["temp", "temp_dew", "ssr", "wind_u", "wind_v", "tcc"]
        temp_correction : bool, optional
             If True temperature are corrected using elevation lapse rate,
             by default False.
        dem_forcing_fn : str, default None
             Elevation data source with coverage of entire meteorological forcing domain.
             If temp_correction is True and dem_forcing_fn is provided this is used in
             combination with elevation at model resolution to correct the temperature.
        """
        self.logger.info(f"Setting dynamic data from climate source {climate_fn}.")
        # TODO function to get times
        freq = pd.to_timedelta(timestepsecs, unit="s")
        mask = self.hydromaps["modelmap"].values > 0

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
        if dem_forcing_fn != None:
            dem_forcing = self.data_catalog.get_rasterdataset(
                dem_forcing_fn,
                geom=ds.raster.box,  # clip dem with forcing bbox for full coverage
                buffer=2,
                variables=["elevtn"],
            ).squeeze()

        ds_out = xr.Dataset()

        # Start with wind
        if "wind10_u" in climate_vars and "wind10_v" in climate_vars:
            ds_out["wind"] = hydromt.workflows.forcing.wind(
                self.hydromaps,
                wind_u=ds["wind10_u"],
                wind_v=ds["wind10_v"],
                altitude_correction=False,
            )
            climate_vars.remove("wind10_u")
            climate_vars.remove("wind10_v")
        elif "wind" in climate_vars:
            ds_out["wind"] = hydromt.workflows.forcing.wind(
                self.hydromaps,
                wind=ds_out["wind"],
                altitude_correction=False,
            )
            climate_vars.remove("wind")
        # Add other variables
        temp_vars = [v for v in climate_vars if v.startswith("temp")]
        for v in climate_vars:
            if v in temp_vars:
                temp_v = hydromt.workflows.forcing.temp(
                    ds[v],
                    dem_model=self.hydromaps["elevtn"],
                    dem_forcing=dem_forcing,
                    lapse_correction=temp_correction,
                    logger=self.logger,
                    freq=freq,
                    reproj_method="nearest_index",
                    lapse_rate=-0.0065,
                    resample_kwargs=dict(label="right", closed="right"),
                )
                ds_out[v] = temp_v
            else:
                da_out = ds[v].raster.reproject_like(
                    self.hydromaps, method="nearest_index"
                )
                da_out = hydromt.workflows.forcing.resample_time(
                    da_out,
                    freq,
                    upsampling="bfill",  # we assume right labeled original data
                    downsampling="mean",
                    conserve_mass=False,
                    logger=self.logger,
                )
                ds_out[v] = da_out
        # Add basin mask
        ds_out.coords["mask"] = xr.DataArray(
            dims=ds_out.raster.dims,
            coords=ds_out.raster.coords,
            data=self.hydromaps["mask"].values,
            attrs=dict(_FillValue=self.hydromaps["mask"].raster.nodata),
        )
        # Add _FillValue to the data attributes
        for dvar in ds_out.data_vars.keys():
            nodata = ds_out[dvar].raster.nodata
            if nodata is not None:
                ds_out[dvar].attrs.update(_FillValue=nodata)
            else:
                ds_out[dvar].attrs.update(_FillValue=-9999.0)
            # Fill with zeros inside mask and keep NaN outside
            ds_out[dvar] = (
                ds_out[dvar].fillna(0).where(ds_out["mask"], ds_out[dvar].raster.nodata)
            )
        # Add zeros to the non surface water compartments
        comp_sfw = segments.sfwcomp(self.compartments, self.config)
        ds_out = segments.extend_comp_with_zeros(
            ds1c=ds_out, comp_ds1c=comp_sfw, compartments=self.compartments
        )

        # Rename and add to forcing
        rmdict = {k: v for k, v in self._FORCING.items() if k in ds_out.data_vars}
        ds_out = ds_out.rename(rmdict)
        self.set_forcing(ds_out)

        lines_ini = {
            "l1": "SEG_FUNCTIONS",
            "l2": " ".join(str(x) for x in ds_out.data_vars),
        }
        for option in lines_ini:
            self.set_config("B7_climate", option, lines_ini[option])

    def setup_emission_raster(
        self,
        emission_fn: str,
        scale_method: str = "average",
        fillna_method: str = "zero",
        fillna_value: int = 0.0,
        area_division: bool = False,
        comp_emi: str = "sfw",
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
            ds_like=self.staticmaps,
            method=scale_method,
            fillna_method=fillna_method,
            fillna_value=fillna_value,
            area_division=area_division,
            logger=self.logger,
        )
        # Attribute to comp_emi compartment and add zeros for the others
        ds_emi = segments.extend_comp_with_zeros(
            ds1c=ds_emi, comp_ds1c=comp_emi, compartments=self.compartments
        )
        self.set_staticmaps(ds_emi.rename(emission_fn))

    def setup_emission_vector(
        self,
        emission_fn: str,
        col2raster: str = "",
        rasterize_method: str = "value",
        comp_emi: str = "sfw",
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
                ds_like=self.staticmaps,
                col_name=col2raster,
                method=rasterize_method,
                mask_name="mask",
                logger=self.logger,
            )
        # Attribute to comp_emi compartment and add zeros for the others
        ds_emi = segments.extend_comp_with_zeros(
            ds1c=ds_emi, comp_ds1c=comp_emi, compartments=self.compartments
        )
        self.set_staticmaps(ds_emi.rename(emission_fn))

    def setup_emission_mapping(
        self,
        region_fn,
        mapping_fn=None,
        comp_emi="sfw",
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
            ds_like=self.staticmaps,
            source_name=region_fn,
            fn_map=mapping_fn,
            logger=self.logger,
        )
        rmdict = {k: v for k, v in self._MAPS.items() if k in ds_admin_maps.data_vars}
        # Attribute to comp_emi compartment and add zeros for the others
        ds_admin_maps = segments.extend_comp_with_zeros(
            ds1c=ds_admin_maps, comp_ds1c=comp_emi, compartments=self.compartments
        )
        self.set_staticmaps(ds_admin_maps.rename(rmdict))

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
        self._hydromaps.coords["mask"] = self._hydromaps["modelmap"].astype(bool)

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
            ds_out[dvar] = ds_out[dvar].where(ds_out["mask"], -9999.0)
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
                f"Writting dynamic data for timestep {i+1}/{len(timesteps)}"
            )
            # Flow
            flname = join(self.root, "dynamicdata", "flow.dat")
            flow_vars = self.get_config("B7_fluxes.l2").split(" ")
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
            vol_vars = self.compartments
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

    @property
    def pointer(self):
        """
        Dictionnary of schematisation attributes of a Delwaq model.

        Contains
        --------
        pointer: np.array
            Model pointer defining exchanges between segments
        compartments: list of str
            List of model compartments names
        boundaries: list of str
            List of model boundaries names
        fluxes: list of str
            List of model fluxes names
        nrofseg: int
            number of segments
        nrofexch: int
            number of exchanges
        nrofcomp: int
            number of compartments
        """
        if not self._pointer:
            # not implemented yet, fix later
            self._pointer = dict()
            # if self._read:
            #    self.read_pointer
        return self._pointer

    def set_pointer(self, attr, name):
        """Add model attribute property to pointer"""
        # Check that pointer attr is a four column np.array
        if name == "pointer":
            if not isinstance(attr, np.ndarray) and attr.shape[1] == 4:
                self.logger.warning(
                    "pointer values in self.pointer should be a np.ndarray with four columns."
                )
                return
        elif np.isin(name, ["compartments", "boundaries", "fluxes"]):
            if not isinstance(attr, list):
                self.logger.warning(
                    f"{name} object in self.pointer should be a list of names."
                )
                return
        elif np.isin(name, ["nrofseg", "nrofexch", "nrofcomp"]):
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
        elif "pointer" in self.pointer:
            nexch = self.pointer["pointer"].shape[0]
            self.set_pointer(nexch, "nrofexch")
        else:
            # from config
            nexch = self.get_config(
                "B4_nrofexch.l1",
                "0 0 0 ; x, y, z direction",
            )
            nexch = int(nexch.split(" ")[0])
            self.set_pointer(nexch, "nrofexch")
        return nexch

    @property
    def nrofcomp(self):
        """Fast accessor to nrofcomp property of pointer"""
        if "nrofcomp" in self.pointer:
            ncomp = self.pointer["nrofcomp"]
        elif "compartments" in self.pointer:
            ncomp = len(self.pointer["compartments"])
            self.set_pointer(ncomp, "nrofcomp")
        else:
            # from config
            ncomp = self.get_config(
                "B3_attributes.l7",
                "     1*01 ; sfw",
            )
            ncells = int(ncomp.split("*")[0])
            ncomp = int(self.nrofseg / ncells)
            self.set_pointer(ncomp, "nrofexch")
        return ncomp

    @property
    def compartments(self):
        """Fast accessor to compartments property of pointer"""
        if "compartments" in self.pointer:
            comp = self.pointer["compartments"]
        else:
            # from config
            nl = 7
            comp = []
            for i in range(self.nrofcomp):
                cp = self.get_config(
                    f"B3_attributes.l{nl}",
                    "     1*01 ; sfw",
                )
                cp = cp.split(";")[0][1:-1]
                comp.append(cp)
                nl += 1
            self.set_pointer(comp, "compartments")
        return comp

    @property
    def fluxes(self):
        """Fast accessor to fluxes property of pointer"""
        if "fluxes" in self.pointer:
            fl = self.pointer["fluxes"]
        else:
            # from config
            fl = self.get_config(
                "B7_fluxes.l2",
                "sfw>sfw inw>sfw",
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
        """Writes Delwaq netCDF geometry file (config/B3_waqgeom.nc)."""
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
        nosegh = int(np.max(np_ptid)) + 1  # because of the zero based
        # Get LDD map
        np_ldd = self.hydromaps["ldd"].squeeze(drop=True).values
        np_ldd[np_ldd == self.hydromaps["ldd"].raster.nodata] = 0

        # print("Input DataArray: ")
        # print(ptid.dims, ptid)
        ptid = xr.where(ptid == -1, np.nan, ptid)
        ptid = ptid.rename({"lat": "y", "lon": "x"})
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
        # uda_waqgeom.to_netcdf("updated_ugrid.nc")
        ##plot pointerId grid
        # da_ptid.ugrid.plot()
        # plt.show()
        ##CHECK resulting DataSet
        # print("CRS: ")
        # print(da_ptid.ugrid.crs)
        # print("Output DataSet: ")
        # print(uda_waqgeom)
        ##CHECK resulting NC file
        # ds = xu.open_dataset(fname)
        # uda = ds["ptid"]
        # uda.ugrid.plot()#uda.plot()
        # plt.show()
