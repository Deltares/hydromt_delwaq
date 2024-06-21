"""Implement demission model class."""

import logging
import os
from os.path import join
from pathlib import Path
from typing import Dict, List, Union

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from hydromt import workflows
from tqdm import tqdm

from . import DATADIR
from .delwaq import DelwaqModel
from .utils import dw_WriteSegmentOrExchangeData
from .workflows import config, forcing, geometry, roads, segments

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
    """Demission model class."""

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
        root: Union[str, Path] = None,
        mode: str = "w",
        config_fn: Union[str, Path] = None,
        hydromodel_name: str = "wflow",
        hydromodel_root: Union[str, Path] = None,
        data_libs: List[Union[str, Path]] = None,
        logger=logger,
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
        hydromodel_name : str, optional
            Name of the hydromodel used to build the emission model. Only useful in
            update mode as this is taken from the ``region`` argument in
            **setup_basemaps** method. By default "wflow".
        hydromodel_root : str or Path, optional
            Root of the hydromodel used to build the emission model. Only useful in
            update mode as this is taken from the ``region`` argument in
            **setup_basemaps** method. By default None.
        data_libs : List[str, Path], optional
            List of data catalog configuration files, by default None
        logger:
            The logger to be used.
        """
        super().__init__(
            root=root,
            mode=mode,
            config_fn=config_fn,
            hydromodel_name=hydromodel_name,
            hydromodel_root=hydromodel_root,
            data_libs=data_libs,
            logger=logger,
        )

        # d-emission specific
        self._geometry = None

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

        self.logger.info("Preparing EM basemaps from hydromodel.")

        ### Select and build hydromaps from model ###
        # Initialise hydromaps
        ds_hydro = segments.hydromaps(hydromodel)
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
            hydromodel, maps=maps, logger=self.logger
        )
        ds_stat["ptiddown"] = self.hydromaps["ptiddown"].squeeze(drop=True)

        rmdict = {k: v for k, v in self._MAPS.items() if k in ds_stat.data_vars}
        self.set_grid(ds_stat.rename(rmdict))
        self.grid.coords["mask"] = self.hydromaps["mask"]

        ### Geometry data ###
        geometry_data = geometry.compute_geometry(
            ds=hydromodel.grid, mask=self.grid["mask"]
        )
        self._geometry = pd.DataFrame(
            geometry_data, columns=(["TotArea", "fPaved", "fUnpaved", "fOpenWater"])
        )

        ### Config ###
        configs = config.base_config(
            nrofseg=nrofseg,
            nrofexch=0,
            layer_name="EM",
            add_surface=True,
        )
        for file in configs:
            for option in configs[file]:
                self.set_config(file, option, configs[file][option])

    def setup_roads(
        self,
        roads_fn: Union[str, gpd.GeoDataFrame],
        highway_list: Union[str, List[str]],
        country_list: Union[str, List[str]],
        non_highway_list: Union[str, List[str]] = None,
        country_fn: Union[str, gpd.GeoDataFrame] = None,
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
            List of countries for the model area in country_fn and optionnally in
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
            gdf_country = self.data_catalog.get_geodataframe(
                country_fn, dst_crs=self.crs
            )
            gdf_country = gdf_country.astype({"country_code": "str"})
            gdf_country = gdf_country.iloc[
                np.isin(gdf_country["country_code"], country_list)
            ]
            # Read the roads data and mask with country geom
            gdf_roads = self.data_catalog.get_geodataframe(
                roads_fn, dst_crs=self.crs, geom=gdf_country, variables=["road_type"]
            )
            # Compute country statistics
            ds_country = roads.roads_emissions_country(
                gdf_roads=gdf_roads,
                highway_list=highway_list,
                non_highway_list=non_highway_list,
                logger=self.logger,
            )
            # Add to grid
            self.set_grid(ds_country)

        ### Road statistics per segment ###
        # Filter road gdf with model mask instead of whole country
        mask = self.grid["mask"].astype("int32").raster.vectorize()
        gdf_roads = self.data_catalog.get_geodataframe(
            roads_fn, dst_crs=self.crs, geom=mask, variables=["road_type"]
        )
        # Compute segment statistics
        ds_segments = roads.roads_emissions_segment(
            gdf_roads=gdf_roads,
            highway_list=highway_list,
            non_highway_list=non_highway_list,
            logger=self.logger,
        )
        self.set_grid(ds_segments)

    def setup_hydrology_forcing(
        self,
        hydro_forcing_fn: Union[str, xr.Dataset],
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

        # Update model timestepsecs attribute
        self.timestepsecs = timestepsecs

        # Compute hydrology forcing
        ds = forcing.hydrology_forcing_em(
            ds=ds_in,
            ds_model=self.hydromaps,
            timestepsecs=timestepsecs,
            include_transport=include_transport,
            logger=self.logger,
        )
        # Rename xdim and ydim
        if ds.raster.x_dim != self.grid.raster.x_dim:
            ds = ds.rename({ds.raster.x_dim: self.grid.raster.x_dim})
        if ds.raster.y_dim != self.grid.raster.y_dim:
            ds = ds.rename({ds.raster.y_dim: self.grid.raster.y_dim})

        # Add to forcing
        self.set_forcing(ds)

        # Add timers info to config
        time_config = config.time_config(
            starttime=starttime,
            endtime=endtime,
            timestepsecs=timestepsecs,
            sysclock_format="days",
        )
        for file in time_config:
            for option in time_config[file]:
                self.set_config(file, option, time_config[file][option])

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
            self.set_config("B7_hydrology", option, lines_ini[option])

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
        self.write_geometry()
        self.write_forcing()

    def read_config(
        self,
        skip: List[str] = [
            "B7_geometry",
            "B2_stations",
            "B2_stations-balance",
            "B2_monareas",
        ],
    ):
        """Read config files in ASCII format at <root/config>."""
        # Skip geometry file (should be read with read_geometry())
        # Skip monitoring files (should be read with read_monitoring())
        super().read_config(skip=skip)

    @property
    def geometry(self):
        """
        Pandas DataFrame containing the geometry of the EM compartment.

        Columns of the DataFrame are:
        * TotArea: area of the segment
        * fPaved: fraction of the segment that is paved
        * fUnpaved: fraction of the segment that is unpaved
        * fOpenWater: fraction of the segment that is open water
        """
        # if not self._geometry:
        #     # not implemented yet, fix later
        #     self._geometry = pd.DataFrame()
        #     if self._read:
        #         self.read_geometry()
        return self._geometry

    def read_geometry(self):
        """Read Delwaq EM geometry file."""
        raise NotImplementedError()

    def write_geometry(self):
        """Write geometry at <root/staticdata> in ASCII and binary format."""
        self._assert_write_mode()
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

    def write_forcing(self, write_nc: bool = False):
        """Write forcing at <root/dynamicdata> in binary format.

        Can also write a netcdf copy if ``write_nc`` is True.
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
        for i in tqdm(timesteps, desc="Writing dynamic data"):
            # hydrology
            fname = join(self.root, "dynamicdata", "hydrology.bin")
            datablock = []
            dvars = ["precip", "runPav", "runUnp", "infilt"]
            if "totflw" in ds_out.data_vars:
                dvars.extend(["totflw"])
            else:
                dvars.extend(["exfilt", "vwcproot", "q_land", "q_ss"])
            for dvar in dvars:
                # Maybe only clim or sed were updated
                if dvar in ds_out.data_vars:
                    nodata = ds_out[dvar].raster.nodata
                    data = ds_out[dvar].isel(time=i).values.flatten()
                    data = data[data != nodata]
                    datablock = np.append(datablock, data)
                else:
                    self.logger.info(f"Variable {dvar} not found in forcing data.")
            if len(datablock) > 0:
                dw_WriteSegmentOrExchangeData(
                    timestepstamp[i], fname, datablock, 1, WriteAscii=False
                )
            else:
                self.logger.info("No hydrology data found.")
            # sediment
            if "B7_sediment" in self.config:
                sedname = join(self.root, "dynamicdata", "sediment.dat")
                sed_vars = self.get_config("B7_sediment.l2").split(" ")
                sedblock = []
                for dvar in sed_vars:
                    # sed maybe not updated or might be present for WQ
                    if dvar in ds_out.data_vars:
                        nodata = ds_out[dvar].raster.nodata
                        data = ds_out[dvar].isel(time=i).values.flatten()
                        data = data[data != nodata]
                        sedblock = np.append(sedblock, data)
                if len(sedblock) > 0:
                    dw_WriteSegmentOrExchangeData(
                        timestepstamp[i], sedname, sedblock, 1, WriteAscii=False
                    )
            # climate
            if "B7_climate" in self.config:
                climname = join(self.root, "dynamicdata", "climate.dat")
                clim_vars = self.get_config("B7_climate.l2").split(" ")
                climblock = []
                for dvar in clim_vars:
                    # clim maybe not updated or might be present for WQ
                    if dvar in ds_out.data_vars:
                        nodata = ds_out[dvar].raster.nodata
                        data = ds_out[dvar].isel(time=i).values.flatten()
                        data = data[data != nodata]
                        climblock = np.append(climblock, data)
                if len(climblock) > 0:
                    dw_WriteSegmentOrExchangeData(
                        timestepstamp[i], climname, climblock, 1, WriteAscii=False
                    )

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
        else:
            nexch = self.nrofseg * 5
            self.set_pointer(nexch, "nrofexch")
        return nexch

    @property
    def fluxes(self):
        """Fast accessor to fluxes property of pointer."""
        if "fluxes" in self.pointer:
            fl = self.pointer["fluxes"]
        else:
            # from config
            fl = self.get_config(
                "B7_hydrology.l2",
                fallback=(
                    "Rainfall RunoffPav RunoffUnp Infiltr "
                    + "Exfiltr Overland Subsurface"
                ),
            )
            fl = fl.split(" ")
            self.set_pointer(fl, "fluxes")
        return fl
