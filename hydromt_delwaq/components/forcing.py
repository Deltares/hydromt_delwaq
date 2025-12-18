"""Delwaq forcing component."""

import logging
import os
from os.path import dirname, join

import numpy as np
from hydromt import hydromt_step
from hydromt.model import Model
from hydromt.model.components import GridComponent
from tqdm import tqdm

from hydromt_delwaq.utils import dw_WriteSegmentOrExchangeData

__all__ = ["DelwaqForcingComponent", "DemissionForcingComponent"]

logger = logging.getLogger(f"hydromt.{__name__}")


class DelwaqForcingComponent(GridComponent):
    """Delwaq forcing component.

    This class is used for setting, creating, writing, and reading the Delwaq forcing.
    The forcing component data stored in the ``data`` property of this class is of the
    hydromt.gis.raster.RasterDataset type which is an extension of xarray.Dataset for
    regular grid.
    """

    def __init__(
        self,
        model: Model,
        *,
        filename: str = "dynamicdata/{name}.dat",
        region_component: str | None = None,
    ):
        """
        Initialize a DelwaqForcingComponent.

        Parameters
        ----------
        model : Model
            HydroMT model instance
        filename : str, optional
            Default path relative to the root where the forcing file(s) will be read
            and written. By default "dynamicdata/{name}.dat".
        region_component : str, optional
            The name of the region component to use as reference for this component's
            region. If None, the region will be set to the grid extent.
        """
        super().__init__(
            model=model,
            filename=filename,
            region_component=region_component,
            region_filename=None,
        )

    def read_forcing(self):
        """Read and forcing at <root/dynamicdata/> and parse to xr.Dataset."""
        self.root._assert_read_mode()
        self._initialize_grid(skip_read=True)

        # raise NotImplementedError()

    @hydromt_step
    def write(
        self,
        filename: str = "dynamicdata/{name}.dat",
        write_nc: bool = False,
    ):
        """
        Write forcing at <root/filename> in binary format.

        Can also write a NetCDF copy if write_nc is True.
        The output files are:

        * **flow.dat** binary: water fluxes [m3/s]
        * **volume.dat** binary: water volumes [m3]
        * **sediment.dat** binary: sediment particles from land erosion [g/timestep]
        * **climate.dat** binary: climate fluxes for climate_vars

        Parameters
        ----------
        filename : str, optional
            filename relative to model root and should contain a {name} placeholder,
            by default 'dynamicdata/{name}.dat'
        write_nc : bool, optional
            If True, write a NetCDF copy of the forcing data, by default False.
        """
        self.root._assert_write_mode()

        if len(self.data) == 0:
            logger.info("No forcing data found, skip writing.")
            return

        # Create and clean dynamicdata folder
        self._create_and_clean_folders(filename)
        # Mask data before writing
        ds_out = self._mask_forcing_before_write()
        # Netcdf format
        if write_nc:
            self._write_nc_copy(ds_out, filename)

        # Binary format
        timestepstamp = np.arange(
            0,
            (len(ds_out.time.values) + 1) * int(self.model.timestepsecs),
            self.model.timestepsecs,
        )

        for i in tqdm(
            np.arange(0, len(ds_out.time.values)), desc="Writing dynamic data"
        ):
            # if last timestep only write volume
            if i != len(ds_out.time.values) - 1:
                # flow
                self._write_flow_forcing(ds_out, filename, i, timestepstamp[i])

            # volume
            self._write_volume_forcing(ds_out, filename, i, timestepstamp[i])

            # sediment
            if (
                "B7_sediment" in self.model.config.data
                and i != len(ds_out.time.values) - 1
            ):
                self._write_sediment_forcing(ds_out, filename, i, timestepstamp[i])

            # climate
            if (
                "B7_climate" in self.model.config.data
                and i != len(ds_out.time.values) - 1
            ):
                self._write_climate_forcing(ds_out, filename, i, timestepstamp[i])

    def _create_and_clean_folders(self, filename):
        """Create and clean dynamicdata folder."""
        # Create output folder if it does not exist
        if not os.path.exists(dirname(join(self.root.path, filename))):
            os.makedirs(dirname(join(self.root.path, filename)))

        # To avoid appending data to existing file, first delete all the .dat files
        dynDir = dirname(join(self.root.path, filename))
        if os.path.exists(dynDir):
            filelist = os.listdir(dynDir)
            for f in filelist:
                os.remove(os.path.join(dynDir, f))

    def _mask_forcing_before_write(self):
        """Mask forcing data before writing to file."""
        # Copy data to avoid modifying original
        ds_out = self.data.copy()

        # Filter data with mask
        for dvar in ds_out.data_vars:
            # nodata = ds_out[dvar].raster.nodata
            # Change the novalue outside of mask for julia compatibility
            ds_out[dvar] = ds_out[dvar].where(ds_out["mask"], -9999.0)
            ds_out[dvar].attrs.update(_FillValue=-9999.0)

        return ds_out

    def _write_nc_copy(self, ds_out, filename):
        """Write NetCDF copy of the forcing data."""
        fname = join(self.root.path, filename.format(name="dynamicdata"))
        fname = os.path.splitext(fname)[0] + ".nc"
        logger.info(f"Writing NetCDF copy of the forcing data to {fname}.")
        ds_out = ds_out.drop_vars(["mask", "spatial_ref"], errors="ignore")
        ds_out.to_netcdf(path=fname)

    def _write_flow_forcing(self, ds_out, filename, timestamp, timestepstamp):
        """Write flow forcing data to binary file."""
        # Flow
        flname = join(self.root.path, filename.format(name="flow"))
        flow_vars = self.model.fluxes
        flowblock = []
        for dvar in flow_vars:
            # Maybe only clim or sed were updated
            if dvar in ds_out.data_vars:
                nodata = ds_out[dvar].raster.nodata
                data = ds_out[dvar].isel(time=timestamp).values.flatten()
                data = data[data != nodata]
                flowblock = np.append(flowblock, data)
            else:
                logger.info(f"Variable {dvar} not found in forcing data.")
        # Maybe do a check on len of flowback
        # based on number of variables,segments,timesteps
        if len(flowblock) > 0:
            dw_WriteSegmentOrExchangeData(
                timestepstamp, flname, flowblock, 1, WriteAscii=False
            )

    def _write_volume_forcing(self, ds_out, filename, timestamp, timestepstamp):
        """Write volume forcing data to binary file."""
        # volume
        voname = join(self.root.path, filename.format(name="volume"))
        vol_vars = [self.model.surface_water]
        volblock = []
        for dvar in vol_vars:
            # Maybe only clim or sed were updated
            if dvar in ds_out.data_vars:
                nodata = ds_out[dvar].raster.nodata
                data = ds_out[dvar].isel(time=timestamp).values.flatten()
                data = data[data != nodata]
                volblock = np.append(volblock, data)
        if len(volblock) > 0:
            dw_WriteSegmentOrExchangeData(
                timestepstamp, voname, volblock, 1, WriteAscii=False
            )

    def _write_sediment_forcing(self, ds_out, filename, timestamp, timestepstamp):
        """Write sediment forcing data to binary file."""
        sedname = join(self.root.path, filename.format(name="sediment"))
        sed_vars = self.model.config.get_value("B7_sediment.l2").split(" ")
        sedblock = []
        for dvar in sed_vars:
            # sed maybe not updated or might be present for EM
            if dvar in ds_out.data_vars:
                nodata = ds_out[dvar].raster.nodata
                data = ds_out[dvar].isel(time=timestamp).values.flatten()
                data = data[data != nodata]
                sedblock = np.append(sedblock, data)
        if len(sedblock) > 0:
            dw_WriteSegmentOrExchangeData(
                timestepstamp, sedname, sedblock, 1, WriteAscii=False
            )

    def _write_climate_forcing(self, ds_out, filename, timestamp, timestepstamp):
        """Write climate forcing data to binary file."""
        climname = join(self.root.path, filename.format(name="climate"))
        clim_vars = self.model.config.get_value("B7_climate.l2").split(" ")
        climblock = []
        for dvar in clim_vars:
            # clim maybe not updated or might be present for EM
            if dvar in ds_out.data_vars:
                nodata = ds_out[dvar].raster.nodata
                data = ds_out[dvar].isel(time=timestamp).values.flatten()
                data = data[data != nodata]
                climblock = np.append(climblock, data)
        if len(climblock) > 0:
            dw_WriteSegmentOrExchangeData(
                timestepstamp, climname, climblock, 1, WriteAscii=False
            )


class DemissionForcingComponent(GridComponent):
    """Demission forcing component.

    Class used for setting, creating, writing, and reading the Demission forcing.
    The forcing component data stored in the ``data`` property of this class is of the
    hydromt.gis.raster.RasterDataset type which is an extension of xarray.Dataset for
    regular grid.
    """

    @hydromt_step
    def write(
        self,
        filename: str = "dynamicdata/{name}.dat",
        write_nc: bool = False,
    ):
        """Write forcing at <root/fn> in binary format.

        Can also write a netcdf copy if ``write_nc`` is True.
        The output files are:

        * **hydrology.bin** binary file containing the hydrology data.
        * **sediment.dat** binary file containing the sediment data.
        * **climate.dat** binary file containing the climate data.

        Parameters
        ----------
        fn : str, optional
            filename relative to model root and should contain a {name} placeholder,
            by default 'dynamicdata/{name}.dat'
        write_nc : bool, optional
            If True, write a NetCDF copy of the forcing data, by default False.
        """
        self.root._assert_write_mode()

        if len(self.data) == 0:
            logger.info("No forcing data found, skip writing.")
            return

        # Create and clean dynamicdata folder
        self._create_and_clean_folders(filename)
        # Mask data before writing
        ds_out = self._mask_forcing_before_write()
        # Netcdf format
        if write_nc:
            self._write_nc_copy(ds_out, filename)

        # Binary format
        timestepstamp = np.arange(
            0,
            (len(ds_out.time.values) + 1) * int(self.model.timestepsecs),
            self.model.timestepsecs,
        )

        for i in tqdm(
            np.arange(0, len(ds_out.time.values)), desc="Writing dynamic data"
        ):
            # hydrology
            self._write_hydrology_forcing(ds_out, filename, i, timestepstamp[i])
            # climate
            if "B7_climate" in self.model.config.data:
                self._write_climate_forcing(ds_out, filename, i, timestepstamp[i])

    def _write_hydrology_forcing(self, ds_out, filename, timestamp, timestepstamp):
        """Write hydrology forcing data to binary file."""
        # hydrology
        fname = join(self.root.path, filename.format(name="hydrology"))
        fname = os.path.splitext(fname)[0] + ".bin"
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
                data = ds_out[dvar].isel(time=timestamp).values.flatten()
                data = data[data != nodata]
                datablock = np.append(datablock, data)
            else:
                logger.info(f"Variable {dvar} not found in forcing data.")
        if len(datablock) > 0:
            dw_WriteSegmentOrExchangeData(
                timestepstamp, fname, datablock, 1, WriteAscii=False
            )
        else:
            logger.info("No hydrology data found.")
