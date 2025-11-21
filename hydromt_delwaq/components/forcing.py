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

__all__ = ["DelwaqForcingComponent"]

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

        # Create output folder if it does not exist
        if not os.path.exists(dirname(join(self.root.path, filename))):
            os.makedirs(dirname(join(self.root.path, filename)))

        # Copy data to avoid modifying original
        ds_out = self.data.copy()

        # To avoid appending data to existing file, first delete all the .dat files
        dynDir = dirname(join(self.root.path, filename))
        if os.path.exists(dynDir):
            filelist = os.listdir(dynDir)
            for f in filelist:
                os.remove(os.path.join(dynDir, f))

        # Filter data with mask
        for dvar in ds_out.data_vars:
            # nodata = ds_out[dvar].raster.nodata
            # Change the novalue outside of mask for julia compatibility
            ds_out[dvar] = ds_out[dvar].where(ds_out["mask"], -9999.0)
            ds_out[dvar].attrs.update(_FillValue=-9999.0)

        logger.info("Writing dynamicmap files.")
        # Netcdf format
        if write_nc:
            fname = join(self.root.path, filename.format(name="dynamicdata"))
            fname = os.path.splitext(fname)[0] + ".nc"
            logger.info(f"Writing NetCDF copy of the forcing data to {fname}.")
            ds_out = ds_out.drop_vars(["mask", "spatial_ref"], errors="ignore")
            ds_out.to_netcdf(path=fname)

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
                # Flow
                flname = join(self.root.path, filename.format(name="flow"))
                flow_vars = self.model.fluxes
                flowblock = []
                for dvar in flow_vars:
                    # Maybe only clim or sed were updated
                    if dvar in ds_out.data_vars:
                        nodata = ds_out[dvar].raster.nodata
                        data = ds_out[dvar].isel(time=i).values.flatten()
                        data = data[data != nodata]
                        flowblock = np.append(flowblock, data)
                    else:
                        logger.info(f"Variable {dvar} not found in forcing data.")
                # Maybe do a check on len of flowback
                # based on number of variables,segments,timesteps
                if len(flowblock) > 0:
                    dw_WriteSegmentOrExchangeData(
                        timestepstamp[i], flname, flowblock, 1, WriteAscii=False
                    )

            # volume
            voname = join(self.root.path, filename.format(name="volume"))
            vol_vars = [self.model.surface_water]
            volblock = []
            for dvar in vol_vars:
                # Maybe only clim or sed were updated
                if dvar in ds_out.data_vars:
                    nodata = ds_out[dvar].raster.nodata
                    data = ds_out[dvar].isel(time=i).values.flatten()
                    data = data[data != nodata]
                    volblock = np.append(volblock, data)
            if len(volblock) > 0:
                dw_WriteSegmentOrExchangeData(
                    timestepstamp[i], voname, volblock, 1, WriteAscii=False
                )

            # sediment
            if (
                "B7_sediment" in self.model.config.data
                and i != len(ds_out.time.values) - 1
            ):
                sedname = join(self.root.path, filename.format(name="sediment"))
                sed_vars = self.model.config.get_value("B7_sediment.l2").split(" ")
                sedblock = []
                for dvar in sed_vars:
                    # sed maybe not updated or might be present for EM
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
            if (
                "B7_climate" in self.model.config.data
                and i != len(ds_out.time.values) - 1
            ):
                climname = join(self.root.path, filename.format(name="climate"))
                clim_vars = self.model.config.get_value("B7_climate.l2").split(" ")
                climblock = []
                for dvar in clim_vars:
                    # clim maybe not updated or might be present for EM
                    if dvar in ds_out.data_vars:
                        nodata = ds_out[dvar].raster.nodata
                        data = ds_out[dvar].isel(time=i).values.flatten()
                        data = data[data != nodata]
                        climblock = np.append(climblock, data)
                if len(climblock) > 0:
                    dw_WriteSegmentOrExchangeData(
                        timestepstamp[i], climname, climblock, 1, WriteAscii=False
                    )
