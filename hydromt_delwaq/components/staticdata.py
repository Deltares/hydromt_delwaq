"""Custom Delwaq staticdata component module."""

import logging
import os
from os.path import dirname, join, splitext

import numpy as np
from hydromt import hydromt_step
from hydromt.model import Model
from hydromt.model.components import GridComponent

from hydromt_delwaq.utils import dw_WriteSegmentOrExchangeData

__all__ = ["DelwaqStaticdataComponent"]

logger = logging.getLogger(f"hydromt.{__name__}")


class DelwaqStaticdataComponent(GridComponent):
    """Delwaq staticdata component.

    Inherits from the HydroMT-core GridComponent model-component.
    It is used for setting, creating, writing, and reading static and cyclic data for a
    Delwaq model on a regular grid. The component data, stored in the ``data``
    property of this class, is of the hydromt.gis.raster.RasterDataset type which
    is an extension of xarray.Dataset for regular grid.
    """

    def __init__(
        self,
        model: Model,
        *,
        filename: str = "staticdata/{name}.dat",
        region_filename: str = "geoms/region.geojson",
    ):
        """Initialize a DelwaqStaticdataComponent.

        Parameters
        ----------
        model : Model
            HydroMT model instance
        filename : str
            The path to use for reading and writing of component data by default.
            By default "staticdata/{name}.dat".
        region_filename : str
            The path to use for reading and writing of the region data by default.
            By default "geoms/region.geojson".
        """
        super().__init__(
            model,
            filename=filename,
            region_component=None,
            region_filename=region_filename,
        )

    ## I/O methods
    @hydromt_step
    def read(self, filename: str = "staticdata/{name}.dat", **kwargs):
        """
        Read staticdata at <root/fn> and parse to xarray.

        For now, this function only allows to read the netcdf copy of the grid.

        Parameters
        ----------
        filename : str, optional
            filename relative to model root and should contain a {name} placeholder,
            by default 'staticdata/{name}.dat'. For the netcdf file, the placeholder
            will be replaced by the name staticdata.
        """
        fname = join(self.root.path, filename.format(name="staticdata"))
        fname = splitext(fname)[0] + ".nc"
        super().read(fname, **kwargs)

    @hydromt_step
    def write(self, filename: str = "staticdata/{name}.dat"):
        """
        Write grid at <root/fn> in NetCDF and binary format.

        Parameters
        ----------
        fn : str, optional
            filename relative to model root and should contain a {name} placeholder,
            by default 'staticdata/{name}.dat'. For the netcdf file, the placeholder
            will be replaced by the name staticdata.
        """
        self.root._assert_write_mode()

        if len(self.data) == 0:
            logger.info(
                f"{self.model.name}.{self.name_in_model}: "
                "No grid data found, skip writing."
            )
            return

        # Create output folder if it does not exist
        if not os.path.exists(dirname(join(self.root.path, filename))):
            os.makedirs(dirname(join(self.root.path, filename)))

        ds_out = self.data
        # Filter data with mask
        for dvar in ds_out.data_vars:
            ds_out[dvar] = ds_out[dvar].raster.mask(mask=ds_out["mask"])

        logger.info("Writing staticmap files.")
        # Netcdf format
        fname = join(self.root.path, filename.format(name="staticdata"))
        fname = splitext(fname)[0] + ".nc"
        # Update attributes for gdal compliance
        # ds_out = ds_out.raster.gdal_compliant(rename_dims=False)
        # Not ideal but to avoid hanging issues in r+ mode
        if self.root.is_reading_mode() and self.root.is_writing_mode():
            ds_out.load()
        ds_out.to_netcdf(path=fname)

        # Binary format
        mask = ds_out["mask"].values.flatten()
        for dvar in ds_out.data_vars:
            if dvar == "monpoints" or dvar == "monareas":
                continue
            # Check if data is 3D
            if len(ds_out[dvar].shape) == 3:
                dim0 = ds_out[dvar].raster.dim0
                for i in range(len(ds_out[dim0])):
                    dim0_val = ds_out[dim0][i].item()
                    fname = join(
                        self.root.path,
                        filename.format(name=f"{dvar}_{dim0}_{dim0_val}"),
                    )
                    data = ds_out[dvar].sel({dim0: dim0_val}).values.flatten()
                    data = data[mask]
                    dw_WriteSegmentOrExchangeData(
                        0, fname, data, 1, WriteAscii=False, mode="w"
                    )
            else:
                fname = join(self.root.path, filename.format(name=dvar))
                data = ds_out[dvar].values.flatten()
                data = data[mask]
                dw_WriteSegmentOrExchangeData(
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

    def write_monitoring(self, monpoints, monareas):
        """
        Write monitoring files and config in ASCII format.

        Input:
            - monpoints - xr.DataArray of monitoring points location
            - monareas - xr.DataArray of monitoring areas location
        """
        if not os.path.isdir(join(self.root.path, "config")):
            os.makedirs(join(self.root.path, "config"))

        ptid = self.model.hydromaps.data["ptid"].values.flatten()
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
                fname = join(self.root.path, "config", "B2_" + name + ".inc")
                with open(fname, "w") as exfile:
                    print(";Written by hydroMT", file=exfile)
                    if name == "stations":
                        np.savetxt(exfile, stations, fmt="%.20s")
                    else:
                        np.savetxt(exfile, stations_balance, fmt="%.20s")
        else:
            fname = join(self.root.path, "config", "B2_stations.inc")
            with open(fname, "w") as exfile:
                print(
                    ";Written by hydroMT: no monitoring points were set.", file=exfile
                )

        # Monitoring areas
        if monareas is not None:
            mv = monareas.raster.nodata
            areas = monareas.values.flatten()
            id_areas = ptid[areas != mv]
            areas = areas[areas != mv]
            # Write to file
            fname = join(self.root.path, "config", "B2_monareas.inc")
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
            fname = join(self.root.path, "config", "B2_monareas.inc")
            exfile = open(fname, "w")
            print(";Written by hydroMT: no monitoring areas were set.", file=exfile)
            exfile.close()
