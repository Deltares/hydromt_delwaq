"""Utilities for hydromt_delwaq."""
import logging
import os
from os.path import abspath, dirname, join
from pathlib import Path

import numpy as np
import xarray as xr
import xugrid as xu

logger = logging.getLogger(__name__)


__all__ = [
    "DATADIR",
    "dw_WriteSegmentOrExchangeData",
    "write_waqgeomfile",
]

DATADIR = join(dirname(abspath(__file__)), "data")


def dw_WriteSegmentOrExchangeData(
    ttime: int,
    fname: str | Path,
    datablock: np.ndarray,
    boundids: int,
    WriteAscii: bool = True,
    mode: str = "a",
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
    # Supress potential NaN values to avoid error (replaced by -999.0)
    if datablock.dtype == np.float64 or datablock.dtype == np.float32:
        datablock[np.isnan(datablock)] = -999.0
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


def write_waqgeomfile(hydromaps: xr.Dataset, fname: str | Path):
    """
    Derive and write WAQ geometry file from hydromaps dataset.

    Parameters
    ----------
    hydromaps : xr.Dataset
        Hydromaps dataset containing the geometry information.
        Required_variables: 'ptid', 'ldd'
    fname : str | Path
        Output file path for the WAQ geometry file.
    """
    # TODO: Update for several layers
    logger.info("Writing waqgeom.nc file")
    ptid = hydromaps["ptid"].copy()
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
    np_ldd = hydromaps["ldd"].squeeze(drop=True).values
    np_ldd[np_ldd == hydromaps["ldd"].raster.nodata] = 0

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
    da_ptid.ugrid.set_crs(crs=hydromaps.raster.crs)  # "EPSG:4326"
    uda_waqgeom = da_ptid.ugrid.to_dataset(
        optional_attributes=True
    )  # .to_netcdf("updated_ugrid.nc")
    uda_waqgeom.coords["projected_coordinate_system"] = -2147483647

    epsg_nb = int(hydromaps.raster.crs.to_epsg())
    uda_waqgeom["projected_coordinate_system"].attrs.update(
        dict(
            epsg=epsg_nb,
            grid_mapping_name="Unknown projected",
            longitude_of_prime_meridian=0.0,
            inverse_flattening=298.257223563,
            epsg_code=f"{hydromaps.raster.crs}",
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
    uda_waqgeom.to_netcdf(path=fname, mode="w")
