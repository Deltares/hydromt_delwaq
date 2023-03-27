# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import xarray as xr
import logging

from hydromt import gis_utils


logger = logging.getLogger(__name__)


__all__ = ["emission_raster", "admin", "emission_vector"]

RESAMPLING = {"EM_level1": "nearest", "EM_fact_X": "average"}
DTYPES = {"EM_level1": np.int16, "EM_fact_X": np.float32}


def gridarea(ds):
    """Returns a DataArray cointaining the area in m2 of the reference grid ds

    Parameters
    ----------
    da : xarray.DataArray or xarray.DataSet
        DataArray containing reference grid.

    Returns
    -------
    da_out : xarray.DataArray
        DataArray containing area in m2 of the reference grid.

    """

    realarea = gis_utils.reggrid_area(
        ds.raster.ycoords.values, ds.raster.xcoords.values
    )
    da_out = xr.DataArray(
        data=realarea.astype("float32"),
        coords=ds.raster.coords,
        dims=ds.raster.dims,
    )

    return da_out


def gridlength_gridwidth(ds):
    """Returns two DataArray for length and width of the grid"""
    lons = ds.raster.xcoords.values
    lats = ds.raster.ycoords.values
    xres = np.abs(np.mean(np.diff(lons)))
    yres = np.abs(np.mean(np.diff(lats)))
    ones = np.ones((lats.size, lons.size), dtype=lats.dtype)

    w, l = gis_utils.cellres(lats, xres, yres)
    width = w[:, None] * ones
    length = l[:, None] * ones

    da_len = xr.DataArray(
        data=length.astype("float32"),
        coords=ds.raster.coords,
        dims=ds.raster.dims,
    )
    da_len.raster.set_nodata(-9999.0)
    da_len.raster.set_crs(ds.raster.crs)
    da_len.attrs.update(unit="m")
    da_width = xr.DataArray(
        data=width.astype("float32"),
        coords=ds.raster.coords,
        dims=ds.raster.dims,
    )
    da_width.raster.set_nodata(-9999.0)
    da_width.raster.set_crs(ds.raster.crs)
    da_width.attrs.update(unit="m")

    return da_len.rename("length"), da_width.rename("width")


def emission_raster(
    da,
    ds_like,
    method="average",
    fillna_method="nearest",
    fillna_value=0.0,
    area_division=False,
    logger=logger,
):
    """Returns emission map.

    The following emission maps are calculated:\
    - da
    
    Parameters
    ----------
    da : xarray.DataArray
        DataArray containing emission map.
    ds_like : xarray.DataArray
        Dataset at model resolution.
    method : str {'average', 'nearest', 'mode'}
        Method for resampling.
    fillna_method : str {'nearest', 'zero', 'value'}
        Method to fill NaN values.
    fillna_value : float
        If fillna_method is set to 'value', NaNs in the emission maps will be replaced by this value.
    area_division : boolean
        If needed do the resampling in count/m2 (True) instead of count (False)

    Returns
    -------
    da_out : xarray.Dataset
        Dataset containing gridded emission map at model resolution.
    """
    nodata = da.raster.nodata
    if nodata is not None:
        if fillna_method == "zero":
            da = da.where(da.values != nodata).fillna(0)
        elif fillna_method == "value":
            da = da.where(da.values != nodata).fillna(fillna_value)
        else:
            da = da.raster.interpolate_na(method="nearest")
    else:
        if np.issubdtype(da.dtype, np.signedinteger):
            nodata = -999
        elif np.issubdtype(da.dtype, np.unsignedinteger):
            nodata = 255
        else:
            nodata = -999.0
        da.raster.set_nodata(nodata)

    if area_division:
        da_area = gridarea(da)
        da = da / da_area

    logger.info(f"Deriving {da.name} using {method} resampling (nodata={nodata}).")
    # da = da.astype(np.float32)

    da_out = da.raster.reproject_like(ds_like, method=method)
    if area_division:
        da_area = gridarea(da_out)
        da_out = da_out * da_area
    da_out.attrs.update(_FillValue=nodata)

    return da_out


def emission_vector(
    gdf,
    ds_like,
    col_name="",
    method="value",
    mask_name="basmsk",
    logger=logger,
):
    """Returns gridded emission data from vector.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        GeoDataFrame containing emission data.
    ds_like : xarray.DataArray
        Dataset at model resolution.
    method : str {'value', 'fraction'}
        Method for rasterizing.
    mask_name : str
        Name of a basin mask array in ds_like.

    Returns
    -------
    da_out : xarray.Dataarray
        Dataarray containing gridded emission map at model resolution.
    """
    if method == "value":
        da_out = ds_like.raster.rasterize(
            gdf,
            col_name=col_name,
            nodata=0,
            all_touched=True,
            dtype=None,
            sindex=False,
        )
    elif method == "fraction":
        # Using the vectors directly (long computation time but most accurate)
        # Create vector grid (for calculating fraction and storage per grid cell)
        logger.debug(
            "Creating vector grid for calculating coverage fraction per grid cell"
        )
        gdf["geometry"] = gdf.geometry.buffer(0)  # fix potential geometry errors
        msktn = ds_like[mask_name]
        idx_valid = np.where(msktn.values.flatten() != msktn.raster.nodata)[0]
        gdf_grid = ds_like.raster.vector_grid().loc[idx_valid]
        gdf_grid["coverfrac"] = np.zeros(len(idx_valid))
        gdf_grid["area"] = gdf_grid.to_crs(
            3857
        ).area  # area calculation in projected crs

        # Calculate fraction per (vector) grid cell
        # Looping over each vector shape
        for i in range(len(gdf)):
            shape = gdf.iloc[i]
            gridded_shape = gdf_grid.intersection(shape.geometry)
            gridded_shape = gridded_shape.loc[~gridded_shape.is_empty]
            idxs = gridded_shape.index
            if np.any(idxs):
                # area calculation needs projected crs
                sharea_cell = gridded_shape.to_crs(3857).area
                gdf_grid.loc[idxs, "coverfrac"] += (
                    sharea_cell / gdf_grid.loc[idxs, "area"]
                )
        # Create the rasterized coverage fraction map
        da_out = ds_like.raster.rasterize(
            gdf_grid,
            col_name="coverfrac",
            nodata=0,
            all_touched=False,
            dtype=None,
            sindex=False,
        )

    return da_out


def admin(da, ds_like, source_name, fn_map, logger=logger):
    """Returns administrative boundaries map and related parameter maps.
    The parameter maps are prepared based on administrative boundaries and
    mapping table as provided in the generic data folder of hydromt.

    The following maps are calculated:\
    - TODO
    
    Parameters
    ----------
    da : xarray.DataArray
        DataArray containing Admin classes.
    ds_like : xarray.DataArray
        Dataset at model resolution.

    Returns
    -------
    ds_out : xarray.Dataset
        Dataset containing gridded emission factor based maps
    """

    # read csv with remapping values
    df = pd.read_csv(fn_map, index_col=0, sep="[,;]", engine="python", dtype=DTYPES)
    # limit dtypes to avoid gdal errors downstream
    ddict = {"float64": np.float32, "int64": np.int32}
    dtypes = {c: ddict.get(str(df[c].dtype), df[c].dtype) for c in df.columns}
    df = pd.read_csv(fn_map, index_col=0, sep="[,;]", engine="python", dtype=dtypes)
    keys = df.index.values
    # include index_col as parameter (for output to map)
    df["EM_ID"] = df.index
    df.EM_ID = df.EM_ID.astype(np.int32)
    # define list of parameters for which output maps are created
    params = [p for p in df.columns if p.startswith("EM_")]
    # setup ds out
    ds_out = xr.Dataset(coords=ds_like.raster.coords)

    # setup reclass method
    def reclass(x):
        return np.vectorize(d.get)(x, nodata)

    # apply for each parameter
    for param in params:
        # TODO change average into value with highest occurence
        method = RESAMPLING.get(param, "average")
        values = df[param].values
        nodata = values[-1]  # NOTE values is set in last row
        d = dict(zip(keys, values))  # NOTE global param in reclass method
        da = da.raster.interpolate_na(method="nearest")
        logger.info(f"Deriving {param} using {method} resampling (nodata={nodata}).")
        da_param = xr.apply_ufunc(
            reclass, da, dask="parallelized", output_dtypes=[values.dtype]
        )
        da_param.attrs.update(_FillValue=nodata)  # first set new nodata values
        if (
            param == "EM_ID"
        ):  # rename ID map to source_name (to avoid overwriting when several sources use this admin mapping function)
            ds_out[source_name] = da_param.raster.reproject_like(
                ds_like, method=method
            )  # then resample
        else:
            ds_out[param] = da_param.raster.reproject_like(
                ds_like, method=method
            )  # then resample
        # ds_out = ds_out.rename({"EM_ID":source_name})

    return ds_out
