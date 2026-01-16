"""Workflows for monitoring points and areas."""

import logging

import geopandas as gpd
import numpy as np
import xarray as xr

logger = logging.getLogger(__name__)

__all__ = [
    "monitoring_points_from_dataarray",
    "monitoring_points_from_geodataframe",
    "monitoring_areas",
]

mv = -999  # missing value


def _get_nb_points(monpoints: xr.DataArray) -> int:
    """Get number of monitoring points."""
    points = monpoints.values.flatten()
    points = points[points != mv]

    return len(points)


def monitoring_points_from_dataarray(
    monpoints: xr.DataArray,
) -> tuple[int, xr.DataArray]:
    """
    Prepare Delwaq monitoring points.

    Parameters
    ----------
    monpoints : {'segments', 'data_source', 'path to station location'}
        Either source from DataCatalog, path to a station location dataset
        or if segments, all segments are monitored.
    """
    if not isinstance(monpoints, xr.DataArray):
        raise ValueError("monpoints must be an xarray DataArray.")

    return _get_nb_points(monpoints), monpoints


def monitoring_points_from_geodataframe(
    monpoints_gdf: gpd.GeoDataFrame,
    ds_like: xr.Dataset,
) -> tuple[int, xr.DataArray, gpd.GeoDataFrame | None]:
    """
    Prepare Delwaq monitoring points from a GeoDataFrame.

    Parameters
    ----------
    monpoints_gdf : gpd.GeoDataFrame
        GeoDataFrame with monitoring point locations.
    ds_like : xr.Dataset
        Dataset to use as reference for rasterizing.
    """
    if monpoints_gdf.index.size == 0:
        logger.warning("No monitoring points found within domain")
        monpoints_gdf = None
    else:
        monpoints_gdf = monpoints_gdf.to_crs(ds_like.raster.crs)
        monpoints_gdf.index.name = "index"
        monpoints = ds_like.raster.rasterize(monpoints_gdf, col_name="index", nodata=mv)

    return _get_nb_points(monpoints), monpoints, monpoints_gdf


def monitoring_areas(
    area_name: str,
    ds: xr.Dataset,
) -> tuple[int, xr.DataArray]:
    """
    Prepare Delwaq monitoring areas.

    Parameters
    ----------
    area_name : {'subcatch', 'riverland'}
        Type of monitoring areas to create.
    ds : xr.Dataset
        Dataset to use as reference for creating monitoring areas.
    """
    if area_name == "subcatch":  # subcatch
        basins = xr.where(
            ds["basins"] > 0,
            ds["basins"],
            mv,
        ).astype(np.int32)
        # Number or monitoring areas
        areas = basins.values.flatten()
        areas = areas[areas != mv]
        nb_areas = len(np.unique(areas))
        monareas = basins
    elif area_name == "riverland":  # riverland
        # seperate areas for land cells (1) and river cells (2)
        lr_areas = xr.where(
            ds["river"],
            2,
            xr.where(ds["basins"], 1, mv),
        ).astype(np.int32)
        # Apply the current model mask
        lr_areas = lr_areas.where(ds["mask"], mv)
        lr_areas.raster.set_nodata(mv)
        # Number or monitoring areas
        areas = lr_areas.values.flatten()
        areas = areas[areas != mv]
        nb_areas = len(np.unique(areas))
        monareas = lr_areas
    else:
        raise ValueError(
            f"Unknown monitoring area type {area_name}. "
            "Valid options are 'subcatch' or 'riverland'."
        )

    # Update attributes
    monareas.attrs.update(_FillValue=mv)
    monareas.attrs["mon_areas"] = area_name

    return nb_areas, monareas
