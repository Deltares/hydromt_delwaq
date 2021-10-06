# -*- coding: utf-8 -*-

import numpy as np
import xarray as xr
import geopandas as gpd
import logging

from hydromt import flw


logger = logging.getLogger(__name__)

__all__ = ["zonal_stats", "zonal_stats_grid"]


def zonal_stats(gdf, zones, variables=[], stats=[], method="overlay"):
    """Calculates zonal statisctics of vector samples aggregated for geometries.

    Adds new columns variables to the zones GeoDataFrame:
        - New columns: {variable}_{stat} (for variable in variables and stat in stats).

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        GeoDataFrame containing roads shape.
    zones : geopandas.GeoDataFrame
        Zones to compute the roads statistics.
    variables: list of str, optional.
        List of variables to derive statistics from.
        Either columns names from gdf (dtype must be either float or integer) or geometric properties such as ['length', 'area'].
        If geometric properties of the gdf are not available, they are computed on the spot.
    stats: list of str, callable
        Statistics to compute from raster values, options include
        {'count', 'min', 'max', 'sum', 'mean', 'std', 'median', 'q##'}.
         Multiple percentiles can be calculated using comma-seperated values, e.g.: 'q10,50,90'
         Statistics ignore the nodata value and are applied along the x and y dimension.
         By default ['mean'].
    method: str, optional
        Method for the spatial interpolation of gdf and zones. Either 'sjoin' for spatial join for
        gdf within zones, or 'overlay' for a union overlay of gdf with zones (default).


    Returns
    -------
    zones : geopandas.GeoDataFrame
        Zones with the computed zonal statistics.
    """
    # statistics
    _ST = ["count", "min", "max", "sum", "mean", "std", "median"]

    if isinstance(stats, str):
        stats = stats.split()
    elif callable(stats):
        stats = list([stats])

    if isinstance(variables, str):
        variables = variables.split()
    elif callable(stats):
        variables = list([variables])

    # zones and gdf should have the same crs
    if gdf.crs is not None and zones.crs is not None and gdf.crs != zones.crs:
        gdf = gdf.to_crs(zones.crs)

    # create a cell_id column for grouping
    zones["zone_id"] = np.arange(1, len(zones) + 1)

    # Spatial join of gdf and zones
    if method == "sjoin":
        gdf = gpd.sjoin(gdf, zones, how="inner", op="within")
    # Overlay with the zones
    elif method == "overlay":
        gdf = (
            gpd.overlay(gdf, zones, how="union", keep_geom_type=True)
            .explode()
            .reset_index(drop=True)
        )
    else:
        raise ValueError(
            f"Method {method} not valid, choose either 'sjoin' or 'overlay'."
        )

    # Loop over variables
    for var in variables:
        if var == "length" or var == "area":
            if gdf.crs.is_geographic:
                gdf = gdf.to_crs(3857)
            if var == "length":
                gdf["length"] = gdf.length
            else:
                gdf["area"] = gdf.area
            gdf = gdf.to_crs(zones.crs)
        elif var not in gdf.columns:
            raise ValueError(f"Variable {var} not found in gdf.")
        # Keep zone_id and variables columns only
        gdf_sub = gdf[["zone_id", var]]

        # Compute stats
        for stat in stats:
            if stat in _ST:
                # Group by using stats
                subset = getattr(gdf_sub.groupby(["zone_id"]), stat)()
            elif isinstance(stat, str) and stat.startswith("q"):
                qs = np.array([float(q) for q in stat.strip("q").split(",")])
                subset = getattr(gdf_sub.groupby(["zone_id"]), stat)(qs / 100)
            else:
                raise ValueError(f"Stat {stat} not valid.")
            # Rename subset columns
            subset = subset.rename(columns={var: f"{var}_{stat}"})
            # Add stat to zones
            zones = zones.merge(subset, on="zone_id")

    return zones.drop(["zone_id"], axis=1)


def zonal_stats_grid(
    gdf,
    ds_like,
    variables=[],
    stats=["mean"],
    mask_name="mask",
    method="overlay",
):
    """Calculates zonal statisctics of vector samples aggregated per raster cell.

    Returns a xr.Dataset with gridded road information :
        - Variables: {name_}{variable}_{stat} (name is taken from feature_filter, variable from variables and stat from stats).

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        GeoDataFrame containing roads shape.
    ds_like : xr.DataArray or xr.Dataset
        Dataset at model resolution.
    variables: list of str, optional.
        List of variables to derive statistics from.
        Either columns names from gdf (dtype must be either float or integer) or geometric properties such as ['length', 'area'].
        If geometric properties of the gdf are not available, they are computed on the spot.
    stats: list of str, callable
        Statistics to compute from raster values, options include
        {'count', 'min', 'max', 'sum', 'mean', 'std', 'median', 'q##'}.
         Multiple percentiles can be calculated using comma-seperated values, e.g.: 'q10,50,90'
         Statistics ignore the nodata value and are applied along the x and y dimension.
         By default ['mean'].
    mask_name : str
        Name of a basin mask array in ds_like.
    method: str, optional
        Method for the spatial interpolation of gdf and zones. Either 'sjoin' for spatial join for
        gdf within zones, or 'overlay' for a union overlay of gdf with zones (default).

    Returns
    -------
    ds_out : xarray.DataSet
        Dataset contains gridded zonal statistics at model resolution.
    """
    ds_out = None

    # Create model grid
    msktn = ds_like[mask_name]
    idx_valid = np.where(msktn.values.flatten() != msktn.raster.nodata)[0]
    gdf_grid = ds_like.raster.vector_grid().loc[idx_valid]

    gdf_stats = zonal_stats(
        gdf=gdf,
        zones=gdf_grid,
        variables=variables,
        stats=stats,
        method=method,
    )

    for var in variables:
        for stat in stats:
            da_out = ds_like.raster.rasterize(
                gdf_stats,
                col_name=f"{var}_{stat}",
                nodata=0,
                all_touched=False,
                dtype=None,
                sindex=False,
            )
            if ds_out is None:
                ds_out = da_out.to_dataset()
            else:
                dvar = f"{var}_{stat}"
                ds_out[dvar] = da_out

    return ds_out
