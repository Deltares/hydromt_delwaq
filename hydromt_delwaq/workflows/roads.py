"""Workflow for roads data."""

import logging
from typing import Dict, List

import geopandas as gpd
import numpy as np
import xarray as xr

from .emissions import emission_vector

logger = logging.getLogger(__name__)

__all__ = [
    "zonal_stats",
    "zonal_stats_grid",
    "roads_emissions_country",
    "roads_emissions_segments",
]


def zonal_stats(
    gdf: gpd.GeoDataFrame,
    zones: gpd.GeoDataFrame,
    variables: List[str] = [],
    stats: List[str] = [],
    method: str = "overlay",
) -> gpd.GeoDataFrame:
    """Calculate zonal statisctics of vector samples aggregated for geometries.

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
        Either columns names from gdf (dtype must be either float or integer) or
        geometric properties such as ['length', 'area'].
        If geometric properties of the gdf are not available, they are computed on
        the spot.
    stats: list of str, callable
        Statistics to compute from raster values, options include
        {'count', 'min', 'max', 'sum', 'mean', 'std', 'median', 'q##'}.
        Multiple percentiles can be calculated using comma-seperated values,
        e.g.: 'q10,50,90'
        Statistics ignore the nodata value and are applied along the x and y dimension.
        By default ['mean'].
    method: str, optional
        Method for the spatial interpolation of gdf and zones. Either 'sjoin' for
        spatial join for gdf within zones, or 'overlay' for a union overlay of gdf with
        zones (default).


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
        gdf = gpd.sjoin(gdf, zones, how="inner", predicate="within")
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
) -> xr.Dataset:
    """Calculate zonal statisctics of vector samples aggregated per raster cell.

    Returns a xr.Dataset with gridded road information :
        - Variables: {name_}{variable}_{stat} (name is taken from feature_filter,
        variable from variables and stat from stats).

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        GeoDataFrame containing roads shape.
    ds_like : xr.DataArray or xr.Dataset
        Dataset at model resolution.
    variables: list of str, optional.
        List of variables to derive statistics from.
        Either columns names from gdf (dtype must be either float or integer) or
        geometric properties such as ['length', 'area'].
        If geometric properties of the gdf are not available, they are computed on
        the spot.
    stats: list of str, callable
        Statistics to compute from raster values, options include
        {'count', 'min', 'max', 'sum', 'mean', 'std', 'median', 'q##'}.
         Multiple percentiles can be calculated using comma-seperated values,
         e.g.: 'q10,50,90'
         Statistics ignore the nodata value and are applied along the x and y dimension.
         By default ['mean'].
    mask_name : str
        Name of a basin mask array in ds_like.
    method: str, optional
        Method for the spatial interpolation of gdf and zones. Either 'sjoin' for
        spatial join for gdf within zones, or 'overlay' for a union overlay of gdf with
        zones (default).

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


def _preprocess_roads_emissions(
    gdf_roads: gpd.GeoDataFrame,
    highway_list: List[str],
    non_highway_list: List[str] = None,
) -> Dict:
    """
    Preprocess arguments for roads_emissions_country and roads_emissions_segments.

    Parameters
    ----------
    gdf_roads : geopandas.GeoDataFrame
        GeoDataFrame containing roads shape.

        * Required columns: ``road_type``.
    highway_list : list of str
        List of highway road types.
    non_highway_list : list of str, optional
        List of non-highway road types. If None, all road types not in highway_list
        are considered as non-highway.

    Returns
    -------
    feature_filter : dict
        Dictionary with keys 'hwy' and 'nnhwy' for highway and non-highway roads
        respectively. Each key contains a dictionary with the column name and the
        list of road types.
    """
    # Convert string to lists
    if not isinstance(highway_list, list):
        highway_list = [highway_list]

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

    return feature_filter


def roads_emissions_country(
    gdf_roads: gpd.GeoDataFrame,
    gdf_country: gpd.GeoDataFrame,
    ds_like: xr.Dataset,
    highway_list: List[str],
    non_highway_list: List[str] = None,
) -> xr.Dataset:
    """
    Compute roads statistics per country of interest.

    Parameters
    ----------
    gdf_roads : geopandas.GeoDataFrame
        GeoDataFrame containing roads shape.

        * Required columns: ``road_type``.
    gdf_country : geopandas.GeoDataFrame
        GeoDataFrame containing country shape.
    ds_like : xarray.Dataset
        Dataset at model resolution.

        * Required variables: ``mask``.
    highway_list : list of str
        List of highway road types.
    non_highway_list : list of str, optional
        List of non-highway road types. If None, all road types not in highway_list
        are considered as non-highway.

    Returns
    -------
    ds_country : xarray.Dataset
        Dataset contains roads statistics per country of interest.
        Contains the following variables:

        * ``hwy_length_sum_country``: highway road length per country of interest.
        * ``nnhwy_length_sum_country``: non-highway road length per country of interest.
    """
    # To ds_like crs
    gdf_roads = gdf_roads.to_crs(ds_like.raster.crs)
    gdf_country = gdf_country.to_crs(ds_like.raster.crs)
    # Common preprocessing of arguments with segments
    feature_filter = _preprocess_roads_emissions(
        gdf_roads=gdf_roads,
        highway_list=highway_list,
        non_highway_list=non_highway_list,
    )
    # Loop over feature_filter
    ds_country = xr.Dataset()
    for name, colfilter in feature_filter.items():
        logger.info(f"Computing {name} roads statistics per country of interest")
        # Filter gdf
        colname = [k for k in colfilter.keys()][0]
        subset_roads = gdf_roads.iloc[
            np.isin(gdf_roads[colname], colfilter.get(colname))
        ]
        gdf_country = zonal_stats(
            gdf=subset_roads,
            zones=gdf_country,
            variables=["length"],
            stats=["sum"],
            method="sjoin",
        )
        gdf_country = gdf_country.rename(columns={"length_sum": f"{name}_length_sum"})
        # Convert from m to km
        gdf_country[f"{name}_length_sum"] = gdf_country[f"{name}_length_sum"] / 1000

        # Rasterize statistics
        da_emi = emission_vector(
            gdf=gdf_country,
            ds_like=ds_like,
            col_name=f"{name}_length_sum",
            method="value",
            mask_name="mask",
        )
        ds_country[f"{name}_length_sum_country"] = da_emi

    return ds_country


def roads_emissions_segments(
    gdf_roads: gpd.GeoDataFrame,
    ds_like: xr.Dataset,
    highway_list: List[str],
    non_highway_list: List[str] = None,
) -> xr.Dataset:
    """
    Compute roads statistics per segment/grid cell.

    Parameters
    ----------
    gdf_roads : geopandas.GeoDataFrame
        GeoDataFrame containing roads shape.

        * Required columns: ``road_type``.
    ds_like : xarray.Dataset
        Dataset at model resolution.

        * Required variables: ``mask``.
    highway_list : list of str
        List of highway road types.
    non_highway_list : list of str, optional
        List of non-highway road types. If None, all road types not in highway_list
        are considered as non-highway.

    Returns
    -------
    ds_segments : xarray.Dataset
        Dataset contains roads statistics per segment/grid cell.
        Contains the following variables:

        * ``hwy_length``: highway road length per segment/grid cell.
        * ``nnhwy_length``: non-highway road length per segment/grid cell.

    """
    # To ds_like crs
    gdf_roads = gdf_roads.to_crs(ds_like.raster.crs)
    # Common preprocessing of arguments with country
    feature_filter = _preprocess_roads_emissions(
        gdf_roads=gdf_roads,
        highway_list=highway_list,
        non_highway_list=non_highway_list,
    )

    # Loop over feature_filter
    ds_segments = xr.Dataset()
    for name, colfilter in feature_filter.items():
        logger.info(f"Computing {name} roads statistics per segment")
        # Filter gdf
        colname = [k for k in colfilter.keys()][0]
        subset_roads = gdf_roads.iloc[
            np.isin(gdf_roads[colname], colfilter.get(colname))
        ]
        ds_roads = zonal_stats_grid(
            gdf=subset_roads,
            ds_like=ds_like,
            variables=["length"],
            stats=["sum"],
            mask_name="mask",
            method="overlay",
        )
        # Convert from m to km and rename
        ds_roads[f"{name}_length"] = ds_roads["length_sum"] / 1000
        ds_roads[f"{name}_length"].attrs.update(_FillValue=0)

        # Add to ds_segments
        ds_segments[f"{name}_length"] = ds_roads[f"{name}_length"]

    return ds_segments
