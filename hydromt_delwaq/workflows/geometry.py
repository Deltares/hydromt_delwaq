"""Prepare geometry data."""

import logging

import numpy as np
import xarray as xr

from .emissions import gridarea

logger = logging.getLogger(__name__)


__all__ = ["compute_geometry"]


def compute_geometry(
    ds: xr.Dataset,
    mask: xr.DataArray,
    fpaved_name: str = "soil_compacted_fraction",
    fopenwater_name: str = "land_water_fraction",
) -> np.ndarray:
    """Compute geometry data for demission.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset containing the geometry data.

        * Required variables: fraction of paved area "soil_compacted_fraction", fraction
        of open water "land_water_fraction"
    mask : xr.DataArray
        Mask to select the active cells (segments) in ds.
    fpaved_name : str
        Name of the variable in ds representing the fraction of paved area.
        By default 'soil_compacted_fraction'.
    fopenwater_name : str
        Name of the variable in ds representing the fraction of open water area.
        By default 'land_water_fraction'.

    Returns
    -------
    geometry : np.ndarray
        Array containing the geometry data.
    """
    surface = gridarea(ds)
    surface.raster.set_nodata(np.nan)
    geom = surface.rename("surface").to_dataset()
    # For EM type build a pointer like object and add to self.geometry
    geom["fPaved"] = ds[fpaved_name]
    geom["fOpenWater"] = ds[fopenwater_name]
    geom["fUnpaved"] = (
        (geom["fPaved"] * 0.0 + 1.0) - geom["fPaved"] - geom["fOpenWater"]
    )
    geom["fUnpaved"] = xr.where(geom["fUnpaved"] < 0.0, 0.0, geom["fUnpaved"])

    mask = mask.values.flatten()
    for dvar in ["surface", "fPaved", "fUnpaved", "fOpenWater"]:
        data = geom[dvar].values.flatten()
        data = data[mask].reshape(len(data[mask]), 1)
        if dvar == "surface":
            geometry = data
        else:
            geometry = np.hstack((geometry, data))

    return geometry
