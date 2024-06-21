"""Worflows dealing with pointer and Delwaq segments."""

import logging
from typing import List

import numpy as np
import xarray as xr
from hydromt import flw

from .emissions import gridarea, gridlength_gridwidth

logger = logging.getLogger(__name__)


__all__ = [
    "hydromaps",
    "geometrymaps",
    "pointer",
]


def hydromaps(
    hydromodel,
    mask: str = "basins",
):
    """Return base information maps from hydromodel.

    The following basemaps are extracted:\
    - flwdir
    - basins
    - basmsk
    - elevtn
    - rivmsk

    Parameters
    ----------
    hydromodel : hydromt.model
        HydroMT Model class containing the hydromodel to build DelwaqModel from.
    mask: str, optional
        Name of the mask to use to mask the maps. Either 'rivers' or 'basins' (default).

    Returns
    -------
    ds_out : xarray.Dataset
        Dataset containing gridded emission map at model resolution.
    """
    ds_out = hydromodel.grid[hydromodel._MAPS["flwdir"]].rename("flwdir").to_dataset()
    ds_out["basins"] = hydromodel.grid[hydromodel._MAPS["basins"]]
    ds_out["river"] = hydromodel.grid[hydromodel._MAPS["rivmsk"]]
    ds_out["rivlen"] = hydromodel.grid[hydromodel._MAPS["rivlen"]]
    ds_out["rivwth"] = hydromodel.grid[hydromodel._MAPS["rivwth"]]
    ds_out["elevtn"] = hydromodel.grid[hydromodel._MAPS["elevtn"]]

    # Surface area maps
    ds_out["rivarea"] = ds_out["rivlen"] * ds_out["rivwth"]
    ds_out["rivarea"].raster.set_nodata(ds_out["rivlen"].raster.nodata)
    if "LakeArea" in hydromodel.grid:
        ds_out["lakearea"] = hydromodel.grid["LakeArea"]
    if "ResSimpleArea" in hydromodel.grid:
        ds_out["resarea"] = hydromodel.grid["ResSimpleArea"]

    basins_mv = ds_out["basins"].raster.nodata
    ds_out["basmsk"] = xr.Variable(
        dims=ds_out.raster.dims,
        data=(ds_out["basins"] != basins_mv),
        attrs=dict(_FillValue=False),
    )
    river_mv = ds_out["river"].raster.nodata
    ds_out["rivmsk"] = xr.Variable(
        dims=ds_out.raster.dims,
        data=(ds_out["river"] != river_mv),
        attrs=dict(_FillValue=False),
    )

    # Add mask
    if mask == "rivers":
        da_mask = ds_out["rivmsk"]
    else:
        da_mask = ds_out["basmsk"]
    ds_out = ds_out.drop_vars(["rivmsk", "basmsk"])
    ds_out.coords["mask"] = da_mask
    ds_out["modelmap"] = da_mask.astype(np.int32)
    ds_out["modelmap"].raster.set_nodata(0)

    return ds_out


def maps_from_hydromodel(
    hydromodel,
    maps=["rivmsk", "lndslp", "strord", "N", "SoilThickness", "thetaS"],
    logger=logger,
):
    """Return maps from hydromodel.

    Parameters
    ----------
    hydromodel : HydroMT Model
        HydroMT Model class containing the hydromodel to get geometry data from from.
    maps: list of str
        List of variables from hydromodel to extract and extend.
        By default ['rivmsk', 'lndslp', 'strord', 'N', 'SoilThickness', 'thetaS'].

    Returns
    -------
    ds_out : xarray.Dataset
        Dataset containing gridded maps at model resolution.
    """
    ds_out = xr.Dataset()
    for m in maps:
        if f"{m}_River" in hydromodel.grid:
            ds_out[m] = hydromodel.grid[f"{m}_River"].where(
                hydromodel.grid[hydromodel._MAPS["rivmsk"]],
                hydromodel.grid[m],
            )
        elif m in hydromodel._MAPS:
            if hydromodel._MAPS[m] in hydromodel.grid:
                ds_out[m] = hydromodel.grid[hydromodel._MAPS[m]]
        elif m in hydromodel.grid:
            ds_out[m] = hydromodel.grid[m]
        else:
            logger.warning(f"Map {m} not found in hydromodel, skipping.")
        # Check if m is 3D and split into several variables
        if m in ds_out and len(ds_out[m].shape) == 3:
            # Find the name of the extra dimension
            dim0 = ds_out[m].raster.dim0
            dim0_values = ds_out[dim0].values
            for i in range(ds_out[m].shape[0]):
                ds_out[f"{m}_{dim0}_{dim0_values[i]}"] = ds_out[m][i]
            ds_out = ds_out.drop_vars([m, dim0])

    return ds_out


def geometrymaps(
    hydromodel,
):
    """Return geometry information maps from hydromodel.

    The following basemaps are extracted:\
    - surface
    - length
    - width
    - latitude_grid
    - optional: bankfull_volume if bankfull depth `rivdph` is available

    Parameters
    ----------
    hydromodel : HydroMT Model
        HydroMT Model class containing the hydromodel to get geometry data from from.

    Returns
    -------
    ds_out : xarray.Dataset
        Dataset containing gridded geometry map at model resolution.
    """
    ### Geometry data ###
    surface = gridarea(hydromodel.grid)
    surface.raster.set_nodata(-9999.0)
    surface = surface.rename("surface")
    length, width = gridlength_gridwidth(hydromodel.grid)

    rivlen = hydromodel.grid[hydromodel._MAPS["rivlen"]]
    rivwth = hydromodel.grid[hydromodel._MAPS["rivwth"]]
    rivmsk = hydromodel.grid[hydromodel._MAPS["rivmsk"]]
    # surface
    surface = surface.where(rivmsk == False, rivlen * rivwth)
    surface = surface.rename("surface")
    # Add waterbodies to surface
    for wb in ["LakeArea", "ResSimpleArea"]:
        if wb in hydromodel.grid:
            wb_surface = hydromodel.grid[wb]
            surface = surface.where(wb_surface == wb_surface.raster.nodata, wb_surface)
    # length
    length = rivlen.where(rivmsk, length)
    length = length.rename("length")
    # width
    width = rivwth.where(rivmsk, width)
    width = width.rename("width")

    ds_out = surface.to_dataset()
    ds_out["length"] = length
    ds_out["width"] = width

    # Derive latitude of each cell
    y_dim = ds_out.raster.y_dim
    x_dim = ds_out.raster.x_dim
    # Expend latitude dimension values to all cells of ds_out
    ds_out["latitude_grid"] = xr.DataArray(
        data=ds_out[y_dim].values.reshape(-1, 1).repeat(ds_out.dims[x_dim], axis=1),
        coords=ds_out.raster.coords,
        dims=ds_out.raster.dims,
    )
    ds_out["latitude_grid"].raster.set_nodata(-9999.0)

    # Bankfull volume
    if hydromodel._MAPS["rivdph"] in hydromodel.grid:
        ds_out["bankfull_volume"] = (
            rivlen * rivwth * hydromodel.grid[hydromodel._MAPS["rivdph"]]
        )
        ds_out["bankfull_volume"].raster.set_nodata(ds_out["length"].raster.nodata)

    return ds_out


def pointer(
    ds_hydro: xr.Dataset,
    build_pointer: bool = False,
    surface_water: str = "sfw",
    boundaries: List[str] = ["bd"],
    fluxes: List[str] = ["sfw>sfw", "bd>sfw"],
    logger=logger,
):
    """Return map with Delwaq segment ID and pointer.

    The pointer matrix is built only if ``build_pointer`` is True.

    Parameters
    ----------
    ds_hydro : xr.Dataset
        Dataset of the hydromaps, contains 'basins', 'ldd', 'modelmap'.
    build_pointer: boolean, optional
        Boolean to build a pointer file (delwaq) or not (demission).
        If True, compartments, boundaries and fluxes lists must be provided.
    surface_water : str, optional
        Name of the surface water layer. By default 'sfw'.
    boundaries: list of str, optional
        List of names of boundaries to include. By default a unique boundary called
        'bd'.
    fluxes: list of str
        List of fluxes to include between surface water/boundaries. Name convention
        is '{surface_water_name}>{boundary_name}' for a flux from the surface water to a
        boundary, ex 'sfw>bd'. By default ['sfw>sfw', 'bd>sfw'] for runoff and inwater.
        Names in the fluxes list should match name in the hydrology_fn source in
        setup_hydrology_forcing.

    Returns
    -------
    nrofseg : int
        Number of segments.
    da_ptid : xarray.DataArray
        DataArray containing the Delwaq segment IDs.
    da_ptiddown : xarray.DataArray
        DataArray containing the Delwaq downstream segment IDs.
    pointer: numpy.ndarray, optional
        Delwaq pointer array (4 columns) for exchanges.
    bd_id: numpy.array, optional
        Array with Delwaq boundary IDs.
    bd_type: numpy.array, optional
        Array ith Delwaq boundary names.
    """
    # Prepare segments ID layer ptid based on mask of active cells
    ptid_mv = 0
    np_ptid = ds_hydro["mask"].values.flatten().astype(np.int32)
    ptid = np_ptid[np_ptid != ptid_mv]
    ptid = np.arange(1, len(ptid) + 1)
    nrofseg = np.amax(ptid)
    np_ptid[np_ptid != ptid_mv] = ptid
    np_ptid = np_ptid.reshape(
        np.size(ds_hydro["mask"], 0), np.size(ds_hydro["mask"], 1)
    )
    da_ptid = xr.DataArray(
        data=np_ptid,
        coords=ds_hydro.raster.coords,
        dims=ds_hydro.raster.dims,
        attrs=dict(_FillValue=ptid_mv),
    )

    ### Downstream IDs ###
    # Start with searching for the ID of the downstream cells for lateral fluxes
    flwdir = flw.flwdir_from_da(ds_hydro["ldd"], ftype="infer", mask=None)
    # Boundaries
    bd_id = []
    bd_type = []
    # Keep track of the lowest boundary id value
    lowerid = 0

    # Apply mask to ldd
    da_ldd = ds_hydro["ldd"].where(ds_hydro["mask"], ds_hydro["ldd"].raster.nodata)
    np_ldd = da_ldd.values
    # Number of outlets
    nb_out = len(np_ldd[np_ldd == 5])

    ptiddown = flwdir.downstream(da_ptid).astype(np.int32)
    # Remask cells draining to rivers
    ptiddown[np_ptid == ptid_mv] = ptid_mv
    # Outlets are boundaries and ptiddown should be negative
    outid = np.arange((-lowerid) + 1, (-lowerid) + nb_out + 1) * -1
    ptiddown[np_ldd == 5] = outid
    lowerid = outid[-1]
    bd_id = np.append(bd_id, (outid * (-1)))
    bd_type = np.append(bd_type, [f"{surface_water}>out{id}" for id in bd_id])
    # Add ptiddown to xarray
    da_ptiddown = xr.DataArray(
        data=ptiddown,
        coords=ds_hydro.raster.coords,
        dims=ds_hydro.raster.dims,
        attrs=dict(_FillValue=ptid_mv),
    )

    # Build pointer
    if build_pointer:
        nbound = len(boundaries)
        nflux = len(fluxes)
        logger.info(f"Preparing pointer with {nbound} boundaries and {nflux} fluxes.")
        ### Add fluxes ###
        zeros = np.zeros((nrofseg, 1))
        pointer = None
        for flux in fluxes:
            flux0 = flux.split(">")[0]
            flux1 = flux.split(">")[-1]
            # Lateral flux (runoff)
            if flux0 == flux1:
                # Start building pointer with lateral fluxes (runoff)
                ptid = da_ptid.values
                ptid = ptid[ptid != ptid_mv].reshape(nrofseg, 1)
                ptiddown = da_ptiddown.values
                ptiddown = ptiddown[ptiddown != ptid_mv].reshape(nrofseg, 1)
                if pointer is None:
                    pointer = np.hstack((ptid, ptiddown, zeros, zeros))
                else:
                    pointer = np.vstack(
                        (pointer, np.hstack((ptid, ptiddown, zeros, zeros)))
                    )
            # Flux from/to boundaries
            elif flux0 != surface_water or flux1 != surface_water:
                # The boundary cells all have the same ID
                boundid = lowerid - 1
                lowerid = boundid
                bd_id = np.append(bd_id, ([boundid * (-1)]))
                bd_type = np.append(bd_type, ([flux]))
                boundid = np.repeat(boundid, nrofseg).reshape(nrofseg, 1)
                # Flux from boundaries
                ptid = da_ptid.values
                ptid = ptid[ptid != ptid_mv].reshape(nrofseg, 1)
                if flux0 != surface_water:
                    pointerbd = np.hstack((boundid, ptid, zeros, zeros))
                else:
                    pointerbd = np.hstack((ptid, boundid, zeros, zeros))
                if pointer is None:
                    pointer = pointerbd
                else:
                    pointer = np.vstack((pointer, pointerbd))
            else:
                raise ValueError(
                    f"Flux {flux} should be of type{surface_water}>{surface_water}"
                    f"or {surface_water}>boundary or boundary>{surface_water}."
                )

        return nrofseg, da_ptid, da_ptiddown, pointer, bd_id, bd_type
    else:
        return nrofseg, da_ptid, da_ptiddown
