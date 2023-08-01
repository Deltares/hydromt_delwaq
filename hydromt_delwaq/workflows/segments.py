# -*- coding: utf-8 -*-

import numpy as np
import xarray as xr
import pandas as pd
import logging
from typing import List

from hydromt import flw
from hydromt.raster import full_like
from . import emissions


logger = logging.getLogger(__name__)


__all__ = [
    "hydromaps",
    "geometrymaps",
    "pointer",
    "extend_comp_with_zeros",
    "extend_comp_with_duplicates",
    "sfwcomp",
]


def hydromaps(
    hydromodel,
    logger=logger,
):
    """Returns base information maps from hydromodel.

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

    Returns
    -------
    ds_out : xarray.Dataset
        Dataset containing gridded emission map at model resolution.
    """
    ds_out = (
        hydromodel.staticmaps[hydromodel._MAPS["flwdir"]].rename("flwdir").to_dataset()
    )
    ds_out["basins"] = hydromodel.staticmaps[hydromodel._MAPS["basins"]]
    ds_out["river"] = hydromodel.staticmaps[hydromodel._MAPS["rivmsk"]]
    ds_out["rivlen"] = hydromodel.staticmaps[hydromodel._MAPS["rivlen"]]
    ds_out["rivwth"] = hydromodel.staticmaps[hydromodel._MAPS["rivwth"]]
    ds_out["elevtn"] = hydromodel.staticmaps[hydromodel._MAPS["elevtn"]]

    # Surface area maps
    ds_out["rivarea"] = ds_out["rivlen"] * ds_out["rivwth"]
    if "LakeArea" in hydromodel.staticmaps:
        ds_out["lakearea"] = hydromodel.staticmaps["LakeArea"]
    if "ResSimpleArea" in hydromodel.staticmaps:
        ds_out["resarea"] = hydromodel.staticmaps["ResSimpleArea"]

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

    return ds_out


def maps_from_hydromodel(
    hydromodel,
    compartments,
    maps=["rivmsk", "lndslp", "strord", "N", "SoilThickness", "thetaS"],
    logger=logger,
):
    """Returns maps from hydromodel and extend to all compartments.

    Parameters
    ----------
    hydromodel : HydroMT Model
        HydroMT Model class containing the hydromodel to get geometry data from from.
    compartments : list of str, optional
        List of names of compartments to include. By default one for surface waters called 'sfw'.
    maps: list of str
        List of variables from hydromodel to extract and extend.
        By default ['rivmsk', 'lndslp', 'strord', 'N', 'SoilThickness', 'thetaS'].

    Returns
    -------
    ds_out : xarray.Dataset
        Dataset containing gridded geometry map at model resolution for all compartments.
    """
    ds_out = xr.Dataset()
    for m in maps:
        if f"{m}_River" in hydromodel.staticmaps:
            ds_out[m] = hydromodel.staticmaps[f"{m}_River"].where(
                hydromodel.staticmaps[hydromodel._MAPS["rivmsk"]],
                hydromodel.staticmaps[m],
            )
        elif m in hydromodel._MAPS:
            if hydromodel._MAPS[m] in hydromodel.staticmaps:
                ds_out[m] = hydromodel.staticmaps[hydromodel._MAPS[m]]
        else:
            ds_out[m] = hydromodel.staticmaps[m]

    # surface water comp
    ds_out = extend_comp_with_duplicates(ds1c=ds_out, compartments=compartments)

    return ds_out


def geometrymaps(
    hydromodel,
    compartments,
    comp_attributes,
    logger=logger,
):
    """Returns geometry information maps for all compartments from hydromodel.

    The following basemaps are extracted:\
    - surface
    - length
    - width
    - slope
    - river
    - streamorder
    - manning
    
    Parameters
    ----------
    hydromodel : HydroMT Model
        HydroMT Model class containing the hydromodel to get geometry data from from.
    compartments : list of str, optional
        List of names of compartments to include. By default one for surface waters called 'sfw'.
    comp_attributes: list of int
        Attribute 1 value of the B3_attributes config file. 1 or 0 for surface water. Also used to compute surface variable.
 
    Returns
    -------
    ds_out : xarray.Dataset
        Dataset containing gridded geometry map at model resolution for all compartments.
    """
    ### Geometry data ###
    surface = emissions.gridarea(hydromodel.staticmaps)
    surface.raster.set_nodata(-9999.0)
    surface = surface.rename("surface")
    length, width = emissions.gridlength_gridwidth(hydromodel.staticmaps)
    surface_tot = []
    length_tot = []
    width_tot = []
    for i in range(len(compartments)):
        if comp_attributes[i] == 0:  # surface water
            rivlen = hydromodel.staticmaps[hydromodel._MAPS["rivlen"]]
            rivwth = hydromodel.staticmaps[hydromodel._MAPS["rivwth"]]
            rivmsk = hydromodel.staticmaps[hydromodel._MAPS["rivmsk"]]
            # surface
            rivsurface = surface.where(rivmsk == False, rivlen * rivwth)
            rivsurface = rivsurface.rename("surface")
            surface_tot.append(rivsurface)
            # length
            rivlength = rivlen.where(rivmsk, rivlen)
            rivlength = rivlength.rename("length")
            length_tot.append(rivlength)
            # width
            rivwidth = rivwth.where(rivmsk, width)
            rivwidth = rivwidth.rename("width")
            width_tot.append(rivwidth)
        else:  # other use direct cell surface
            surface_tot.append(surface)
            length_tot.append(length)
            width_tot.append(width)
    ds_out = (
        xr.concat(
            surface_tot,
            pd.Index(np.arange(1, len(compartments) + 1, dtype=int), name="comp"),
        ).transpose("comp", ...)
    ).to_dataset()
    ds_out["length"] = xr.concat(
        length_tot,
        pd.Index(np.arange(1, len(compartments) + 1, dtype=int), name="comp"),
    ).transpose("comp", ...)
    ds_out["width"] = xr.concat(
        width_tot,
        pd.Index(np.arange(1, len(compartments) + 1, dtype=int), name="comp"),
    ).transpose("comp", ...)

    return ds_out


def pointer(
    ds_hydro: xr.Dataset,
    build_pointer: bool = False,
    compartments: List[str] = None,
    boundaries: List[str] = None,
    fluxes: List[str] = None,
    logger=logger,
):
    """Returns map with Delwaq segment ID. As well as the pointer matrix if ``build_pointer`` is True.

    Parameters
    ----------
    ds_hydro : xr.Dataset
        Dataset of the hydromaps, contains 'basins', 'ldd', 'modelmap'.
    build_pointer: boolean, optional
        Boolean to build a pointer file (delwaq) or not (demission).
        If True, compartments, boundaries and fluxes lists must be provided.
    compartments : list of str, optional
        List of names of compartments to include. By default one for surface waters called 'sfw'.
    boundaries: list of str, optional
        List of names of boundaries to include. By default a unique boundary called 'bd'.
    fluxes: list of str
        List of fluxes to include between compartments/boundaries. Name convention is '{compartment_name}>{boundary_name}'
        for a flux from a compartment to a boundary, ex 'sfw>bd'. By default ['sfw>sfw', 'bd>sfw'] for runoff and inwater.
        Names in the fluxes list should match name in the hydrology_fn source in setup_hydrology_forcing.

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
    if compartments is None:
        ncomp = 1
        compartments = ["em"]
    else:
        ncomp = len(compartments)
    comp_ids = np.arange(1, ncomp + 1, dtype=int)

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

    # Add other compartments
    da_tot = [da_ptid]
    for i in np.arange(1, ncomp):
        da_comp = xr.where(da_ptid == ptid_mv, 0, da_ptid + (i * nrofseg))
        da_tot.append(da_comp)
    da_ptid = xr.concat(
        da_tot,
        pd.Index(comp_ids, name="comp"),
    ).transpose("comp", ...)
    da_ptid.assign_coords(comp_labels=("comp", np.array(compartments)))
    # Update cells/segments/ptid based on ncomp
    nb_cell = nrofseg  # len(ptid)
    nrofseg = nrofseg * ncomp

    ### Downstream IDs ###
    # Start with searching for the ID of the downstream cells for lateral fluxes
    flwdir = flw.flwdir_from_da(ds_hydro["ldd"], ftype="infer", mask=None)
    # Boundaries
    bd_id = []
    bd_type = []
    # Keep track of the lowest boundary id value
    lowerid = 0

    da_tot = []
    # Apply mask to ldd
    da_ldd = ds_hydro["ldd"].where(ds_hydro["mask"], ds_hydro["ldd"].raster.nodata)
    np_ldd = da_ldd.values
    # Number of outlets
    nb_out = len(np_ldd[np_ldd == 5])
    for i in range(1, ncomp + 1):
        comp_label = compartments[i - 1]
        ptiddown = flwdir.downstream(da_ptid.sel(comp=i)).astype(np.int32)
        # Remask cells draining to rivers
        ptiddown[np_ptid == ptid_mv] = ptid_mv
        # Outlets are boundaries and ptiddown should be negative
        outid = np.arange((-lowerid) + 1, (-lowerid) + nb_out + 1) * -1
        ptiddown[np_ldd == 5] = outid
        lowerid = outid[-1]
        bd_id = np.append(bd_id, (outid * (-1)))
        bd_type = np.append(bd_type, [f"{comp_label}>out{id}" for id in bd_id])
        # Add ptiddown to xarray
        da_ptiddown = xr.DataArray(
            data=ptiddown,
            coords=ds_hydro.raster.coords,
            dims=ds_hydro.raster.dims,
            attrs=dict(_FillValue=ptid_mv),
        )
        da_tot.append(da_ptiddown)
    da_ptiddown = xr.concat(
        da_tot,
        pd.Index(np.arange(1, ncomp + 1, dtype=int), name="comp"),
    ).transpose("comp", ...)
    da_ptiddown.assign_coords(comp_labels=("comp", np.array(compartments)))

    # Build pointer
    if build_pointer:
        nbound = len(boundaries)
        nflux = len(fluxes)
        logger.info(
            f"Preparing pointer with {ncomp} compartments, {nbound} boundaries and {nflux} fluxes."
        )
        ### Add fluxes ###
        zeros = np.zeros((nb_cell, 1))
        pointer = None
        for flux in fluxes:
            flux0 = flux.split(">")[0]
            flux1 = flux.split(">")[-1]
            # Lateral flux (runoff)
            if flux0 == flux1:
                # Start building pointer with lateral fluxes (runoff)
                comp_id = comp_ids[compartments.index(flux0)]
                ptid = da_ptid.sel(comp=comp_id).values
                ptid = ptid[ptid != ptid_mv].reshape(nb_cell, 1)
                ptiddown = da_ptiddown.sel(comp=comp_id).values
                ptiddown = ptiddown[ptiddown != ptid_mv].reshape(nb_cell, 1)
                if pointer is None:
                    pointer = np.hstack((ptid, ptiddown, zeros, zeros))
                else:
                    pointer = np.vstack(
                        (pointer, np.hstack((ptid, ptiddown, zeros, zeros)))
                    )
            # Flux from/to boundaries
            elif flux0 not in compartments or flux1 not in compartments:
                # The boundary cells all have the same ID
                boundid = lowerid - 1
                lowerid = boundid
                bd_id = np.append(bd_id, ([boundid * (-1)]))
                bd_type = np.append(bd_type, ([flux]))
                boundid = np.repeat(boundid, nb_cell).reshape(nb_cell, 1)
                # Flux from boundaries
                if flux0 not in compartments:
                    comp_id = comp_ids[compartments.index(flux1)]
                # Flux to boundaries
                else:
                    comp_id = comp_ids[compartments.index(flux0)]
                ptid = da_ptid.sel(comp=comp_id).values
                ptid = ptid[ptid != ptid_mv].reshape(nb_cell, 1)
                if flux0 not in compartments:
                    pointerbd = np.hstack((boundid, ptid, zeros, zeros))
                else:
                    pointerbd = np.hstack((ptid, boundid, zeros, zeros))
                if pointer is None:
                    pointer = pointerbd
                else:
                    pointer = np.vstack((pointer, pointerbd))
            # Flux from comp to comp
            else:
                comp_id0 = comp_ids[compartments.index(flux0)]
                ptid0 = da_ptid.sel(comp=comp_id0).values
                ptid0 = ptid0[ptid0 != ptid_mv].reshape(nb_cell, 1)
                comp_id1 = comp_ids[compartments.index(flux1)]
                ptid1 = da_ptid.sel(comp=comp_id1).values
                ptid1 = ptid1[ptid1 != ptid_mv].reshape(nb_cell, 1)
                if pointer is None:
                    pointer = np.hstack((ptid0, ptid1, zeros, zeros))
                else:
                    pointer = np.vstack(
                        (pointer, np.hstack((ptid0, ptid1, zeros, zeros)))
                    )

        return nrofseg, da_ptid, da_ptiddown, pointer, bd_id, bd_type
    else:
        return nrofseg, da_ptid, da_ptiddown


def extend_comp_with_zeros(ds1c: xr.Dataset, comp_ds1c: str, compartments: list):
    """
    Transform 2D data (y_dim, x_dim) into 3D data (comp, y_dim, x_dim).
    Data from ds1c is assign to comp named comp_ds1c and other compartments are filled with zeros.

    Parameters
    ----------
    ds1c: xr.Dataset
        2D gridded data
    comp_ds1c: str
        Compartment name in compartments to which data in ds1c is assigned to.
    compartemnts: list
        List with compartments names.

    Returns
    -------
    3D gridded data with zeros for compartments other than comp_ds1c

    """
    for var in ds1c.data_vars:
        da_zeros = full_like(ds1c[var], lazy=True)
        da_zeros = da_zeros.where(ds1c[var] == ds1c[var].raster.nodata, 0)
        da_tot = []
        for i in range(len(compartments)):
            if compartments[i] == comp_ds1c:  # compartment with input
                da_tot.append(ds1c[var])
            else:  # other compartments to fill with zeros
                da_tot.append(da_zeros)
        da = xr.concat(
            da_tot,
            pd.Index(np.arange(1, len(compartments) + 1, dtype=int), name="comp"),
        ).transpose("comp", ...)
        ds1c[var] = da

    return ds1c


def extend_comp_with_duplicates(ds1c, compartments):
    """
    Transform 2D data (y_dim, x_dim) into 3D data (comp, y_dim, x_dim).
    Data from ds1c is duplicated for all compartments.

    Parameters
    ----------
    ds1c: xr.Dataset
        2D gridded data
    compartemnts: list
        List with compartments names.

    Returns
    -------
    3D gridded data with ds1c duplicated for all compartments

    """
    for var in ds1c.data_vars:
        da_tot = []
        for i in range(len(compartments)):
            da_tot.append(ds1c[var])
        da = xr.concat(
            da_tot,
            pd.Index(np.arange(1, len(compartments) + 1, dtype=int), name="comp"),
        ).transpose("comp", ...)
        ds1c[var] = da

    return ds1c


def sfwcomp(compartments: list, config: dict):
    """Finds and return surface water compartment based on B3_attributes config"""
    nl = 7
    sfw = None
    attributes = config.get("B3_attributes")
    for i in range(len(compartments)):
        cp = attributes.get(
            f"l{nl}",
            "     1*01 ; sfw",
        )
        issfw = cp.split("*")[1][0:2]
        if int(issfw) == 1:
            sfw = cp.split(";")[1][1:]
        nl += 1

    return sfw
