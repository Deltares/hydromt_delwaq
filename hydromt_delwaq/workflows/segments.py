# -*- coding: utf-8 -*-

import numpy as np
import xarray as xr
import pandas as pd
import logging

from hydromt import flw


logger = logging.getLogger(__name__)


__all__ = ["hydromaps", "pointer"]


def hydromaps(
    hydromodel,
    logger=logger,
):
    """Returns base information maps from hydromodel.

    The following basemaps are extracted:\
    - flwdir
    - basins
    - basmsk
    
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
    ds_out["rivlen"] = hydromodel.staticmaps[hydromodel._MAPS["rivlen"]]
    ds_out["rivwth"] = hydromodel.staticmaps[hydromodel._MAPS["rivwth"]]

    basins_mv = ds_out["basins"].raster.nodata
    ds_out["basmsk"] = xr.Variable(
        dims=ds_out.raster.dims,
        data=(ds_out["basins"] != basins_mv),
        attrs=dict(_FillValue=0),
    )

    return ds_out


def pointer(
    ds_hydro,
    build_pointer=False,
    compartments=None,
    boundaries=None,
    fluxes=None,
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
    da_out : xarray.DataArray
        DataArray containing the Delwaq segment IDs.
    nrofseg : int
        Number of segments.
    """
    if compartments is None:
        ncomp = 1
        compartments = ["em"]
    else:
        ncomp = len(compartments)
    comp_ids = np.arange(1, ncomp + 1, dtype=int)

    ptid_mv = ds_hydro["basins"].raster.nodata
    np_ptid = ds_hydro["basins"].values.flatten()
    ptid = np_ptid[np_ptid != ptid_mv]
    ptid = np.arange(1, len(ptid) + 1)
    nrofseg = np.amax(ptid)
    np_ptid[np_ptid != ptid_mv] = ptid
    np_ptid = np_ptid.reshape(
        np.size(ds_hydro["basins"], 0), np.size(ds_hydro["basins"], 1)
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
        da_tot.append((da_ptid + i * nrofseg))
    da_ptid = xr.concat(
        da_tot,
        pd.Index(comp_ids, name="comp"),
    ).transpose("comp", ...)
    da_ptid.assign_coords(comp_labels=("comp", np.array(compartments)))
    # Update cells/segments/ptid based on ncomp
    nb_cell = len(ptid)
    nrofseg = nrofseg * ncomp

    # Build pointer
    if build_pointer:
        nbound = len(boundaries)
        nflux = len(fluxes)
        logger.info(
            f"Preparing pointer with {ncomp} compartments, {nbound} boundaries and {nflux} fluxes."
        )

        ### Downstream IDs ###
        # Start with searching for the ID of the downstream cells for lateral fluxes
        flwdir = flw.flwdir_from_da(ds_hydro["ldd"], ftype="infer", mask=None)
        # Boundaries
        bd_id = []
        bd_type = []
        # Keep track of the lowest boundary id value
        lowerid = 0

        da_tot = []
        np_ldd = ds_hydro["ldd"].values
        nb_out = len(np_ldd[np_ldd == 5])
        for i in range(1, ncomp + 1):
            comp_label = compartments[i - 1]
            ptiddown = flwdir.downstream(da_ptid.sel(comp=i)).astype(np.int32)
            # Outlets are boundaries and ptiddown should be negative
            outid = np.arange((-lowerid) + 1, nb_out + 1) * -1
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

        ### Add fluxes ###
        zeros = np.zeros((nb_cell, 1))
        pointer = None
        for flux in fluxes:
            flux0 = flux.split(">")[0]
            flux1 = flux.split(">")[-1]
            # Lateral flux (runoff)
            if flux0 == flux1:
                # Start building pointer with lateral fluxes (runoff)
                comp_id = comp_ids[compartments == flux0]
                ptid = da_ptid.sel(comp=comp_id).values
                ptid = ptid[ptid != ptid_mv].reshape(nb_cell, 1)
                ptiddown = da_ptiddown.sel(comp=comp_id).values
                ptiddown = ptiddown[ptiddown != ptid_mv].reshape(nb_cell, 1)
                if pointer is None:
                    pointer = np.hstack((ptid, ptiddown, zeros, zeros))
                else:
                    pointer = np.vstack(
                        (pointer, np.hstack(ptid, ptiddown, zeros, zeros))
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
                    comp_id = comp_ids[compartments == flux0]
                # Flux to boundaries
                else:
                    comp_id = comp_ids[compartments == flux1]
                ptid = da_ptid.sel(comp=comp_id).values
                ptid = ptid[ptiddown != ptid_mv].reshape(nb_cell, 1)
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
                comp_id0 = comp_ids[compartments == flux0]
                ptid0 = da_ptid.sel(comp=comp_id0).values
                ptid0 = ptid0[ptid0 != ptid_mv].reshape(nb_cell, 1)
                comp_id1 = comp_ids[compartments == flux1]
                ptid1 = da_ptid.sel(comp=comp_id1).values
                ptid1 = ptid1[ptid1 != ptid_mv].reshape(nb_cell, 1)
                if pointer is None:
                    pointer = np.hstack((ptid0, ptid1, zeros, zeros))
                else:
                    pointer = np.vstack(
                        (pointer, np.hstack(ptid0, ptid1, zeros, zeros))
                    )

        return nrofseg, da_ptid, da_ptiddown, pointer, bd_id, bd_type
    else:
        return nrofseg, da_ptid
