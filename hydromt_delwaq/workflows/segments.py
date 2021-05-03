# -*- coding: utf-8 -*-

import numpy as np
import xarray as xr
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
    build_pointer=True,
    logger=logger,
):
    """Returns map with Delwaq segment ID.

    Parameters
    ----------
    ds_hydro : xr.Dataset
        Dataset of the hydromaps, contains 'basins', 'ldd', 'modelmap'.

    Returns
    -------
    da_out : xarray.DataArray
        DataArray containing the Delwaq segment IDs.
    nrofseg : int
        Number of segments.
    """
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
    nb_cell = len(ptid)

    if build_pointer:
        logger.info(f"Preparing pointer with surface runoff and inwater.")
        # Start with searching for the ID of the downstream cells
        flwdir = flw.flwdir_from_da(ds_hydro["ldd"], ftype="infer", mask=None)
        ptiddown = flwdir.downstream(da_ptid).astype(np.int32)
        # Add boundaries
        bd_id = []
        bd_type = []
        # Outlets are boundaries and ptiddown should be negative
        np_ldd = ds_hydro["ldd"].values
        nb_out = len(np_ldd[np_ldd == 5])
        outid = np.arange(1, nb_out + 1) * -1
        ptiddown[np_ldd == 5] = outid
        # Keep track of the lowest boundary id value
        lowerid = outid[-1]
        bd_id = np.append(bd_id, (outid * (-1)))
        bd_type = np.append(bd_type, np.repeat("Sfw2Outflow", len(outid)))
        # Add ptiddown to xarray
        da_ptiddown = xr.DataArray(
            data=ptiddown,
            coords=ds_hydro.raster.coords,
            dims=ds_hydro.raster.dims,
            attrs=dict(_FillValue=ptid_mv),
        )

        # Start building pointer with lateral fluxes (runoff)
        ptid = ptid.reshape(nb_cell, 1)
        ptiddown = ptiddown[ptiddown != ptid_mv].reshape(nb_cell, 1)
        zeros = np.zeros((nb_cell, 1))
        pointer = np.hstack((ptid, ptiddown, zeros, zeros))

        # Add the inwater flux for mass balance conservation
        # The inwater boundaries all have the same ID
        boundid = lowerid - 1
        lowerid = boundid
        bd_id = np.append(bd_id, ([boundid * (-1)]))
        bd_type = np.append(bd_type, (["Inw2Sfw"]))
        boundid = np.repeat(boundid, nb_cell).reshape(nb_cell, 1)

        pointer = np.vstack((pointer, np.hstack((boundid, ptid, zeros, zeros))))

        return nrofseg, da_ptid, da_ptiddown, pointer, bd_id, bd_type

    else:
        return nrofseg, da_ptid
