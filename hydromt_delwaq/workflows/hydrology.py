import pandas as pd
import xarray as xr
import logging
from typing import List
from xclim.indices.stats import frequency_analysis

import hydromt

from . import emissions

logger = logging.getLogger(__name__)

__all__ = [
    "hydrology_forcing",
    "bankfull_volume",
    "bankfull_volume_from_discharge",
]


def hydrology_forcing(
    ds: xr.Dataset,
    ds_model: xr.Dataset,
    time_tuple: tuple,
    timestepsecs: int,
    fluxes: List,
    surface_water: str = "sfw",
    add_volume_offset: bool = True,
    min_volume: float = 0.1,
    override: List = [],
    logger: logging.Logger = logger,
):
    """Calculate hydrology forcing data.

    As the fluxes order should precisely macth the pointer defined in setup_basemaps,
    the variables names in ``hydro_forcing_fn`` should match names defined in the
    ``fluxes`` argument of setup_basemaps. These names should also have been saved in
    the file config/B7_fluxes.inc.

    If several sub-variables in ``hydro_forcing_fn`` need to be summed up to get the
    expected flux in pointer, they can be named {flux_name_in_pointer}_{number}
    (eg "sfw>sfw_1" and "sfw>sfw_2" to get "sfw>sfw") and the function will sum them on
    the fly. To remove (- instead of +) use unit_mult attribute of the data catalog with
    -1 for the sub-variables of interest.
    To override rather than sum fluxes, use the ``override`` argument and name of the
    concerned flux (eg "sfw>sfw"). The flux are overwritten (excluding nodata) in the
    order of the list, so the last one will be the one used.

    Unit conversions are possible from mm/area to m3/s for fluxes, volumes should be
    provided directly in m3. For conversion from mm to m3/s, it is possible to specify
    over wich surface area the mm are calculated. If 'mm' (default), the cellarea is
    assumed. Else, you can use 'mm/{surfacearea}' where {surfacearea} should be a map
    available in hydromaps (rivarea, lakearea, resarea).

    Parameters
    ----------
    ds : xr.Dataset
        Dataset containing the input data
    ds_model : xr.Dataset
        Dataset containing the model grid

        * Required variables: [mask]

        * Optional variables for unit conversion: [rivarea, lakearea, resarea]

    time_tuple : tuple
        Tuple containing the start and end time of the simulation
    timestepsecs : int
        Timestep of the simulation in seconds
    fluxes : list
        List of fluxes to be calculated
    surface_water : str, optional
        Name of the surface water layer to calculate volume. By default "sfw".
    add_volume_offset : bool, optional
        Add time offset of one timestepsecs to the volumes, by default True.
    min_volume : float, optional
        Minimum volume to be considered, by default 0.1 to avoid zero volumes.
    override : list, optional
        List of fluxes to be overriden if several fluxes varibales are found, by default
        [].
    logger : logging.Logger, optional
        Logger object, by default logger

    Returns
    -------
    ds_out : xr.Dataset
        Dataset containing the input data
    """
    # If needed reproject forcing data to model grid
    # As forcing should come from hydromodel the grids should already be identical
    if not ds.raster.identical_grid(ds_model):
        logger.warning(
            "hydro_forcing_fn and model grid are not identical. Reprojecting."
        )
        ds = ds.raster.reproject_like(ds_model)

    # Add _FillValue to the data attributes
    for dvar in ds.data_vars.keys():
        nodata = ds[dvar].raster.nodata
        if nodata is not None:
            ds[dvar].attrs.update(_FillValue=nodata)
        else:
            ds[dvar].attrs.update(_FillValue=-9999.0)

    # Copy of ds to be filled
    dsvar = [v for v in ds.data_vars if v.startswith(fluxes[0])][0]
    ds_out = hydromt.raster.full_like(ds[dsvar], lazy=True).to_dataset()
    ds_out = ds_out.sel(time=slice(*time_tuple))
    # Array of zeros
    da_zeros = ds[dsvar] * 0.0
    da_zeros.attrs.update(unit="m3/s")
    da_zeros.raster.set_nodata(-999.0)

    ### Fluxes ###
    for flux in fluxes:
        # Check if the flux is split into several variables
        fl_vars = [v for v in ds.data_vars if v.split("_")[0] == flux]
        # Sort the list of flux by name, important for override type of sum
        fl_vars.sort()
        # Unit conversion (from mm to m3/s)
        for fl in fl_vars:
            # Assume by default mm/cellarea
            unit = ds[fl].attrs.get("unit")
            if unit == "mm":
                surface = emissions.gridarea(ds)
                ds[fl] = ds[fl] * surface / (1000 * timestepsecs)
            # For other area eg river, lake, reservoir
            elif unit.startswith("mm/"):
                surfacemap = unit.split("/")[1]
                if surfacemap in ds_model:
                    ds[fl] = ds[fl] * ds_model[surfacemap] / (1000 * timestepsecs)
                else:
                    surface_fns = [f for f in ds_model.keys() if f.endswith("area")]
                    logger.error(
                        f"Map {surfacemap} not found in hydromaps to convert unit {unit}. Allowed names are {surface_fns}."
                    )
        if len(fl_vars) > 0:  # flux not in ds
            attrs = ds[fl_vars[0]].attrs.copy()
        else:
            logger.warning(f"Flux {flux} not found in hydro_forcing_fn. Using zeros.")
            ds[flux] = da_zeros
            attrs = da_zeros.attrs.copy()
        if len(fl_vars) > 1:  # need to sum or override
            if flux not in override:
                ds[flux] = ds[fl_vars].fillna(0).to_array().sum("variable")
            else:
                ds[flux] = ds[fl_vars[0]]
                for fl in fl_vars[1:]:
                    ds[flux] = ds[flux].where(ds[fl].isnull(), ds[fl])
        ds_out[flux] = ds[flux].sel(time=slice(*time_tuple))
        ds_out[flux].attrs.update(attrs)
        ds_out[flux].attrs.update(unit="m3/s")

    ### Volumes ###
    # Add offset for the volumes if needed
    if add_volume_offset:
        # Get the freq of ds and add + 1 offset
        times = pd.to_datetime(ds["time"].values)
        times.freq = pd.infer_freq(times)
        times_offset = times + times.freq
    vol = surface_water
    # Check if the flux is split into several variables
    vol_vars = [v for v in ds.data_vars if v.split("_")[0] == vol]
    # Unit conversion (from mm to m3)
    for vl in vol_vars:
        unit = ds[vl].attrs.get("unit")
        if unit == "mm":
            surface = emissions.gridarea(ds)
            ds[vl] = ds[vl] * surface / (1000 * timestepsecs)
        # For other area eg river, lake, reservoir
        elif unit.startswith("mm/"):
            surfacemap = unit.split("/")[1]
            if surfacemap in ds_model:
                ds[vl] = ds[vl] * ds_model[surfacemap] / (1000 * timestepsecs)
            else:
                surface_fns = [f for f in ds_model.keys() if f.endswith("area")]
                logger.error(
                    f"Map {surfacemap} not found in hydromaps to convert unit {unit}."
                    f"Allowed names are {surface_fns}."
                )
        attrs = ds[vol_vars[0]].attrs.copy()
        if len(vol_vars) > 1:  # need to sum
            if vol not in override:
                ds[vol] = ds[vol_vars].fillna(0).to_array().sum("variable")
            else:
                ds[vol] = ds[vol_vars[0]]
                for vl in vol_vars[1:]:
                    ds[vol] = ds[vol].where(ds[vl].isnull(), ds[vl])
        # In order to avoid zero volumes
        # a basic minimum value of 0.0001 m3 is added to all volumes
        ds[vol] = ds[vol] + min_volume
        # Add offset for the volumes if needed
        if add_volume_offset:
            da_vol = ds[vol].copy()
            ds = ds.drop_vars(vol)
            da_vol["time"] = times_offset
        else:
            da_vol = ds[vol]
        ds_out[vol] = da_vol.sel(time=slice(*time_tuple))
        ds_out[vol].attrs.update(attrs)
        ds_out[vol].attrs.update(unit="m3")

    # Select variables, needed??
    variables = fluxes.copy()
    variables.extend([surface_water])
    ds_out = ds_out[variables]

    ds_out.coords["mask"] = xr.DataArray(
        dims=ds_out.raster.dims,
        coords=ds_out.raster.coords,
        data=ds_model["mask"].values,
        attrs=dict(_FillValue=ds_model["mask"].raster.nodata),
    )

    return ds_out


def bankfull_volume_from_discharge(
    da_q: xr.DataArray,
    ds_model: xr.Dataset = None,
    bankfull_rp: int = 2,
    distribution: str = "genextreme",
):
    """
    Derive bankfull discharge/depth/volume from a time series of discharge.

    The bankfull discharge is defined as the discharge corresponding to a given return period.
    The power law method of hydromt.workflows.river_depth is used to get bankfull depth.

    Parameters
    ----------
    da_q : xr.DataArray
        Time series of discharge [m3/s].
    ds_model: xr.Dataset
        Model dataset with width [m] and length [m] to derive bankfull volume.

        * Required variables: ['width', 'length']

    bankfull_rp : int, optional
        Return period of the bankfull discharge, by default 2 years
    distribution : str, optional
        Extreme value distribution, by default "genextreme" for Generalized Extreme Value.

    Returns
    -------
    ds_out: xr.Dataset
        Dataset containing bankfull_discharge [m3/s], bankfull_depth [m] and bankfull_volume [m3].
    """
    # Make sure ds_model is 2D and not 3D in case one comp
    ds_model = ds_model.squeeze()
    # Obtain bankfull dicharge using xclim
    da_rps = frequency_analysis(
        da_q, t=bankfull_rp, dist=distribution, mode="max", freq="YS"
    ).squeeze()

    # Renaming and simplifying
    da_rps.name = "qbankfull"
    ds_out = da_rps.to_dataset()

    # Derive bankfull depth
    # power law method - might be better to use manning / wflow kin wave method
    ds_out["bankfull_depth"] = hydromt.workflows.river_depth(
        data=ds_out, method="powlaw", min_rivdph=1.0
    )

    # Derive bankfull volume
    ds_out["bankfull_volume"] = (
        ds_out["bankfull_depth"] * ds_model["width"] * ds_model["length"]
    )

    return ds_out.rename({"qbankfull": "bankfull_discharge"})


def bankfull_volume(
    da_v: xr.DataArray,
    bankfull_rp: int = 2,
    distribution: str = "genextreme",
    method: str = "ML",
):
    """
    Derive bankfull volume from a time series of volume.

    The bankfull volume is defined as the volume corresponding to a given return period.
    xclim.indices.stats.frequency_analysis is used to get bankfull volume based on annual maxima method.

    Parameters
    ----------
    da_v : xr.DataArray
        Time series of volume [m3].
    bankfull_rp : int, optional
        Return period of the bankfull volume, by default 2 years
    distribution : str, optional
        Extreme value distribution, by default "genextreme" for Generalized Extreme Value.
        Available values: `beta`, `expon`, `genextreme`, `gamma`, `gumbel_r`, `lognorm`, `norm`.
    method : str, optional
        Fitting method, either maximum likelihood (ML or MLE), method of moments (MOM),
        probability weighted moments (PWM), also called L-Moments, or approximate method (APP).
        The PWM method is usually more robust to outliers. ML by default.

    Returns
    -------
    da_out: xr.Dataarray
        Datarray containing bankfull_volume [m3].
    """
    # Obtain bankfull volume using xclim frequency analysis method
    da_out = frequency_analysis(
        da_v, t=bankfull_rp, dist=distribution, mode="max", freq="YS", method=method
    ).squeeze()
    da_out.name = "bankfull_volume"
    # Remove return period coordinate and add to attributes
    da_out.attrs.update({"return_period": da_out.return_period.values})
    da_out = da_out.drop_vars("return_period")

    return da_out
