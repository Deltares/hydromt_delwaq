"""Forcing data workflows."""

import logging
from typing import List

import hydromt
import pandas as pd
import xarray as xr

from .emissions import gridarea

logger = logging.getLogger(__name__)

__all__ = [
    "hydrology_forcing",
    "hydrology_forcing_em",
    "sediment_forcing",
    "climate_forcing",
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
    # Delwaq needs one more timestep for the volume
    starttime, endtime = pd.to_datetime(time_tuple)
    endtime = endtime + pd.to_timedelta(timestepsecs, unit="s")
    ds_out_times = pd.date_range(starttime, endtime, freq=f"{timestepsecs}s")
    coords = {
        "time": ds_out_times,
        ds.raster.y_dim: ds[ds.raster.y_dim],
        ds.raster.x_dim: ds[ds.raster.x_dim],
    }
    ds_out = hydromt.raster.full(
        coords=coords,
        nodata=-999.0,
        dtype="float32",
        name=dsvar,
        crs=ds.raster.crs,
        lazy=True,
    ).to_dataset()
    # Array of zeros
    da_zeros = xr.zeros_like(ds_out[dsvar])
    da_zeros.raster.set_nodata(-999.0)
    da_zeros.attrs.update(unit="m3/s")

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
                surface = gridarea(ds)
                ds[fl] = ds[fl] * surface / (1000 * timestepsecs)
            # For other area eg river, lake, reservoir
            elif unit.startswith("mm/"):
                surfacemap = unit.split("/")[1]
                if surfacemap in ds_model:
                    ds[fl] = ds[fl] * ds_model[surfacemap] / (1000 * timestepsecs)
                else:
                    surface_fns = [f for f in ds_model.keys() if f.endswith("area")]
                    logger.error(
                        f"Map {surfacemap} not found in hydromaps to convert unit "
                        f"{unit}. Allowed names are {surface_fns}."
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
            surface = gridarea(ds)
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
        ds_out[vol] = da_vol.sel(time=slice(starttime, endtime))
        ds_out[vol].attrs.update(attrs)
        ds_out[vol].attrs.update(unit="m3")

    # Select variables (drop dsvar if needed)
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


def hydrology_forcing_em(
    ds: xr.Dataset,
    ds_model: xr.Dataset,
    timestepsecs: int,
    include_transport: bool = True,
    logger: logging.Logger = logger,
) -> xr.Dataset:
    """
    Calculate hydrology forcing data for emission model.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset containing the EM forcing to resample.

        * Required variables without transport: ['time', 'precip', 'runPav',
            'runUnp', 'infilt']
        * Required variables with transport: ['time', 'precip', 'runPav',
            'runUnp', 'infilt', 'exfilt*', 'vwcproot', 'q_land', 'q_ss']
    ds_model : xr.Dataset
        Dataset containing the model grid.

        * Required variables: ['mask', 'modelmap']
    timestepsecs : int
        Model timestep in seconds.
    include_transport : bool, optional
        If False (default), only use the vertical fluxes for emission [precip,
        runPav, runUnp, infilt, totflw].
        If True, includes additional fluxes for land and subsurface transport
        [precip, runPav, runUnp, infilt, exfilt, vwcproot, q_land, q_ss].

    Returns
    -------
    ds_out : xr.Dataset
        Dataset containing the processed hydrology data.
    """
    # Select variables based on model type
    vars = ["precip", "runPav", "runUnp", "infilt"]
    if include_transport:
        vars.extend(["vwcproot", "q_land", "q_ss"])
        ex_vars = [v for v in ds.data_vars if v.startswith("exfilt")]
        vars.extend(ex_vars)
    ds = ds[vars]

    # Unit conversion (from mm to m3/s)
    for dvar in ds.data_vars.keys():
        if ds[dvar].attrs.get("unit") == "mm":
            attrs = ds[dvar].attrs.copy()
            surface = gridarea(ds)
            ds[dvar] = ds[dvar] * surface / (1000 * timestepsecs)
            ds[dvar].attrs.update(attrs)  # set original attributes
            ds[dvar].attrs.update(unit="m3/s")

    # Sum up exfiltrattion or add totflw depending on include_transport
    if include_transport:
        # Add exfiltration (can be split into several variables)
        # Check if the flux is split into several variables
        ex_vars = [v for v in ds.data_vars if v.startswith("exfilt")]
        if len(ex_vars) > 1:  # need to sum
            attrs = ds[ex_vars[0]].attrs.copy()
            nodata = ds[ex_vars[0]].raster.nodata
            # Cover exfilt with zeros (some negative zeros in wflow outputs?)
            ds["exfilt"] = ds[ex_vars].fillna(0).to_array().sum("variable")
            ds["exfilt"] = ds["exfilt"].where(ds["exfilt"] > 0.0, 0.0)
            ds["exfilt"].attrs.update(attrs)
            ds["exfilt"].raster.set_nodata(nodata)
            ds = ds.drop_vars(ex_vars)
        elif ex_vars[0] != "exfilt":
            ds = ds.rename({ex_vars[0]: "exfilt"})
    else:
        # Add total flow
        ds["totflw"] = ds["precip"].copy()

    # align forcing file with hydromaps
    # as hydro forcing comes from hydro model, it should be aligned with hydromaps
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

    # Add mask
    ds.coords["mask"] = xr.DataArray(
        dims=ds.raster.dims,
        coords=ds.raster.coords,
        data=ds_model["modelmap"].values,
        attrs=dict(_FillValue=ds_model["mask"].raster.nodata),
    )

    return ds


def sediment_forcing(
    ds: xr.Dataset,
    ds_model: xr.Dataset,
    timestepsecs: int,
    logger: logging.Logger = logger,
) -> xr.Dataset:
    """Calculate sediment forcing data.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset containing several particule sizes (IM1, IM2 etc)

        * Required variables: ['ErodIM*']
    ds_model : xr.Dataset
        Dataset containing the model grid

        * Required variables: [mask]
    timestepsecs : int
        Timestep of the simulation in seconds
    logger : logging.Logger, optional
        Logger object, by default logger

    Returns
    -------
    ds_out : xr.Dataset
        Dataset containing the processed sediment data.
    """
    # Select time
    freq = pd.to_timedelta(timestepsecs, unit="s")

    # Resample in time if needed
    ds_out = xr.Dataset()
    for var in ds.data_vars:
        ds_out[var] = hydromt.workflows.forcing.resample_time(
            ds[var],
            freq,
            upsampling="bfill",  # we assume right labeled original data
            downsampling="sum",
            conserve_mass=True,
            logger=logger,
        )

    # If needed reproject forcing data to model grid
    # As forcing should come from hydromodel the grids should already be identical
    if not ds_out.raster.identical_grid(ds_model):
        logger.warning(
            "hydro_forcing_fn and model grid are not identical. Reprojecting."
        )
        ds_out = ds_out.raster.reproject_like(ds_model)
    # Add basin mask
    ds_out.coords["mask"] = xr.DataArray(
        dims=ds_out.raster.dims,
        coords=ds_out.raster.coords,
        data=ds_model["mask"].values,
        attrs=dict(_FillValue=ds_model["mask"].raster.nodata),
    )
    # Add _FillValue to the data attributes
    for dvar in ds_out.data_vars.keys():
        nodata = ds_out[dvar].raster.nodata
        if nodata is not None:
            ds_out[dvar].attrs.update(_FillValue=nodata)
        else:
            ds_out[dvar].attrs.update(_FillValue=-9999.0)
        # Fill with zeros inside mask and keep NaN outside
        ds_out[dvar] = (
            ds_out[dvar].fillna(0).where(ds_out["mask"], ds_out[dvar].raster.nodata)
        )

    return ds_out


def climate_forcing(
    ds: xr.Dataset,
    ds_model: xr.Dataset,
    timestepsecs: int,
    temp_correction: bool = False,
    dem_forcing: xr.DataArray = None,
) -> xr.Dataset:
    """Prepare Delwaq climate fluxes.

    Parameters
    ----------
    ds : xr.Dataset
        Climate data

        * Required variables: variables listed in climate_vars
    ds_model : xr.Dataset
        Dataset containing the model grid

        * Required variables: [mask]

        * Optional variables for tempearture correction: [elevtn]
    timestepsecs: int
        Delwaq model timestep in seconds.
    temp_correction : bool, optional
        If True temperature are corrected using elevation lapse rate,
        by default False.
    dem_forcing : xr.DataArray, default None
        Elevation data source with coverage of entire meteorological forcing domain.
        If temp_correction is True and dem_forcing is provided this is used in
        combination with elevation at model resolution to correct the temperature.

    Returns
    -------
    ds_out : xr.Dataset
        Dataset containing the processed climate data.
    """
    ds_out = xr.Dataset()
    climate_vars = [v for v in ds.data_vars]
    freq = pd.to_timedelta(timestepsecs, unit="s")

    # Start with wind
    if "wind10_u" in climate_vars and "wind10_v" in climate_vars:
        ds_out["wind"] = hydromt.workflows.forcing.wind(
            ds_model,
            wind_u=ds["wind10_u"],
            wind_v=ds["wind10_v"],
            altitude_correction=False,
        )
        climate_vars.remove("wind10_u")
        climate_vars.remove("wind10_v")
    elif "wind" in climate_vars:
        ds_out["wind"] = hydromt.workflows.forcing.wind(
            ds_model,
            wind=ds_out["wind"],
            altitude_correction=False,
        )
        climate_vars.remove("wind")
    # Add other variables
    temp_vars = [v for v in climate_vars if v.startswith("temp")]
    for v in climate_vars:
        if v in temp_vars:
            temp_v = hydromt.workflows.forcing.temp(
                ds[v],
                dem_model=ds_model["elevtn"],
                dem_forcing=dem_forcing,
                lapse_correction=temp_correction,
                logger=logger,
                freq=freq,
                reproj_method="nearest_index",
                lapse_rate=-0.0065,
                resample_kwargs=dict(label="right", closed="right"),
            )
            ds_out[v] = temp_v
        else:
            da_out = ds[v].raster.reproject_like(ds_model, method="nearest_index")
            da_out = hydromt.workflows.forcing.resample_time(
                da_out,
                freq,
                upsampling="bfill",  # we assume right labeled original data
                downsampling="mean",
                conserve_mass=False,
                logger=logger,
            )
            ds_out[v] = da_out
    # Add basin mask
    ds_out.coords["mask"] = xr.DataArray(
        dims=ds_out.raster.dims,
        coords=ds_out.raster.coords,
        data=ds_model["mask"].values,
        attrs=dict(_FillValue=ds_model["mask"].raster.nodata),
    )
    # Add _FillValue to the data attributes
    for dvar in ds_out.data_vars.keys():
        nodata = ds_out[dvar].raster.nodata
        if nodata is not None:
            ds_out[dvar].attrs.update(_FillValue=nodata)
        else:
            ds_out[dvar].attrs.update(_FillValue=-9999.0)
        # Fill with zeros inside mask and keep NaN outside
        ds_out[dvar] = (
            ds_out[dvar].fillna(0).where(ds_out["mask"], ds_out[dvar].raster.nodata)
        )

    return ds_out
