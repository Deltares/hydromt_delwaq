import xarray as xr
import logging
from xclim.indices.stats import frequency_analysis

import hydromt

logger = logging.getLogger(__name__)

__all__ = [
    "bankfull_volume",
    "bankfull_volume_from_discharge",
]


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
