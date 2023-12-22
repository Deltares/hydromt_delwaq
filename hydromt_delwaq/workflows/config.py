"""Prepare demission and delwaq config dictionnary."""

import logging
from datetime import datetime
from typing import Dict, List

import pandas as pd

logger = logging.getLogger(__name__)


__all__ = ["base_config", "time_config"]


def base_config(
    nrofseg: int,
    nrofexch: int = 0,
    layer_name: str = "EM",
    add_surface: bool = False,
    boundaries: List[str] = None,
    boundaries_type: List[str] = None,
    fluxes: List[str] = None,
) -> Dict:
    """
    Prepare base config dictionnary.

    Files concerned:
    - B3_nrofseg
    - B3_attributes
    - B4_nrofexch
    - B5_boundlist
    - B7_surf (optional)
    - B7_fluxes (optional)

    Parameters
    ----------
    nrofseg : int
        Number of segments.
    nrofexch : int, optional
        Number of exchanges, by default 0 (for demission).
    layer_name : str, optional
        Layer name for example for surface water or emission substance, by default "EM".
    add_surface : bool, optional
        Add surface config, by default False.
    boundaries : List[str], optional
        List of boundaries IDs (numbers). If None, no boundary data is added to
        B5_boundlist. By default None.
    boundaries_type : List[str], optional
        List of boundaries types matching the order in ``boundaries``. By default None.
    fluxes : List[str], optional
        List of fluxes to add to B7_fluxes if fluxes is not None. By default None.

    Returns
    -------
    config : Dict
        Base Config dictionnary.
    """
    config = dict()

    # B3_nrofseg
    config["B3_nrofseg"] = {
        "l1": f"{nrofseg} ; nr of segments",
    }

    # B3_attributes
    config["B3_attributes"] = {
        "l1": "      ; DELWAQ_COMPLETE_ATTRIBUTES",
        "l2": " 1    ; one block with input",
        "l3": " 2    ; number of attributes, they are :",
        "l4": "     1     2",
        "l5": " 1    ; file option in this file",
        "l6": " 1    ; option without defaults",
        "l7": f"     {int(nrofseg)}*01 ; {layer_name}",
        "l8": " 0    ; no time dependent attributes",
    }

    # B4_nrofexch
    config["B4_nrofexch"] = {
        "l1": f"{nrofexch} 0 0 ; x, y, z direction",
    }

    # B5_boundlist
    config["B5_boundlist"] = {
        "l1": ";'NodeID' 'Number' 'Type'",
    }
    if boundaries is not None:
        for i in range(len(boundaries)):
            lstr = "l" + str(i + 2)
            bd_id = int(boundaries[i])
            bd_type = boundaries_type[i]
            config["B5_boundlist"].update({lstr: f"'BD_{bd_id}' '{bd_id}' '{bd_type}'"})

    # B7_surf
    if add_surface:
        config["B7_surf"] = {
            "l1": f"PARAMETERS Surf ALL DATA {nrofseg}*1.0",
        }

    # B7_fluxes
    if fluxes is not None:
        config["B7_fluxes"] = {
            "l1": "SEG_FUNCTIONS",
            "l2": " ".join(fluxes),
        }

    return config


def time_config(
    starttime: str,
    endtime: str,
    timestepsecs: int,
    sysclock_format: str = "seconds",
) -> Dict:
    """
    Prepare time config dictionnary.

    Files concerned:
    - B1_timestamp
    - B2_outputtimes
    - B2_sysclock
    - B2_timers
    - B2_timers_only

    Parameters
    ----------
    starttime : str
        Start time in format "yyyy-mm-dd hh:mm:ss".
    endtime : str
        End time in format "yyyy-mm-dd hh:mm:ss".
    timestepsecs : int
        Timestep in seconds.
    sysclock_format : str, optional
        Format of the system clock, either "seconds" or "days", by default "seconds".

    Returns
    -------
    config : Dict
        Time Config dictionnary.
    """
    config = dict()

    # Add time info to config
    ST = datetime.strptime(starttime, "%Y-%m-%d %H:%M:%S")
    ET = datetime.strptime(endtime, "%Y-%m-%d %H:%M:%S")

    # B1_timestamp
    config["B1_timestamp"] = {
        "l1": f"'T0: {ST.strftime('%Y.%m.%d %H:%M:%S')}  (scu=       1s)'",
    }

    # B2_outputtimes
    STstr = ST.strftime("%Y/%m/%d-%H:%M:%S")
    ETstr = ET.strftime("%Y/%m/%d-%H:%M:%S")
    # Get timedelta from timestepsecs
    timestep = pd.Timedelta(seconds=timestepsecs)
    hours = int(timestep.seconds / 3600)
    minutes = int(timestep.seconds / 60)
    seconds = int(timestep.seconds - minutes * 60)
    minutes -= hours * 60
    timestepstring = "%03d%02d%02d%02d" % (timestep.days, hours, minutes, seconds)
    config["B2_outputtimes"] = {
        "l1": f"  {STstr}  {ETstr}  {timestepstring} ; mon/bal",
        "l2": f"  {STstr}  {ETstr}  {timestepstring} ; map",
        "l3": f"  {STstr}  {ETstr}  {timestepstring} ; his",
    }

    # B2_sysclock
    if sysclock_format == "seconds":
        timestepstring = f"{timestepsecs:7d}"
        # timestepsec = timestep.days * 86400 + timestep.seconds
    else:
        # days
        timestepstring = timestep.days
    config["B2_sysclock"] = {
        "l1": f" {timestepstring} 'DDHHMMSS' 'DDHHMMSS'  ; system clock",
    }

    # B2_timers
    config["B2_timers"] = {
        "l1": f"  {STstr} ; start time",
        "l2": f"  {ETstr} ; stop time",
        "l3": "  0 ; timestep constant",
        "l4": "; dddhhmmss format for timestep",
        "l5": f"{timestepstring} ; timestep",
    }

    # B2_timers_only
    config["B2_timers_only"] = {
        "l1": f"  {STstr} ; start time",
        "l2": f"  {ETstr} ; stop time",
        "l3": "  0 ; timestep constant",
    }

    return config
