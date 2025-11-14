"""Utilities for hydromt_delwaq."""
import logging
import os
from os.path import abspath, dirname, join
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)


__all__ = ["DATADIR", "dw_WriteSegmentOrExchangeData"]

DATADIR = join(dirname(abspath(__file__)), "data")


def dw_WriteSegmentOrExchangeData(
    ttime: int,
    fname: str | Path,
    datablock: np.ndarray,
    boundids: int,
    WriteAscii: bool = True,
    mode: str = "a",
):
    """
    Write a timestep to a segment/exchange data file.

    Either appends to an existing file or creates a new one.

    Input:
        - time - timestep number for this timestep
        - fname - File path of the segment/exchange data file</param>
        - datablock - array with data
        - boundids to write more than 1 block
        - WriteAscii - if True to make a copy in an ascii checkfile
        - mode - {"a", "w"} Force the writting mode, append or overwrite existing
            files.

    """
    # Supress potential NaN values to avoid error (replaced by -999.0)
    if datablock.dtype == np.float64 or datablock.dtype == np.float32:
        datablock[np.isnan(datablock)] = -999.0
    # Convert the array to a 32 bit float
    totareas = datablock
    for i in range(boundids - 1):
        totareas = np.vstack((totareas, datablock))

    artow = np.array(totareas, dtype=np.float32).copy()
    timear = np.array(ttime, dtype=np.int32)

    if os.path.isfile(fname) and mode == "a":  # append to existing file
        fp = open(fname, "ab")
        tstr = timear.tobytes() + artow.tobytes()
        fp.write(tstr)
        if WriteAscii:
            fpa = open(fname + ".asc", "a")
            timear.tofile(fpa, format="%d\t", sep=":")
            artow.tofile(fpa, format="%10.8f", sep="\t")
            fpa.write("\n")
    else:
        fp = open(fname, "wb")
        tstr = timear.tobytes() + artow.tobytes()
        fp.write(tstr)
        if WriteAscii:
            fpa = open(fname + ".asc", "w")
            timear.tofile(fpa, format="%d\t", sep=":")
            artow.tofile(fpa, format="%10.8f", sep="\t")
            fpa.write("\n")

    fp.close()
    if WriteAscii:
        fpa.close()
