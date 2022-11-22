"""Unit tests for hydromt_delwaq methods and workflows"""

import pytest
from os.path import join, dirname, abspath
import warnings
import pdb
import numpy as np
from hydromt_delwaq.delwaq import DelwaqModel

import logging

TESTDATADIR = join(dirname(abspath(__file__)), "data")
EXAMPLEDIR = join(dirname(abspath(__file__)), "..", "examples")


def test_setup_staticmaps():
    logger = logging.getLogger(__name__)
    # read model from examples folder
    root = join(EXAMPLEDIR, "WQ_piave")

    # Initialize model and read results
    mod = DelwaqModel(root=root, mode="r", data_libs="artifact_data", logger=logger)

    # Tests on setup_staticmaps_from_raster
    mod.setup_emission_vector(
        emission_fn="hydro_lakes",
        rasterize_method="fraction",
    )

    assert "hydro_lakes" in mod.staticmaps
    assert np.round(mod.staticmaps["hydro_lakes"].values.max(), 4) == 0.8608
