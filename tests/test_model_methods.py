"""Unit tests for hydromt_delwaq methods and workflows"""

import pytest
from os.path import join, dirname, abspath
import warnings
import pdb
import numpy as np
from hydromt_delwaq.delwaq import DelwaqModel
import hydromt

import logging

TESTDATADIR = join(dirname(abspath(__file__)), "data")
EXAMPLEDIR = join(dirname(abspath(__file__)), "..", "examples")


def test_setup_grid():
    logger = logging.getLogger(__name__)
    # read model from examples folder
    root = join(EXAMPLEDIR, "WQ_piave")

    # Initialize model and read results
    mod = DelwaqModel(root=root, mode="r", data_libs="artifact_data", logger=logger)

    # Tests on setup_emission_vector
    mod.setup_emission_vector(
        emission_fn="hydro_lakes",
        rasterize_method="fraction",
    )

    assert "hydro_lakes" in mod.grid
    assert np.round(mod.grid["hydro_lakes"].values.max(), 4) == 0.8609

    mod.setup_emission_vector(
        emission_fn="hydro_reservoirs",
        rasterize_method="area",
    )

    gdf_grid = mod.grid.raster.vector_grid()
    crs_utm = hydromt.gis_utils.parse_crs("utm", gdf_grid.to_crs(4326).total_bounds)
    gdf_grid = gdf_grid.to_crs(crs_utm)

    assert "hydro_reservoirs" in mod.grid
    assert mod.grid["hydro_reservoirs"].values.max() <= gdf_grid.area.max()
