"""Unit tests for hydromt_delwaq methods and workflows."""

from os.path import abspath, dirname, isfile, join

import hydromt
import numpy as np
import xarray as xr

TESTDATADIR = join(dirname(abspath(__file__)), "data")


def test_setup_grid(example_delwaq_model):
    # Initialize model and read results
    mod = example_delwaq_model

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


def test_setup_3dgrid(tmpdir, example_delwaq_model):
    """Test adding 3D grid to model and writing."""
    # 3D grid
    grid_fn = join(TESTDATADIR, "INM_INM-CM5-0_ssp585_far.nc")
    grid = xr.open_dataset(grid_fn, mask_and_scale=False).squeeze()

    # Set root to tmpdir
    example_delwaq_model.read()
    example_delwaq_model.set_root(tmpdir, mode="w")

    # Add to grid
    example_delwaq_model.setup_grid_from_rasterdataset(
        raster_fn=grid,
        variables=["temp"],
        reproject_method="nearest",
        rename={"temp": "temp_INM_INM-CM5-0_ssp585_far"},
    )

    # Checks on grid
    assert "temp_INM_INM-CM5-0_ssp585_far" in example_delwaq_model.grid
    assert (
        example_delwaq_model.grid["temp_INM_INM-CM5-0_ssp585_far"].raster.dim0
        == "month"
    )

    # Write grid
    example_delwaq_model.write_grid(fn="staticdata_CC/{name}.dat")

    # Check on files
    assert isfile(
        join(tmpdir, "staticdata_CC", "temp_INM_INM-CM5-0_ssp585_far_month_1.dat")
    )
