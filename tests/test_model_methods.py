"""Unit tests for hydromt_delwaq methods and workflows."""

from os.path import abspath, dirname, isfile, join

import hydromt
import numpy as np
import xarray as xr

from hydromt_delwaq.workflows.emissions import gridarea

TESTDATADIR = join(dirname(abspath(__file__)), "data")


def test_setup_grid(example_demission_model):
    # Initialize model and read results
    mod = example_demission_model

    # Tests on setup_emission_vector
    mod.setup_emission_vector(
        emission_fn="hydro_lakes",
        rasterize_method="fraction",
    )

    assert "hydro_lakes" in mod.staticdata.data
    assert np.round(mod.staticdata.data["hydro_lakes"].values.max(), 4) == 0.8609

    mod.setup_emission_vector(
        emission_fn="hydro_reservoirs",
        rasterize_method="area",
    )

    gdf_grid = mod.staticdata.data.raster.vector_grid()
    crs_utm = hydromt.gis.gis_utils.parse_crs("utm", gdf_grid.to_crs(4326).total_bounds)
    gdf_grid = gdf_grid.to_crs(crs_utm)

    assert "hydro_reservoirs" in mod.staticdata.data
    assert mod.staticdata.data["hydro_reservoirs"].values.max() <= gdf_grid.area.max()


def test_setup_3dgrid(tmpdir, example_delwaq_model):
    """Test adding 3D grid to model and writing."""
    # 3D grid
    grid_fn = join(TESTDATADIR, "INM_INM-CM5-0_ssp585_far.nc")
    grid = xr.open_dataset(grid_fn, mask_and_scale=False).squeeze()

    # Set root to tmpdir
    example_delwaq_model.read()
    example_delwaq_model.root.set(tmpdir, mode="w")

    # Add to grid
    example_delwaq_model.setup_staticdata_from_rasterdataset(
        raster_data=grid,
        variables=["temp"],
        reproject_method="nearest",
        rename={"temp": "temp_INM_INM-CM5-0_ssp585_far"},
    )

    # Checks on grid
    assert "temp_INM_INM-CM5-0_ssp585_far" in example_delwaq_model.staticdata.data
    assert (
        example_delwaq_model.staticdata.data[
            "temp_INM_INM-CM5-0_ssp585_far"
        ].raster.dim0
        == "month"
    )

    # Write grid
    example_delwaq_model.staticdata.write(filename="staticdata_CC/{name}.dat")

    # Check on files
    assert isfile(
        join(tmpdir, "staticdata_CC", "temp_INM_INM-CM5-0_ssp585_far_month_1.dat")
    )


def test_setup_roads(example_demission_model):
    # Initialize model and read results
    mod = example_demission_model

    # Tests on setup_roads
    mod.setup_roads(
        roads_fn="grip_roads",
        highway_list=["1"],
        country_list=["380"],
        country_fn="wb_countries",
    )

    # Check maps and values
    ds = mod.staticdata.data
    assert "hwy_length_sum_country" in ds
    assert "nnhwy_length_sum_country" in ds
    assert "hwy_length" in ds
    assert "nnhwy_length" in ds

    assert len(np.unique(ds["hwy_length_sum_country"].values)) == 2
    assert np.isclose(ds["hwy_length_sum_country"].values.max(), 1247.7513, atol=1e-4)
    assert np.isclose(ds["hwy_length"].values.max(), 7.6356, atol=1e-4)


def test_setup_emission_raster(example_demission_model):
    # Initialize model and read results
    mod = example_demission_model

    # Tests on setup_emission_raster with classfraction
    mod.setup_emission_raster(
        emission_fn="vito_2015",
        scale_method="classfraction",
        classnumber=40,  # agriculture
    )
    assert "vito_2015" in mod.staticdata.data
    da = mod.staticdata.data["vito_2015"]
    assert da.values.max() <= 1.0
    assert np.isclose(da.values.mean(), 0.1458, atol=1e-4)
    assert da.raster.nodata == -9999.0
    assert da.dtype == "float32"

    # Tests on setup_emission_raster with classarea
    mod.setup_emission_raster(
        emission_fn="vito_2015",
        scale_method="classarea",
        classnumber=40,
        output_name="agriculture_area",
    )
    assert "agriculture_area" in mod.staticdata.data
    da = mod.staticdata.data["agriculture_area"]
    assert np.isclose(da.values.mean(), 348199.78, atol=1e-4)
    assert da.values.max() <= gridarea(mod.staticdata.data).values.max()
