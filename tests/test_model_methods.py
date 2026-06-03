"""Unit tests for hydromt_delwaq methods and workflows."""

from os.path import abspath, dirname, isfile, join

import hydromt
import numpy as np
import xarray as xr

from hydromt_delwaq.workflows.emissions import gridarea

TESTDATADIR = join(dirname(abspath(__file__)), "data")


def _make_grid(crs, x0, y0, res, nx=4, ny=3):
    x = x0 + (np.arange(nx) + 0.5) * res
    y = y0 - (np.arange(ny) + 0.5) * res  # N -> S
    da = xr.DataArray(
        np.ones((ny, nx), dtype="float32"),
        coords={"y": y, "x": x},
        dims=("y", "x"),
        name="dummy",
    )
    da.raster.set_crs(crs)
    return da


def test_gridarea_projected_returns_uniform_positive_area():
    """Projected CRS (e.g. SVY21 / EPSG:3414) must produce uniform positive areas.

    Regression test for the gridarea bug that fed SVY21 Easting metres into
    ``_reggrid_area`` (which interprets them as longitude degrees). At
    Easting=10500-18000 m, ``sin(radians(Easting))`` is a periodic function
    of large angles and produces partly-negative garbage, which propagated
    through the mm -> m3/s unit conversion and yielded negative precip values
    in dynamicdata.nc. The patched gridarea delegates to
    ``ds.raster.area_grid()`` for projected CRS, which uses res*res*ucf**2.
    """
    res = 30.0
    da = _make_grid("EPSG:3414", x0=10490.0, y0=39280.0, res=res)

    area = gridarea(da)

    assert area.shape == da.shape
    assert float(area.min()) > 0, "projected cell area must be strictly positive"
    expected = res * res  # SVY21 linear unit factor is 1.0 (metres)
    np.testing.assert_allclose(area.values, expected, rtol=1e-6)


def test_gridarea_geographic_path_unchanged():
    """Geographic CRS path must still go through _reggrid_area (spherical cap)."""
    da = _make_grid("EPSG:4326", x0=4.0, y0=52.0, res=0.01)

    area = gridarea(da)

    assert float(area.min()) > 0
    # At ~52N, 0.01deg cells are ~700-800 m on a side -> ~5-7e5 m2.
    assert 4e5 < float(area.mean()) < 9e5


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
