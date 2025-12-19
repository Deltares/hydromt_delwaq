"""Test plugin model class against hydromt.models.model_api."""

from os.path import abspath, dirname, join

import pytest
import xarray as xr
from hydromt.readers import read_workflow_yaml
from hydromt_wflow.components import utils

from hydromt_delwaq.delwaq import DelwaqModel
from hydromt_delwaq.demission import DemissionModel

TESTDATADIR = join(dirname(abspath(__file__)), "data")
EXAMPLEDIR = join(dirname(abspath(__file__)), "..", "examples")

_models = {
    "EM": {
        "example": "EM_piave",
        "ini": "delwaq_build_EM_full.yml",
        "model_type": DemissionModel,
    },
    "WQ": {
        "example": "WQ_piave",
        "ini": "delwaq_build_WQ.yml",
        "model_type": DelwaqModel,
    },
}


@pytest.mark.parametrize("model", list(_models.keys()))
def test_model_build(tmpdir, model):
    # test build method
    # compare results with model from examples folder
    _model = _models[model]

    # Root
    root = str(tmpdir.join(model))
    # Config
    config = join(EXAMPLEDIR, _model["ini"])
    _, _, steps = read_workflow_yaml(config)

    # Add local sources for wflow outputs
    model_init = {"data_libs": ["artifact_data", join(EXAMPLEDIR, "local_sources.yml")]}

    # Build model
    model_type = _model["model_type"]
    mod1 = model_type(root=root, mode="w", **model_init)
    mod1.build(steps=steps)
    if model == "WQ":
        # Re-write forcing to also get netcdf copy
        mod1.forcing.write(write_nc=True)

    # Compare with model from examples folder
    root0 = join(EXAMPLEDIR, _model["example"])
    mod0 = model_type(root=root0, mode="r")
    mod0.read()
    mod1 = model_type(root=root, mode="r")
    mod1.read()
    # compare models
    eq, errors = mod1.test_equal(mod0)
    assert eq, f"models not equal: {errors}"


def test_model_update_WQ_forcing(tmpdir):
    # test update method
    root = str(tmpdir.join("WQ_update"))
    config = join(TESTDATADIR, "delwaq_update_WQ_sed.yml")
    _, model_init, steps = read_workflow_yaml(config)
    if len(model_init) > 0:
        if "data_libs" in model_init:
            cats = model_init["data_libs"]
            cats.extend(["artifact_data"])
            model_init["data_libs"] = cats
    else:
        model_init = {"data_libs": ["artifact_data"]}
    delwaq_path = join(EXAMPLEDIR, "WQ_piave")
    # Transform the path to be processed by CLI runner and json.load
    delwaq_path = str(delwaq_path).replace("\\", "/")
    # Update model
    mod1 = DelwaqModel(root=delwaq_path, mode="r", **model_init)
    mod1.update(model_out=root, steps=steps)

    # Check the created forcing
    mod1 = DelwaqModel(root=root, mode="r")
    mod1.read()

    # Check that sediment was added in config
    assert "B7_sediment" in mod1.config.data
    assert mod1.config.data["B7_sediment"] == {
        "l1": "SEG_FUNCTIONS",
        "l2": "ErodIM1 ErodIM2 ErodIM3",
    }

    # Check the forcing file
    forcing_test = xr.open_dataset(
        join(TESTDATADIR, "dynamicdata_WQ.nc"),
        mask_and_scale=False,
    )
    eq, errors = utils.test_equal_grid_data(mod1.forcing.data, forcing_test)

    assert eq, f"Forcing data not equal: {errors}"
