"""Test plugin model class against hydromt.models.model_api."""

import warnings
from os.path import abspath, dirname, join

import numpy as np
import pytest
from hydromt.readers import read_workflow_yaml

from hydromt_delwaq.delwaq import DelwaqModel
from hydromt_delwaq.demission import DemissionModel

TESTDATADIR = join(dirname(abspath(__file__)), "data")
EXAMPLEDIR = join(dirname(abspath(__file__)), "..", "examples")

_models = {
    "EM": {
        "example": "EM_piave",
        "ini": "delwaq_build_EM.yml",
        "ini_update": "delwaq_update_EM_forcing.yml",
        "model_type": DemissionModel,
    },
    "WQ": {
        "example": "WQ_piave",
        "ini": "delwaq_build_WQ.yml",
        "ini_update": "delwaq_update_WQ_sed.yml",
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

    if model == "EM":
        model_init = {"data_libs": ["artifact_data"]}
    elif model == "WQ":
        model_init = {
            "data_libs": ["artifact_data", join(EXAMPLEDIR, "local_sources.yml")]
        }

    # Build model
    model_type = _model["model_type"]
    mod1 = model_type(root=root, mode="w", **model_init)
    mod1.build(steps=steps)

    # Compare with model from examples folder
    root0 = join(EXAMPLEDIR, _model["example"])
    mod0 = model_type(root=root0, mode="r")
    mod0.read()
    mod1 = model_type(root=root, mode="r")
    mod1.read()
    # check maps
    invalid_maps = []
    if len(mod0.staticdata.data) > 0:
        maps = mod0.staticdata.data.raster.vars
        assert mod0.crs == mod1.crs, "map crs mismatch"
        for name in maps:
            map0 = mod0.staticdata.data[name].fillna(0)
            if name not in mod1.staticdata.data:
                invalid_maps.append(f"{name} missing")
            else:
                map1 = mod1.staticdata.data[name].fillna(0)
                if not np.allclose(map0, map1):
                    invalid_maps.append(f"{name} changed")
    invalid_map_str = ", ".join(invalid_maps)
    assert len(invalid_maps) == 0, f"invalid maps: {invalid_map_str}"
    # check geoms
    if mod0.geoms.data:
        for name in mod0.geoms.data:
            geom0 = mod0.geoms.data[name]
            assert name in mod1.geoms.data, f"geom {name} missing"
            geom1 = mod1.geoms.data[name]
            assert geom0.index.size == geom1.index.size, f"geom index {name}"
            assert np.all(geom0.index == geom1.index), f"geom index {name}"
            assert geom0.columns.size == geom1.columns.size, f"geom columns {name}"
            assert np.all(geom0.columns == geom1.columns), f"geom columns {name}"
            assert geom0.crs == geom1.crs, f"geom crs {name}"
            if not np.all(geom0.geometry == geom1.geometry):
                warnings.warn(f"New geom {name} different than the example one.")


@pytest.mark.parametrize("model", list(_models.keys()))
def test_model_update(tmpdir, model):
    # test update method
    # just sees if the CLI run
    # TODO: also check results
    _model = _models[model]
    root = str(tmpdir.join(model + "_update"))
    config = join(TESTDATADIR, _model["ini_update"])
    _, model_init, steps = read_workflow_yaml(config)
    if len(model_init) > 0:
        if "data_libs" in model_init:
            cats = model_init["data_libs"]
            cats.extend(["artifact_data"])
            model_init["data_libs"] = cats
    else:
        model_init = {"data_libs": ["artifact_data"]}
    delwaq_path = join(EXAMPLEDIR, _model["example"])
    # Transform the path to be processed by CLI runner and json.load
    delwaq_path = str(delwaq_path).replace("\\", "/")
    # Update model
    mod1 = _model["model_type"](root=delwaq_path, mode="r", **model_init)
    mod1.update(model_out=root, steps=steps)
    # Check if model is api compliant
    non_compliant_list = mod1._test_model_api()
    assert len(non_compliant_list) == 0
