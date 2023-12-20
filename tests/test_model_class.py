"""Test plugin model class against hydromt.models.model_api"""

import pytest
from os.path import join, dirname, abspath
import numpy as np
import pdb
import warnings
from click.testing import CliRunner

import hydromt
from hydromt_delwaq.delwaq import DelwaqModel
from hydromt.cli.cli_utils import parse_config
from hydromt.cli.main import main as hydromt_cli

TESTDATADIR = join(dirname(abspath(__file__)), "data")
EXAMPLEDIR = join(dirname(abspath(__file__)), "..", "examples")

_models = {
    "EM": {
        "example": "EM_piave",
        "ini": "delwaq_build_EM.ini",
        "ini_update": "delwaq_update_EM_forcing.ini",
    },
    "WQ": {
        "example": "WQ_piave",
        "ini": "delwaq_build_WQ.ini",
        "ini_update": "delwaq_update_WQ_sed.ini",
    },
}


@pytest.mark.parametrize("model", list(_models.keys()))
def test_model_class(model):
    _model = _models[model]
    # read model in examples folder
    root = join(EXAMPLEDIR, _model["example"])
    mod = DelwaqModel(root=root, mode="r")
    mod.read()
    # run test_model_api() method
    non_compliant_list = mod.test_model_api()
    assert len(non_compliant_list) == 0
    # pass


@pytest.mark.parametrize("model", list(_models.keys()))
def test_model_build(tmpdir, model):
    # test build method
    # compare results with model from examples folder
    _model = _models[model]
    root = str(tmpdir.join(model))
    config = join(EXAMPLEDIR, _model["ini"])
    wflow_path = join(EXAMPLEDIR, "wflow_piave")
    # Transform the path to be processed by CLI runner and json.load
    wflow_path = str(wflow_path).replace("\\", "/")
    region = "{'wflow': " + f"'{wflow_path}'" + "}"
    # region = "{'wflow': '../examples/wflow_piave'}"
    # Build model
    r = CliRunner().invoke(
        hydromt_cli,
        [
            "build",
            "delwaq",
            root,
            "-r",
            region,
            "-i",
            config,
            "-d",
            "artifact_data",
            "-vv",
        ],
    )
    assert r.exit_code == 0

    # Compare with model from examples folder
    root0 = join(EXAMPLEDIR, _model["example"])
    mod0 = DelwaqModel(root=root0, mode="r")
    mod0.read()
    mod1 = DelwaqModel(root=root, mode="r")
    mod1.read()
    # check maps
    invalid_maps = []
    if len(mod0._grid) > 0:
        maps = mod0.grid.raster.vars
        assert np.all(mod0.crs == mod1.crs), f"map crs mismatch"
        for name in maps:
            map0 = mod0.grid[name].fillna(0)
            map1 = mod1.grid[name].fillna(0)
            if not np.allclose(map0, map1):
                invalid_maps.append(name)
    invalid_map_str = ", ".join(invalid_maps)
    assert len(invalid_maps) == 0, f"invalid maps: {invalid_map_str}"
    # check geoms
    if mod0._geoms:
        for name in mod0.geoms:
            geom0 = mod0.geoms[name]
            geom1 = mod1.geoms[name]
            assert geom0.index.size == geom1.index.size and np.all(
                geom0.index == geom1.index
            ), f"geom index {name}"
            assert geom0.columns.size == geom1.columns.size and np.all(
                geom0.columns == geom1.columns
            ), f"geom columns {name}"
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
    delwaq_path = join(EXAMPLEDIR, _model["example"])
    # Transform the path to be processed by CLI runner and json.load
    delwaq_path = str(delwaq_path).replace("\\", "/")
    # Update model
    r = CliRunner().invoke(
        hydromt_cli,
        [
            "update",
            "delwaq",
            delwaq_path,
            "-o",
            root,
            "-i",
            config,
            "-d",
            "artifact_data",
            "-vv",
        ],
    )
    assert r.exit_code == 0
