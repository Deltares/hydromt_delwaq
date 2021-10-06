"""Test plugin model class against hydromt.models.model_api"""

import pytest
from os.path import join, dirname, abspath
import numpy as np
import pdb
import warnings
from click.testing import CliRunner

import hydromt
from hydromt.models import MODELS
from hydromt.cli.cli_utils import parse_config
from hydromt.cli.main import main as hydromt_cli

TESTDATADIR = join(dirname(abspath(__file__)), "data")
EXAMPLEDIR = join(dirname(abspath(__file__)), "..", "examples")

_models = {
    "EM": {
        "example": "EM_piave",
        "ini": "delwaq_build_EM.ini",
    },
    "WQ": {
        "example": "WQ_piave",
        "ini": "delwaq_build_WQ.ini",
    },
}


@pytest.mark.parametrize("model", list(_models.keys()))
def test_model_class(model):
    _model = _models[model]
    # read model in examples folder
    root = join(EXAMPLEDIR, _model["example"])
    mod = MODELS.get("delwaq")(root=root, mode="r")
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
        hydromt_cli, ["build", "delwaq", root, region, "-i", config, "-vv"]
    )
    assert r.exit_code == 0

    # Compare with model from examples folder
    root0 = join(EXAMPLEDIR, _model["example"])
    mod0 = MODELS.get("delwaq")(root=root0, mode="r")
    mod0.read()
    mod1 = MODELS.get("delwaq")(root=root, mode="r")
    mod1.read()
    # check maps
    invalid_maps = []
    if len(mod0._staticmaps) > 0:
        maps = mod0.staticmaps.raster.vars
        assert np.all(mod0.crs == mod1.crs), f"map crs mismatch"
        for name in maps:
            map0 = mod0.staticmaps[name].fillna(0)
            map1 = mod1.staticmaps[name].fillna(0)
            if not np.allclose(map0, map1):
                invalid_maps.append(name)
    invalid_map_str = ", ".join(invalid_maps)
    assert len(invalid_maps) == 0, f"invalid maps: {invalid_map_str}"
    # check geoms
    if mod0._staticgeoms:
        for name in mod0.staticgeoms:
            geom0 = mod0.staticgeoms[name]
            geom1 = mod1.staticgeoms[name]
            assert geom0.index.size == geom1.index.size and np.all(
                geom0.index == geom1.index
            ), f"geom index {name}"
            assert geom0.columns.size == geom1.columns.size and np.all(
                geom0.columns == geom1.columns
            ), f"geom columns {name}"
            assert geom0.crs == geom1.crs, f"geom crs {name}"
            if not np.all(geom0.geometry == geom1.geometry):
                warnings.warn(f"New geom {name} different than the example one.")
