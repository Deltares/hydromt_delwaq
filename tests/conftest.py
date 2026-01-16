"""add global fixtures."""
from os.path import abspath, dirname, join

import pytest

from hydromt_delwaq.delwaq import DelwaqModel
from hydromt_delwaq.demission import DemissionModel

TESTDATADIR = join(dirname(abspath(__file__)), "data")
EXAMPLEDIR = join(dirname(abspath(__file__)), "..", "examples")


@pytest.fixture()
def example_delwaq_model():
    root = join(EXAMPLEDIR, "WQ_piave")
    mod = DelwaqModel(
        root=root,
        mode="r",
        data_libs=["artifact_data"],
    )
    return mod


@pytest.fixture()
def example_demission_model():
    root = join(EXAMPLEDIR, "EM_piave")
    mod = DemissionModel(
        root=root,
        mode="r",
        data_libs=["artifact_data"],
    )
    return mod
