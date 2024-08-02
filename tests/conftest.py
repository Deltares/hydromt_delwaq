"""add global fixtures."""
import logging
from os.path import abspath, dirname, join

import pytest

from hydromt_delwaq.delwaq import DelwaqModel

TESTDATADIR = join(dirname(abspath(__file__)), "data")
EXAMPLEDIR = join(dirname(abspath(__file__)), "..", "examples")


@pytest.fixture()
def example_delwaq_model():
    logger = logging.getLogger(__name__)
    root = join(EXAMPLEDIR, "WQ_piave")
    mod = DelwaqModel(
        root=root,
        mode="r",
        data_libs=["artifact_data"],
        logger=logger,
    )
    return mod
