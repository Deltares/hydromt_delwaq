"""hydroMT plugin for DELWAQ models."""

from os.path import abspath, dirname, join

__version__ = "0.3.0"

DATADIR = join(dirname(abspath(__file__)), "data")

# Set environment variables (this will be temporary)
# to use shapely 2.0 in favor of pygeos (if installed)
import os

os.environ["USE_PYGEOS"] = "0"

from .delwaq import DelwaqModel
from .demission import DemissionModel
