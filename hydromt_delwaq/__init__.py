"""hydroMT plugin for DELWAQ models."""

from os.path import abspath, dirname, join

__version__ = "0.3.2.dev0"

# Set environment variables (this will be temporary)
# to use shapely 2.0 in favor of pygeos (if installed)
import os

os.environ["USE_PYGEOS"] = "0"

from hydromt_delwaq.delwaq import DelwaqModel
from hydromt_delwaq.demission import DemissionModel

__all__ = ["DelwaqModel", "DemissionModel"]
