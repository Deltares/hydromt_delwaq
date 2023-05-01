"""hydroMT plugin for DELWAQ models."""

from os.path import dirname, join, abspath


__version__ = "0.2.2.dev"

DATADIR = join(dirname(abspath(__file__)), "data")

# Set environment variables (this will be temporary)
# to use shapely 2.0 in favor of pygeos (if installed)
import os

os.environ["USE_PYGEOS"] = "0"

from .delwaq import *
