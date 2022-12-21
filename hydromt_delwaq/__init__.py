"""hydroMT plugin for DELWAQ models."""

from os.path import dirname, join, abspath


__version__ = "0.2.0"

DATADIR = join(dirname(abspath(__file__)), "data")

try:
    import pygeos
    import geopandas as gpd

    gpd.options.use_pygeos = True
except ImportError:
    pass

from .delwaq import *
