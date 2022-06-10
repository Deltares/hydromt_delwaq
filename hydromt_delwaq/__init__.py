"""hydroMT plugin for DELWAQ models."""

from os.path import dirname, join, abspath


__version__ = "0.1.3.dev"

DATADIR = join(dirname(abspath(__file__)), "data")

from .delwaq import *
