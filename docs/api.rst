.. currentmodule:: hydromt_delwaq

.. _api_reference:

=============
API reference
=============

.. _api_model:

DELWAQ model class
==================

Initialize
----------

.. autosummary::
   :toctree: _generated

   DelwaqModel

.. _methods:

Setup methods
-------------

.. autosummary::
   :toctree: _generated

   DelwaqModel.setup_config
   DelwaqModel.setup_basemaps
   DelwaqModel.setup_monitoring
   DelwaqModel.setup_hydrology_forcing
   DelwaqModel.setup_sediment_forcing
   DelwaqModel.setup_emission_raster
   DelwaqModel.setup_emission_vector
   DelwaqModel.setup_emission_mapping

Attributes
----------

.. autosummary::
   :toctree: _generated

   DelwaqModel.region
   DelwaqModel.crs
   DelwaqModel.res
   DelwaqModel.root
   DelwaqModel.config
   DelwaqModel.grid
   DelwaqModel.geoms
   DelwaqModel.forcing
   DelwaqModel.basins
   DelwaqModel.hydromaps

High level methods
------------------

.. autosummary::
   :toctree: _generated

   DelwaqModel.read
   DelwaqModel.write
   DelwaqModel.build
   DelwaqModel.update
   DelwaqModel.set_root

General methods
---------------

.. autosummary::
   :toctree: _generated


   DelwaqModel.setup_config
   DelwaqModel.get_config
   DelwaqModel.set_config
   DelwaqModel.read_config
   DelwaqModel.write_config

   DelwaqModel.set_grid
   DelwaqModel.read_grid
   DelwaqModel.write_grid

   DelwaqModel.set_geoms
   DelwaqModel.read_geoms
   DelwaqModel.write_geoms

   DelwaqModel.set_forcing
   DelwaqModel.read_forcing
   DelwaqModel.write_forcing

   DelwaqModel.set_hydromaps
   DelwaqModel.read_hydromaps
   DelwaqModel.write_hydromaps

.. _workflows:

DELWAQ workflows
================

.. autosummary::
   :toctree: _generated

   workflows.emissions
   workflows.segments
   workflows.roads
