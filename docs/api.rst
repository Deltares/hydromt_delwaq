.. currentmodule:: hydromt_delwaq

.. _api_reference:

=============
API reference
=============

.. _api_model:

Delwaq model class
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
   DelwaqModel.setup_roads

Attributes
----------

.. autosummary::
   :toctree: _generated

   DelwaqModel.region
   DelwaqModel.crs
   DelwaqModel.res
   DelwaqModel.root
   DelwaqModel.config
   DelwaqModel.staticmaps
   DelwaqModel.staticgeoms
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

   DelwaqModel.set_staticmaps
   DelwaqModel.read_staticmaps
   DelwaqModel.write_staticmaps

   DelwaqModel.set_staticgeoms
   DelwaqModel.read_staticgeoms
   DelwaqModel.write_staticgeoms

   DelwaqModel.set_forcing
   DelwaqModel.read_forcing
   DelwaqModel.write_forcing

   DelwaqModel.set_hydromaps
   DelwaqModel.read_hydromaps
   DelwaqModel.write_hydromaps

.. _workflows:

Delwaq workflows
================

.. autosummary::
   :toctree: _generated

   workflows.emissions
   workflows.segments
   workflows.roads
