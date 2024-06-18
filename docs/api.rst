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
   DelwaqModel.setup_climate_forcing
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
   DelwaqModel.pointer
   DelwaqModel.nrofseg
   DelwaqModel.nrofexch
   DelwaqModel.surface_water
   DelwaqModel.fluxes

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

   DelwaqModel.set_pointer
   DelwaqModel.read_pointer

   DelwaqModel.write_waqgeom

.. _api_model_emission:

D-Emission model class
==================

Initialize
----------

.. autosummary::
   :toctree: _generated

   DemissionModel

.. _methods:

Setup methods
-------------

.. autosummary::
   :toctree: _generated

   DemissionModel.setup_config
   DemissionModel.setup_basemaps
   DemissionModel.setup_monitoring
   DemissionModel.setup_hydrology_forcing
   DemissionModel.setup_climate_forcing
   DemissionModel.setup_emission_raster
   DemissionModel.setup_emission_vector
   DemissionModel.setup_emission_mapping
   DemissionModel.setup_roads

Attributes
----------

.. autosummary::
   :toctree: _generated

   DemissionModel.region
   DemissionModel.crs
   DemissionModel.res
   DemissionModel.root
   DemissionModel.config
   DemissionModel.grid
   DemissionModel.geoms
   DemissionModel.forcing
   DemissionModel.basins
   DemissionModel.hydromaps
   DemissionModel.nrofseg
   DemissionModel.nrofexch
   DemissionModel.fluxes

High level methods
------------------

.. autosummary::
   :toctree: _generated

   DemissionModel.read
   DemissionModel.write
   DemissionModel.build
   DemissionModel.update
   DemissionModel.set_root

General methods
---------------

.. autosummary::
   :toctree: _generated


   DemissionModel.setup_config
   DemissionModel.get_config
   DemissionModel.set_config
   DemissionModel.read_config
   DemissionModel.write_config

   DemissionModel.set_grid
   DemissionModel.read_grid
   DemissionModel.write_grid

   DemissionModel.set_geoms
   DemissionModel.read_geoms
   DemissionModel.write_geoms

   DemissionModel.set_forcing
   DemissionModel.read_forcing
   DemissionModel.write_forcing

   DemissionModel.set_hydromaps
   DemissionModel.read_hydromaps
   DemissionModel.write_hydromaps

   DemissionModel.write_geometry

   DemissionModel.write_waqgeom

.. _workflows:

DELWAQ workflows
================

config
------

.. autosummary::
   :toctree: _generated

   workflows.config.base_config
   workflows.config.time_config

emissions
---------

.. autosummary::
   :toctree: _generated

   workflows.emissions.emission_raster
   workflows.emissions.emission_vector
   workflows.emissions.admin

forcing
-------

.. autosummary::
   :toctree: _generated

   workflows.forcing.hydrology_forcing
   workflows.forcing.hydrology_forcing_em
   workflows.forcing.sediment_forcing
   workflows.forcing.climate_forcing

geometry
--------

.. autosummary::
   :toctree: _generated

   workflows.geometry.compute_geometry

roads
-----

.. autosummary::
   :toctree: _generated

   workflows.roads.roads_emissions_country
   workflows.roads.roads_emissions_segments

segments
--------

.. autosummary::
   :toctree: _generated

   workflows.segments.hydromaps
   workflows.segments.geometrymaps
   workflows.segments.pointer
