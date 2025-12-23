.. currentmodule:: hydromt_delwaq

.. _workflows:

Workflows
=========

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

monitoring
----------

.. autosummary::
   :toctree: _generated

   workflows.monitoring.monitoring_points_from_dataarray
   workflows.monitoring.monitoring_points_from_geodataframe
   workflows.monitoring.monitoring_areas

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
   workflows.segments.maps_from_hydromodel
   workflows.segments.geometrymaps
   workflows.segments.pointer
