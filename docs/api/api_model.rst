.. currentmodule:: hydromt_delwaq.delwaq
.. _api_model:

DELWAQ model
============

Initialize
----------

.. autosummary::
   :toctree: ../_generated/

   DelwaqModel

.. _model_components:

Build components
----------------

For DELWAQ, the different components available for building or updating are:

.. autosummary::
   :toctree: ../_generated/

   DelwaqModel.setup_config
   DelwaqModel.setup_basemaps
   DelwaqModel.setup_monitoring
   DelwaqModel.setup_hydrology_forcing
   DelwaqModel.setup_emission_raster
   DelwaqModel.setup_emission_vector
   DelwaqModel.setup_emission_mapping

Model specific attributes
-------------------------

A Delwaq model in HydroMT also has a set of specific attributes, on top of the ones from the model API. These are:

.. autosummary::
   :toctree: ../_generated/

   DelwaqModel.basins
   DelwaqModel.hydromaps

Model specific methods
--------------------------

.. autosummary::
   :toctree: ../_generated/

   DelwaqModel.set_hydromaps

Model specific I/O methods
--------------------------
High level I/O methods

.. autosummary::
   :toctree: ../_generated/

   DelwaqModel.read_hydromaps
   DelwaqModel.write_hydromaps
   DelwaqModel.read_pointer
   DelwaqModel.write_pointer
   DelwaqModel.read_geometry
   DelwaqModel.write_geometry

Intermediate I/O methods

.. autosummary::
   :toctree: ../_generated/

   DelwaqModel.write_monitoring
   DelwaqModel.dw_WriteSegmentOrExchangeData
   DelwaqModel.dw_WriteWaqGeom
