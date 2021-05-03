.. currentmodule:: hydromt_delwaq.delwaq
.. _api_model:

DELWAQ model
============

Initialize
----------

.. autosummary::
   :toctree: ../generated/

   DelwaqModel

Build components
----------------

.. autosummary::
   :toctree: ../generated/

   DelwaqModel.setup_config
   DelwaqModel.setup_basemaps
   DelwaqModel.setup_monitoring
   DelwaqModel.setup_hydrology_forcing
   DelwaqModel.setup_emission_raster
   DelwaqModel.setup_emission_vector
   DelwaqModel.setup_emission_mapping

Model specific attributes
-------------------------

.. autosummary::
   :toctree: ../generated/

   DelwaqModel.basins
   DelwaqModel.hydromaps

Model specific methods
--------------------------

.. autosummary::
   :toctree: ../generated/

   DelwaqModel.set_hydromaps

Model specific I/O methods
--------------------------
High level I/O methods

.. autosummary::
   :toctree: ../generated/

   DelwaqModel.read_hydromaps
   DelwaqModel.write_hydromaps
   DelwaqModel.read_pointer
   DelwaqModel.write_pointer
   DelwaqModel.read_geometry
   DelwaqModel.write_geometry

Intermediate I/O methods

.. autosummary::
   :toctree: ../generated/

   DelwaqModel.write_monitoring
   DelwaqModel.dw_WriteSegmentOrExchangeData
   DelwaqModel.dw_WriteWaqGeom
