.. currentmodule:: hydromt_delwaq

.. _em_model:

===========================
D-Emission model components
===========================

The HydroMT-delwaq plugin helps you preparing or updating several components of a D-Emissions model such as emission information or forcing.
The main interactions are available from the HydroMT Command Line Interface and allow you to configure
HydroMT in order to build or update D-Emissions (EM) models.

When building or updating a model from command line a
`model region <https://deltares.github.io/hydromt/latest/user_guide/model_region>`_, a model setup
:ref:`configuration <config_file_EM>` (.ini file) with model components and options and, optionally,
a `data sources <https://deltares.github.io/hydromt/latest/user_guide/data_main>`_ (.yml) file should be prepared.

.. _em_methods:

D-Emissions model setup methods
===============================

An overview of the available D-Emissions model setup components
is provided in the table below. When using HydroMT from the command line only the
setup components are exposed. Click on
a specific method see its documentation.

.. list-table::
    :widths: 20 55
    :header-rows: 1
    :stub-columns: 1

    * - Method
      - Explanation
    * - :py:func:`~DelwaqModel.setup_config`
      - Update config with a dictionary
    * - :py:func:`~DelwaqModel.setup_basemaps`
      - This component sets the model schematization using the hydromodel region and resolution.
    * - :py:func:`~DelwaqModel.setup_monitoring`
      - This component sets the monitoring points and areas.
    * - :py:func:`~DelwaqModel.setup_hydrology_forcing`
      - This component prepares the hydrological fluxes.
    * - :py:func:`~DelwaqModel.setup_emission_raster`
      - This component generates emission related data from a raster data source.
    * - :py:func:`~DelwaqModel.setup_emission_vector`
      - This component generates emission related data from a vector data source.
    * - :py:func:`~DelwaqModel.setup_emission_mapping`
      - This component generates several emission related data based on a table and distributed over related regions (eg administrative boundaries).
    * - :py:func:`~DelwaqModel.setup_roads`
      - This component sets roads network and statistics.


.. _em_files:

D-Emissions datamodel files
===========================

The following table provides an overview of which :py:class:`~hydromt_delwaq.DelwaqModel`
attribute contains which D-Emissions in- and output files. The files are read and written with the associated
read- and write- methods, i.e. :py:func:`~DelwaqModel.read_config`
and :py:func:`~DelwaqModel.write_config` for the
:py:attr:`~DelwaqModel.config`  attribute.


.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - :py:class:`~hydromt_delwaq.DelwaqModel` attribute
     - D-Emissions files
   * - :py:attr:`~hydromt_delwaq.DelwaqModel.config`
     - several `config/*.inc` files linked to emission.inp
   * - :py:attr:`~hydromt_delwaq.DelwaqModel.grid`
     - several static binary data for emission calculation in `staticdata/*.dat`
   * - :py:attr:`~hydromt_delwaq.DelwaqModel.geoms`
     - geometries from the geoms folder (basins.geojson, monareas.geojson etc.)
   * - :py:attr:`~hydromt_delwaq.DelwaqModel.forcing`
     - hydrology.bin
   * - :py:attr:`~hydromt_delwaq.DelwaqModel.hydromaps`
     - additional `hydromodel/*.tif` files with information on the coupled hydrology/hydrodynamic model.
