.. currentmodule:: hydromt_delwaq

.. _em_model:

===========================
D-Emission model components
===========================

The HydroMT-delwaq plugin helps you preparing or updating several components of a D-Emissions model such as emission information or forcing.
The main interactions are available from the HydroMT Command Line Interface and allow you to configure
HydroMT in order to build or update D-Emissions (demission) models.

When building or updating a model from command line a
`model region <https://deltares.github.io/hydromt/latest/user_guide/model_region>`_, a model setup
:ref:`configuration <config_file_EM>` (.yml file) with model components and options and, optionally,
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
    * - :py:func:`~DemissionModel.setup_config`
      - Update config with a dictionary
    * - :py:func:`~DemissionModel.setup_basemaps`
      - This component sets the model schematization using the hydromodel region and resolution.
    * - :py:func:`~DemissionModel.setup_monitoring`
      - This component sets the monitoring points and areas.
    * - :py:func:`~DemissionModel.setup_hydrology_forcing`
      - This component prepares the hydrological fluxes.
    * - :py:func:`~DemissionModel.setup_climate_forcing`
      - This component prepares the climate forcing (eg for temperature modelling).
    * - :py:func:`~DemissionModel.setup_emission_raster`
      - This component generates emission related data from a raster data source.
    * - :py:func:`~DemissionModel.setup_emission_vector`
      - This component generates emission related data from a vector data source.
    * - :py:func:`~DemissionModel.setup_emission_mapping`
      - This component generates several emission related data based on a table and distributed over related regions (eg administrative boundaries).
    * - :py:func:`~DemissionModel.setup_roads`
      - This component sets roads network and statistics.


.. _em_files:

D-Emissions model components
============================

The following table provides an overview of which :py:class:`~hydromt_delwaq.DemissionModel`
attribute contains which D-Emissions in- and output files. The files are read and written with the associated
read- and write- methods, i.e. :py:func:`~DemissionModel.read_config`
and :py:func:`~DemissionModel.write_config` for the
:py:attr:`~DemissionModel.config`  attribute.


.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - :py:class:`~hydromt_delwaq.DemissionModel` attribute
     - D-Emissions model
   * - :py:attr:`~hydromt_delwaq.DemissionModel.config`
     - several `config/*.inc` files linked to emission.inp
   * - :py:attr:`~hydromt_delwaq.DemissionModel.grid`
     - several static binary data for emission calculation in `staticdata/*.dat`
   * - :py:attr:`~hydromt_delwaq.DemissionModel.geoms`
     - geometries from the geoms folder (basins.geojson, monareas.geojson etc.)
   * - :py:attr:`~hydromt_delwaq.DemissionModel.forcing`
     - hydrology.bin
   * - :py:attr:`~hydromt_delwaq.DemissionModel.hydromaps`
     - additional `hydromodel/*.tif` files with information on the coupled hydrology/hydrodynamic model.
