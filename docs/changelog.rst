==========
What's new
==========
All notable changes to this project will be documented in this page.

The format is based on `Keep a Changelog`_, and this project adheres to
`Semantic Versioning`_.

v0.4.0 (19 January 2026)
========================
Update for HydroMT v1.3.0 and HydroMT-Wflow v1.0.0 compatibility. The main change is that
the ``DelwaqModel`` and ``DEmissionModel`` are now separated and have a component structure.
HydroMT data catalog format and workflow YAML were also  updated. The code was also adapted
to support Wflow.jl models v1.0.0.

Added
-----
- ``DelwaqConfigComponent`` and ``DEmissionConfigComponent`` to handle the configuration files.
- ``DelwaqStaticdataComponent`` to handle gridded static data.
- ``DelwaqHydromapsComponent`` to handle hydromaps from hydromodel needed for processing by HydroMT.
- ``DelwaqPointerComponent`` to handle pointer properties of the delwaq model.
- ``DEmissionGeometryComponent`` to handle geometry properties of the demission model.
- ``DelwaqForcingComponent`` and ``DEmissionForcingComponent`` to handle hydrological, sediment and climate forcing data.
- Readers for ``pointer``, ``forcing`` (netcdf copy) and ``geometry`` components.

Changed
-------
- ``DelwaqModel`` and ``DEmissionModel`` are now separated classes with a component structure.
- ``grid`` is now ``staticdata``.
- All ``setup_emission_*`` methods are now in the ``DEmissionModel`` class only.
- ``setup_sediment_forcing`` is now in the ``DelwaqModel`` class only.

Fixed
-----
- Integration tests for all components (not just geoms and staticdata).

v0.3.1 (16 January 2026)
========================
Minor fixes and improvements.

Added
-----
- Support for 3D data (example for cyclic data) in ``grid``.
- Flexible names for ``grid`` and ``forcing`` filenames in read/write methods.

Fixed
-----
- Errors in ``setup_roads`` method.
- Fix dimension mismatch for waqgeom file by upgrading xugrid dependency.

v0.3.0 (21 June 2023)
=====================
This release introduces a new separate class to prepare D-Emission models **DemissionModel** and adds several new features to the **DelwaqModel** class.
You will need to re-install hydromt-delwaq to use this new version.

Added
-----
- **DemissionModel** class: specific class to prepare D-Emission models. Can be accessed in the command line using **demission** entrypoint.
- **setup_climate_forcing** method to add climate data to a delwaq or demission model.
- **pointer** property for DelwaqModel to easily access pointer properties such as number of segments, exchanges, fluxes etc.
- **read_config** method to read the config files from the config folder.
- **write_waqgeom** method to write a waqgeom file in order to save Delwaq outputs in NetCDF format (not written by default).
- New workflows to prepare Delwaq config **workflows.config** and forcing **workflows.forcing**.

Changed
-------
- **DelwaqModel** (delwaq) only builds classic D-Water Quality models. For D-Emissions model use the new DemissionModel class (hydromt build demission).
- DelwaqModel (and DemissionModel) are subclasses of the HydroMT core GridModel and core GridModel methods are now available.
- **setup_basemaps** method is now more flexible to include any flux and boundary exchange to surface_water that the user wishes to include rather than a few predefined fluxes only.
- **setup_basemaps** method allows to add any variable from hydromodel to the staticdata grid.
- **setup_monitoring**: for monitoring areas, 'compartments' option is not available anymore (only surface_water) and a new option for 'riverland' to distinguish between river and land cells was added.
- **setup_hydrology_forcing**: added option to support summing delwaq flux and volumes using different variables from the hydrology forcing data source. Eg. surface water volume can be river+lake+reservoir volumes. Check the documentation for more details.

Dependencies
------------

- Added dependency to tqdm for progress bars in the write_forcing method.
- Added dependency to xugrid to produce the waqgeom file.

v0.2.1 (27 April 2023)
======================
Small changes to release on conda-forge and update the example notebooks to hydromt version >= 0.7.0 (region argument is optional and artifact_data not read by default).
The documentation was also update to the new light theme style.

v0.2.0 (21 December 2022)
=========================
We now use rioxarray to read raster data. We recommend reinstalling your hydromt and hydromt_delwaq environment including the rioxarray package.
Following an update in xarray, hydromt version should be >= 0.5.0.

Added
-----

- New setup_sediment_forcing method to add eroded soil particles emissions as input to the river in the WQ component.
  Different particles size can be taken into account.

Changed
-------

- Faster method for fraction calculation in setup_emission_vector
- New area method for setup_emission_vector

Fixed
-----

- Added rioxarray dependency for xr.open_rasterio deprecation

v0.1.2 (6 October 2021)
=======================
This release introduces a new component to prepare road related data and adds several fixes for computation of hydrological forcing.

Added
-----

- Option to add a minimum water volumes to avoid zero volume errrors in setup_hydrology_forcing.
- Road network processing in setup_roads.

Changed
-------

- Rewrite and update of the setup_emission_mapping component.

Fixed
-----

- Fixed the waqgeom file for FEWS.
- Reproject the hydrology forcing to the model grid to avoid inverted latitudes.
- Fixes for volume computations.
- Tests and examples linked to updates in the data catalog names.

Documentation
-------------

- Added changelog.
- Update of the template after updae of sphinx_rtd_theme.

v0.1.1 (6 May 2021)
===================

Added
-----

- Compatibility with wflow Julia.

Fixed
-----

- Compatibility of the plugin on Linux.
- Bugfix error when running setup_emission_raster in update mode.

Documentation
-------------

- Additional notebooks.
- Working Binder environment.

v0.1.0 (4 May 2021)
===================
Initial open source release of HydroMT-delwaq plugin, also published on pypi. Noticeable changes are listed below.

Added
-----

- Plugin working in update mode.

Changed
-------

- Consistent setup fonctions arguments for data sources ('_fn').
- Moving segments creation in a new segments.py workflow.

Documentation
-------------

- Initial version of the documentation on github-pages.
- **Latest** and **stable** version of the documentation.
- Add examples notebooks for the documentation.

Tests
-----

- Initial tests for DELWAQ EM and WQ.

.. _Keep a Changelog: https://keepachangelog.com/en/1.0.0/
.. _Semantic Versioning: https://semver.org/spec/v2.0.0.html
