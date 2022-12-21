==========
What's new
==========
All notable changes to this project will be documented in this page.

The format is based on `Keep a Changelog`_, and this project adheres to
`Semantic Versioning`_.

Unreleased
==========

Added
-----

Changed
-------

Fixed
-----

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