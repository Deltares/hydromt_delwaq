What's new
==========
All notable changes to this project will be documented in this page.

The format is based on `Keep a Changelog`_, and this project adheres to
`Semantic Versioning`_.

[Unreleased]
------------

Added
^^^^^

- Option to correct the water volumes to minimize water balances errors in setup_hydrology_forcing.

Fixed
^^^^^

- Fixed the waqgeom file for FEWS.
- Reproject the hydrology forcing to the model grid to avoid inverted latitudes.

Documentation
^^^^^^^^^^^^^

- Added changelog.

v0.1.1 (6 May 2021)
-------------------

Added
^^^^^

- Compatibility with wflow Julia.

Fixed
^^^^^

- Compatibility of the plugin on Linux.
- Bugfix error when running setup_emission_raster in update mode.

Documentation
^^^^^^^^^^^^^

- Additional notebooks.
- Working Binder environment.

v0.1.0 (4 May 2021)
----------------------
Initial open source release of hydroMT delwaq plugin, also published on pypi. Noticeable changes are listed below.

Added
^^^^^

- Plugin working in update mode.

Changed
^^^^^^^

- Consistent setup fonctions arguments for data sources ('_fn').
- Moving segments creation in a new segments.py workflow.

Documentation
^^^^^^^^^^^^^

- Initial version of the documentation on github-pages.
- **Latest** and **stable** version of the documentation.
- Add examples notebooks for the documentation.

Tests
^^^^^

- Initial tests for delwaq EM and WQ.

.. _Keep a Changelog: https://keepachangelog.com/en/1.0.0/
.. _Semantic Versioning: https://semver.org/spec/v2.0.0.html