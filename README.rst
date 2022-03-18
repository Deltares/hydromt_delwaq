HydroMT-delwaq: DELWAQ plugin for HydroMT
#########################################

.. image:: https://codecov.io/gh/Deltares/hydromt_delwaq/branch/main/graph/badge.svg?token=ss3EgmwHhH
    :target: https://codecov.io/gh/Deltares/hydromt_delwaq

.. image:: https://img.shields.io/badge/docs-latest-brightgreen.svg
    :target: https://deltares.github.io/hydromt_delwaq/latest
    :alt: Latest developers docs

.. image:: https://img.shields.io/badge/docs-stable-brightgreen.svg
    :target: https://deltares.github.io/hydromt_delwaq/stable
    :alt: Stable docs last release

.. image:: https://pypip.in/v/hydromt_delwaq/badge.png
    :target: https://pypi.org/project/hydromt_delwaq/
    :alt: Latest PyPI version

.. image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/Deltares/hydromt_delwaq/main?urlpath=lab/tree/examples


`HydroMT <https://github.com/Deltares/hydromt>`_ is a python package, developed by Deltares, to build 
and analyse environmental models. It provides a generic model api with attributes to access the model schematization, 
(dynamic) forcing data, results and states.

This plugin provides an implementation for the `DELWAQ <https://oss.deltares.nl/web/delft3d/delwaq1>`_ water quality engine. 
It details the different steps and explains how to use HydroMT to easily get started and work on your own DELWAQ model. WIth this plugin 
you can interact with both classic **D-Water Quality** models as well as **D-Emission** models.

For detailed information on HydroMT itself, you can visit the `core documentation <https://deltares.github.io/hydromt_plugin/latest/>`_.

Installation
------------

hydroMT-delwaq is availble from pypi and we are working on adding a release from conda-forge (ongoing).

To install hydromt_delwaq using pip do:

.. code-block:: console

  pip install hydromt_delwaq

We recommend installing a hydromt-delwaq environment including the hydromt_delwaq package
based on the environment.yml file. This environment will install all package dependencies 
including the core of hydroMT_ as well as the `hydromt_wflow <https://github.com/Deltares/hydromt_wflow>`_ plugin.

.. code-block:: console

  conda env create -f environment.yml

Documentation
-------------

Learn more about the hydroMT_delwaq plugin in its `online documentation <https://deltares.github.io/hydromt_delwaq/>`_

Contributing
------------

You can find information about contributing to hydroMT at our `Contributing page <https://deltares.github.io/hydromt_plugin/latest/contributing.html>`_.

License
-------

Copyright (c) 2021, Deltares

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General 
Public License as published by the Free Software Foundation, either version 3 of the License, or (at your 
option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the 
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
for more details.

You should have received a copy of the GNU General Public License along with this program. If not, 
see <https://www.gnu.org/licenses/>.
