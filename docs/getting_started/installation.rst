.. _installation_guide:

==================
Installation guide
==================

Prerequisites
=============

For the necessary prerequisites, see `HydroMT core prerequisites <https://deltares.github.io/hydromt/latest/getting_started/installation.html#prerequisites>`_.

Additional dependencies
-----------------------

The HydroMT-delwaq Python package makes use of the HydroMT core.
For a complete list of dependencies, see the pyproject.toml file. 

DELWAQ models are built on top of existing hydrology/hydrodynamic models. 
To build a DELWAQ model, the plugin of the corresponding hydrology/hydrodynamic model 
also needs to be installed (example HydroMT-Wflow).

Installation
============

HydroMT-delwaq is available from pypi and in the near future also installation from conda-forge will be supported.

Install HydroMT-delwaq in a new environment
-------------------------------------------
.. Tip::

    This is our recommended way of installing HydroMT!

To create a new environment called `hydromt-delwaq` from conda-forge with HydroMT and for example HydroMT-Wflow plugin installed, do:

.. code-block:: console

    $ conda create -n hydromt-delwaq -c conda-forge hydromt hydromt_wflow

Then activate your new environment and install HydroMT-delwaq using pip.

.. code-block:: console

    $ conda activate hydromt-delwaq
    $ pip install hydromt_delwaq

Install HydroMT-delwaq in an existing environment
-------------------------------------------------
To install HydroMT Wflow **using mamba or conda** execute the command below after activating the correct environment. 
Note that if some dependencies are not installed from conda-forge the installation may fail. As HydroMT-delwaq is not 
available yet from conda-forge, we recommend to first install HydroMT and the required hydrology/hydrodynamic plugins:

.. code-block:: console

    $ conda install -c conda-forge hydromt hydromt_wflow
    $ pip install hydromt_delwaq

.. Note::

    HydroMT-delwaq is available from pypi and we are working on adding a release from conda-forge (ongoing).

Developer install
-----------------

To be able to test and develop the HydroMT-delwaq package see instructions in the :ref:`Developer installation guide <dev_env>`.

For more information about how to contribute, see `HydroMT contributing guidelines <https://hydromt.readthedocs.io/en/latest/contributing.html>`_.
