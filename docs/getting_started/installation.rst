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

HydroMT-delwaq is available from pypi and conda-forge. We recommend installing using mamba from conda-forge in a new environment.

.. Note::

    In the commands below you can exchange `mamba` for `conda`, see
    `here <https://deltares.github.io/hydromt/latest/getting_started/installation.html#installation-guide>`_
    for the difference between both.

Install HydroMT-delwaq in a new environment
-------------------------------------------
You can install HydroMT-delwaq in a new environment called `hydromt-delwaq` together with
all optional (see above) and a few additional dependencies with:

.. code-block:: console

  $ mamba env create -f https://raw.githubusercontent.com/Deltares/hydromt_delwaq/main/environment.yml

Then, activate the environment (as stated by mamba/conda) to start making use of HydroMT-Delwaq:

.. code-block:: console

  conda activate hydromt-delwaq

.. Tip::

    If you already have this environment with this name either remove it with
    `conda env remove -n hydromt-delwaq` **or** set a new name for the environment
    by adding `-n <name>` to the line below.

Install HydroMT-delwaq in an existing environment
-------------------------------------------------
To install HydroMT-delwaq in an existing environment execute the command below
where you replace `<environment_name>` with the name of the existing environment.
Note that if some dependencies are not installed from conda-forge but from other
channels the installation may fail.

.. code-block:: console

   $ mamba install -c conda-forge hydromt_delwaq -n <environment_name>

Developer install
-----------------

To be able to test and develop the HydroMT-delwaq package see instructions in the :ref:`Developer installation guide <dev_env>`.

For more information about how to contribute, see `HydroMT contributing guidelines <https://hydromt.readthedocs.io/en/latest/contributing.html>`_.
