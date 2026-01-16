.. _installation_guide:

==================
Installation guide
==================

Getting started
===============

HydroMT-DELWAQ is a model plugin for `HydroMT <https://deltares.github.io/hydromt>`_, extending its core functionalities with DELWAQ-specific components and workflows.
It can be installed as a standalone package or alongside other HydroMT model plugins (e.g. HydroMT-SFINCS, HydroMT-Fiat).
We recommend installing HydroMT-DELWAQ in a dedicated Python environment to ensure dependency consistency.

Prerequisite: Python installation
=================================

You will need **Python 3.11 or newer** and a package/environment manager such as pip, conda, mamba or uv.
These tools simplify installing packages and managing isolated environments.

If you do not yet have one installed, we recommend either:

- `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_
- `Miniforge (Mambaforge) <https://conda-forge.org/docs/>`_
- `uv <https://docs.astral.sh/uv/>`_

Both conda variants come preconfigured with the **conda-forge** channel, which provides free and open packages used by HydroMT.

Installing HydroMT-DELWAQ
=========================

HydroMT-DELWAQ is available from both **PyPI** and **conda-forge**.
The simplest and most flexible approach is to install it using **pip** inside a new environment.

Installation in a conda environment using uv / pip
--------------------------------------------------

We recommend creating a clean environment to avoid dependency conflicts. For example:

.. code-block:: console

    $ conda create -n hydromt-delwaq uv python=3.13
    $ conda activate hydromt-delwaq
    $ uv pip install hydromt_delwaq

This will install HydroMT-DELWAQ along with HydroMT core, HydroMT-Wflow and all required dependencies.

To verify the installation, you can list the installed HydroMT plugins:

.. code-block:: console

    $(hydromt-delwaq) hydromt --models
        Model plugins:
            - model (hydromt 1.3.0)
            - example_model (hydromt 1.3.0)
            - delwaq (hydromt_delwaq 0.4.0)
            - demission (hydromt_delwaq 0.4.0)
            - wflow_sbm (hydromt_wflow 1.0.0)
            - wflow_sediment (hydromt_wflow 1.0.0)


Installing optional dependencies
--------------------------------

HydroMT-DELWAQ provides several optional dependencies that extend its capabilities, such as running JUpyter Notebook and plotting maps.
You can install these easily using pip's extras syntax:

.. code-block:: console

    $(hydromt-delwaq) uv pip install "hydromt_delwaq[examples]"

This will install optional packages such as:

- **jupyterlab** - enables running Jupyter Notebooks for examples and tutorials.
- **cartopy** - enables advanced geospatial plotting capabilities.


For a list of all the optional dependency groups and their contents, have a look at the
`pyproject.toml` file. Use `hydromt_delwaq[full]` to install all optional dependencies
(including developer's dependencies).


Installing via conda
--------------------

HydroMT-DELWAQ is also available through the conda-forge channel. You can install it directly with:

.. code-block:: console

    $ conda create -n hydromt-delwaq -c conda-forge hydromt_delwaq
    $ conda activate hydromt-delwaq

Developer installation
======================

If you want to contribute to HydroMT-DELWAQ or modify its source code, see the
:ref:`Developer installation guide <dev_env>`.
