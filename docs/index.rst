.. _main_sections:

=========================================
HydroMT-DELWAQ: DELWAQ plugin for HydroMT
=========================================

|pypi| |conda_forge| |docs_latest| |docs_stable| |codecov| |license| |binder|

HydroMT_ (Hydro Model Tools) is an open-source Python package that facilitates the process of
building and analyzing spatial geoscientific models with a focus on water system models.
It does so by automating the workflow to go from raw data to a complete model instance which
is ready to run and to analyze model results once the simulation has finished.

This plugin provides an implementation of the model API for DELWAQ including:

- **demission**: emission models, `D-Emission <https://www.deltares.nl/en/software/module/D-Emissions/>`_
  (also referred to as EM in this documentation), that predicts the quantity (fluxes) of substances
  being released from sources of pollutants to the surface waters (mass/time).
- **delwaq** fate and transport models, `D-Water Quality <https://www.deltares.nl/en/software/module/d-water-quality/>`_
  (also referred to as WQ in this documentation), that predicts the fate and transport of substances
  in the surface waters (concentrations in mass / volume).

.. grid:: 2
    :gutter: 2

    .. grid-item-card::
        :text-align: center
        :link: getting_started/intro
        :link-type: doc

        :octicon:`rocket;5em;sd-text-icon blue-icon`
        +++
        **Getting Started**

        First time user? Learn how to install, configure, and start using HydroMT-DELWAQ
        effectively.

    .. grid-item-card::
        :text-align: center
        :link: user_guide/intro
        :link-type: doc

        :octicon:`book;5em;sd-text-icon blue-icon`
        +++
        **User Guide**

        Explore detailed guides on model setup, configuration, and workflows.

    .. grid-item-card::
        :text-align: center
        :link: user_guide/WQ_model
        :link-type: doc

        :octicon:`list-unordered;5em;sd-text-icon blue-icon`
        +++
        **DELWAQ methods**

        Regular user? Here is a quick access to the available methods to prepare a DELWAQ
        model.

    .. grid-item-card::
        :text-align: center
        :link: user_guide/EM_model
        :link-type: doc

        :octicon:`list-unordered;5em;sd-text-icon blue-icon`
        +++
        **D-Emission methods**

        Regular user? Here is a quick access to the available methods to prepare a D-Emission
        model.

    .. grid-item-card::
        :text-align: center
        :link: https://deltares.github.io/hydromt/latest/user_guide/models/model_build.html

        :octicon:`device-desktop;5em;sd-text-icon blue-icon`
        +++
        **Building models with HydroMT**

        Learn more about the basics of model building with HydroMT in the core documentation.

    .. grid-item-card::
        :text-align: center
        :link: https://deltares.github.io/hydromt/latest/user_guide/data_catalog/data_prepare_cat.html

        :octicon:`stack;5em;sd-text-icon blue-icon`
        +++
        **Preparing a Data Catalog**

        Learn more about the HydroMT Data Catalog and how to prepare your own in the
        core documentation.


    .. grid-item-card::
        :text-align: center
        :link: api/index
        :link-type: doc

        :octicon:`code-square;5em;sd-text-icon blue-icon`
        +++
        **API Reference**

        Access the full API documentation for HydroMT-DELWAQ's modules and functions.

    .. grid-item-card::
        :text-align: center
        :link: changelog
        :link-type: doc

        :octicon:`graph;5em;sd-text-icon blue-icon`
        +++
        **What's new?**

        Want to learn what are the new developments and features in the new release?
        Check the changelog.

.. toctree::
   :titlesonly:
   :hidden:

   getting_started/intro.rst
   user_guide/intro.rst
   api/index.rst
   dev/intro.rst

.. _Hydromt: https://deltares.github.io/hydromt/latest/

.. |pypi| image:: https://badge.fury.io/py/hydromt_delwaq.svg
    :target: https://pypi.org/project/hydromt_delwaq/
    :alt: Latest PyPI version

.. |conda_forge| image:: https://anaconda.org/conda-forge/hydromt_delwaq/badges/version.svg
    :target: https://anaconda.org/conda-forge/hydromt_delwaq
    :alt: Conda-Forge

.. |docs_latest| image:: https://img.shields.io/badge/docs-latest-brightgreen.svg
    :target: https://deltares.github.io/hydromt_delwaq/latest
    :alt: Latest developers docs

.. |docs_stable| image:: https://img.shields.io/badge/docs-stable-brightgreen.svg
    :target: https://deltares.github.io/hydromt_delwaq/stable
    :alt: Stable docs last release

.. |codecov| image:: https://codecov.io/gh/Deltares/hydromt_delwaq/branch/main/graph/badge.svg?token=ss3EgmwHhH
    :target: https://codecov.io/gh/Deltares/hydromt_delwaq

.. |license| image:: https://img.shields.io/github/license/Deltares/hydromt_delwaq?style=flat
    :alt: License
    :target: https://github.com/Deltares/hydromt_delwaq/blob/main/LICENSE

.. |binder| image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/Deltares/hydromt_delwaq/main?urlpath=lab/tree/examples

.. |doi| image:: https://zenodo.org/badge/348020332.svg
    :alt: Zenodo
    :target: https://zenodo.org/badge/latestdoi/348020332
