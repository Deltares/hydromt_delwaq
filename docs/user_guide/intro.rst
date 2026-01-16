.. _intro_user_guide:

User Guide
==========

The user guide is organised through the following sections:

.. grid:: 3
    :gutter: 1

    .. grid-item-card::
        :text-align: center
        :link: 1_getting_started_hydromt/index
        :link-type: doc

        :octicon:`rocket;5em;sd-text-icon blue-icon`
        +++
        Getting started with HydroMT Core

    .. grid-item-card::
        :text-align: center
        :link: generic_delwaq_EM
        :link-type: doc

        :octicon:`file-directory-open-fill;5em;sd-text-icon blue-icon`
        +++
        Working with D-Emissions models

    .. grid-item-card::
        :text-align: center
        :link: generic_delwaq_WQ
        :link-type: doc

        :octicon:`file-directory-open-fill;5em;sd-text-icon blue-icon`
        +++
        Working with D-Water Quality models

    .. grid-item-card::
        :text-align: center
        :link: process_analyze
        :link-type: doc

        :octicon:`graph;5em;sd-text-icon blue-icon`
        +++
        Processing and Plots

    .. grid-item-card::
        :text-align: center
        :link: coupling_main
        :link-type: doc

        :octicon:`cpu;5em;sd-text-icon blue-icon`
        +++
        Technical Description

    .. grid-item-card::
        :text-align: center
        :link: https://deltares.github.io/hydromt/latest/user_guide/migration_guide/index.html

        :octicon:`file-moved;5em;sd-text-icon blue-icon`
        +++
        Migration Guide

With the **HydroMT-DELWAQ plugin**, users can easily benefit from the rich set of tools of the
`HydroMT package <https://deltares.github.io/hydromt/latest/index.html>`_ to build and update
water quality models from available global and local data.
As a reminder, a full hydrology/hydrodynamic - water quality simulation is in general separated in three parts:

- the hydrological/hydrodynamic model, for example Wflow, that predicts water movements through the catchments.
- an emission model, D-Emissions (also referred to as EM in this documentation), that predicts
  the quantity (fluxes) of substances being released from sources of pollutants to the surface waters (mass/time).
- a fate and transport model, D-Water Quality (also referred to as WQ in this documentation),
  that predicts the fate and transport of substances in the surface waters (concentrations in mass / volume).

This plugin assist the DELWAQ modeller in:

- Quickly setting up a base D-Emissions or D-Water Quality model on top of an hydrologic Wflow model
- Making maximum use of the best available global or local data
- Adjusting and updating components of a DELWAQ model and their associated parameters in a consistent way

The prepared DELWAQ models with HydroMT are almost ready-to-run, and the users only needs to adjust the `.inp` configuration
file to finalise the setup of his water quality model run and substances to model.

.. toctree::
   :caption: Table of Contents
   :maxdepth: 2
   :hidden:
   :titlesonly:

   1_getting_started_hydromt/index.rst
   generic_delwaq_EM.rst
   generic_delwaq_WQ.rst
   process_analyze.rst
   coupling_main.rst
