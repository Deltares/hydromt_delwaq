.. _intro_user_guide:

User Guide
==========

With the **HydroMT-delwaq plugin**, users can easily benefit from the rich set of tools of the 
`HydroMT package <https://deltares.github.io/hydromt/latest/index.html>`_ to build and update
water quality models from available global and local data.
As a reminder, a full hydrology/hydrodynamic - water quality simulation is in general separated in three parts:

- the hydrological/hydrodynamic model, for example Wflow, that predicts water movements through the catchments.
- an emission model, D-Emissions (also referred to as EM in this documentation), that predicts 
  the quantity (fluxes) of substances being released from sources of pollutants to the surface waters (mass/time).
- a fate and transport model, D-Water Quality (also referred to as WQ in this documentation), 
  that predicts the fate and transport of substances in the surface waters (concentrations in mass / volume).

This plugin assist the delwaq modeller in:

- Quickly setting up a base D-Emissions or D-Water Quality model on top of an hydrologic Wflow model
- Making maximum use of the best available global or local data
- Adjusting and updating components of a delwaq model and their associated parameters in a consistent way

The prepared delwaq models with HydroMT are almost ready-to-run, and the users only needs to adjust the `.inp` configuration
file to finalise the setup of his water quality model run and substances to model.

.. toctree::
   :maxdepth: 2
   :hidden:

   generic_delwaq_EM.rst
   generic_delwaq_WQ.rst
   process_analyze.rst
   coupling_main.rst