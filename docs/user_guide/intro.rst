.. _intro_user_guide:

User Guide
==========

HydroMT can support the preparation of water quality models for `DELWAQ <https://oss.deltares.nl/web/delft3d/delwaq1/-/message_boards/category/205375>`_. 
As a reminder, a full hydrology - water quality simulation is in general separated in three parts:

- the hydrological model, for example wflow_sbm, that predicts water movements through the catchments.
- an emission model, D-Emission (also referred to as EM in HydroMT and this documentation), that predicts 
  the quantity (fluxes) of substances being released from sources of pollutants to the surface waters (mass/time).
- a fate and transport model, D-Water Quality (also referred to as WQ in HydroMT and this documentation), 
  that predicts the fate and transport of substances in the surface waters (concentrations in mass / volume).

For an explanation about how to use the HydroMT delwaq plugin for setting up a model, follow one of the links below:

.. toctree::
   :maxdepth: 1

   build_config_main.rst
   generic_delwaq_EM.rst
   generic_delwaq_WQ.rst
   coupling_main.rst