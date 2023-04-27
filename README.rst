HydroMT-delwaq: DELWAQ plugin for HydroMT
#########################################

|pypi| |docs_latest| |docs_stable| |codecov| |license| |binder|

What is HydroMT-delwaq?
-----------------------
HydroMT-delwaq is a plugin of the `HydroMT core <https://deltares.github.io/hydromt/latest/index.html>`_, a python package, developed by Deltares, to build 
and analyse environmental models. This plugin provides an implementation for the `DELWAQ <https://www.deltares.nl/en/software/module/d-water-quality/>`_ water quality engine. 
It details the different steps and explains how to use HydroMT to easily get started and work on your own DELWAQ model. WIth this plugin 
you can interact with both classic **D-Water Quality** models as well as **D-Emission** models.

With the **HydroMT-delwaq plugin**, users can easily benefit from the rich set of tools of the 
`HydroMT package <https://github.com/Deltares/hydromt>`_ to build and update 
DELWAQ models from available global and local data.

This plugin assists the modeller in:

- Quickly setting up a model and default parameter values
- Making maximum use of the best available global or local data
- Adjusting and updating components of a DELWAQ model and their associated parameters in a consistent way

Water quality and emission, fate and transport of pollutants through the landscape, rivers and oceans is very much linked to the water movements. 
For this reason, HydroMT strongly links DELWAQ modelbuilding to an underlying hydrology / hydraulic or hydrodynamic model. As of now, only the link 
between `Wflow and DELWAQ <coupling_wflow>`_ is supported.

This plugin supports both:

- emission models, `D-Emission <https://www.deltares.nl/en/software/module/D-Emissions/>`_ (also referred to as EM in hydroMT), that predicts the quantity (fluxes) of substances being released from sources of pollutants to the surface waters (mass/time).
- fate and transport models, `D-Water Quality <https://www.deltares.nl/en/software/module/d-water-quality/>`_ (also referred to as WQ in hydroMT), that predicts the fate and transport of substances in the surface waters (concentrations in mass / volume).


Why HydroMT-delwaq?
-------------------
Setting up spatial geoscientific models typically requires many (manual) steps 
to process input data and might therefore be time consuming and hard to reproduce. 
Like the other HydroMT model plugins, this package aims to make the model building process **fast**, **modular** and **reproducible** 
by configuring the model building process from a single *ini* configuration file
and **model- and data-agnostic** through a common model and data interface. 


How to use HydroMT-delwaq?
--------------------------
HydroMT-delwaq can be used as a **command line** application, which provides commands to *build*,
and *update* DELWAQ model with a single line, in combination with the HydroMT core's rich data processing capabilities.

How to cite?
------------
For publications, please cite our work using the DOI provided in the Zenodo badge |doi| that points to the latest release.


How to contribute?
------------------
If you find any issues in the code or documentation feel free to leave an issue on the `github issue tracker. <https://github.com/Deltares/hydromt_delwaq/issues>`_
You can find information about contributing to HydroMT at our `Contributing page <https://deltares.github.io/hydromt/latest/dev/contributing>`_.


.. |pypi| image:: https://badge.fury.io/py/hydromt_delwaq.svg
    :target: https://pypi.org/project/hydromt_delwaq/
    :alt: Latest PyPI version

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