.. _model_config:

Model configuration
===================

This HydroMT plugin provides an implementation for DELWAQ in order to build, update or clip from 
command line. Specific details on the HydroMT CLI methods can be found in the
`core documentation <https://deltares.github.io/hydromt_plugin/latest/user_guide/cli.html>`_ 

Introduction
------------
HydroMT can support the preparation of water quality models for `DELWAQ <https://oss.deltares.nl/web/delft3d/delwaq1/-/message_boards/category/205375>`_. 
As a reminder, a full hydrology - water quality simulation is in general separated in three parts:

- the hydrological model, for example wflow_sbm, that predicts water movements through the catchments.
- an emission model, D-Emission (also referred to as EM in HydroMT and this documentation), that predicts 
  the quantity (fluxes) of substances being released from sources of pollutants to the surface waters (mass/time).
- a fate and transport model, D-Water Quality (also referred to as WQ in HydroMT and this documentation), 
  that predicts the fate and transport of substances in the surface waters (concentrations in mass / volume).

The Delwaq model class in HydroMT is able to prepare and interact both with EM or WQ types of Delwaq models.

Configuration file
------------------
Settings to build or update a Delwaq model are managed in a configuration file. In this file, 
every option from each :ref:`model component <model_components>` can be changed by the user 
in its corresponding section.

In addition to the model components, the ini also contains a [global] section 
in order to specify if the model to build/update is either an EM or a WQ model ('mtype' option).

.. code-block:: console

    [global]
    mtype = WQ                     # type of Delwaq model ['EM', 'WQ']

Most of the options are the same for preparing either a EM or WQ type of model. The main difference is with the 
hydrological information needed by each model type. For an **EM** type of model, these are:

- **precip**: precipitation
- **infilt**: the amount of the precipitation that infiltrates into the soil from unpaved areas
- **runPav**: the amount of the precipation that goes directly to surface runoff from paved areas
- **runUnp**: the amount of the precipation that goes directly to surface runoff from unpaved areas

For a **WQ** type of model, these are:

- **run**: surface runoff
- **vol**: surface water volumes
- **inwater** - **inwaterInternal**: total water entering the surface waters (e.g. precipitation, exfiltration from the soil, 
  evaporation…) in order to close the surface water mass balance.

Below is an example of ini file that can be used to build a complete **EM model**
:download:`.ini file <../examples/delwaq_build_EM.ini>`. Each section corresponds 
to a model component with the same name.

.. literalinclude:: ../examples/delwaq_build_EM.ini
   :language: Ini

Below is an example of ini file that can be used to build a complete **WQ model**
:download:`.ini file <../examples/delwaq_build_WQ.ini>`. Each section corresponds 
to a model component with the same name.

.. literalinclude:: ../examples/delwaq_build_WQ.ini
   :language: Ini

Selecting data
--------------
Data sources in HydroMT are provided in one of several yaml libraries. These libraries contain required 
information on the different data sources so that HydroMT can process them for the different models. There 
are three ways for the user to select which data libraries to use:

- If no yaml file is selected, HydroMT will use the data stored in the 
  `hydromt-artifacts <https://github.com/DirkEilander/hydromt-artifacts>`_ 
  which contains an extract of global data for a small region around the Piave river in Northern Italy.
- Another options for Deltares users is to select the deltares-data library (requires access to the Deltares 
  P-drive). In the command lines examples below, this is done by adding either **-dd** or **--deltares-data** 
  to the build / update command line.
- Finally, the user can prepare its own yaml libary (or libraries) (see 
  `HydroMT documentation <https://deltares.github.io/hydromt/latest/user_guide/data.html>`_ to check the guidelines). 
  These user libraries can be added either in the command line using the **-d** option and path/to/yaml or in the **ini file** 
  with the **data_libs** option in the [global] sections.

Building a model
----------------
This plugin allows to build a complete DELWAQ model from available data. Once the configuration and 
data libraries are set, you can build a model by using:

.. code-block:: console

    activate hydromt-delwaq
    hydromt build delwaq path/to/built_model "{'wflow': path/to/wflow_model}" -i delwaq_build.ini -d data_sources.yml -vv

The recommended `region options <https://deltares.github.io/hydromt/latest/user_guide/cli.html#region-options>`_ 
for a proper implementation of this model are:

- model

Alternatively, do start from a complete new region, you can start by first using HydroMT to build the linked hydrological/hydrualic 
model for your case and then build delwaq.

.. warning::

  As of now, DELWAQ models can only be built on top of existing wflow models.

Updating a model
----------------
This plugin allows to update any components from a DELWAQ model. To do so, list the components to update in a configuration file,
if needed edit your data library with new data sources required for the update and use the command:

.. code-block:: console

    activate hydromt-delwaq
    hydromt update delwaq path/to/model_to_update -o path/to/updated_model -i delwaq_update.ini -d data_sources.yml -vv