.. _build_config_EM:

Building a model
================

This plugin allows to build a complete D-Emissions model from available data. Once the configuration and 
data libraries are set, you can build a model by using:

.. code-block:: console

    activate hydromt-delwaq
    hydromt build delwaq path/to/built_model "{'wflow': path/to/wflow_model}" --opt global.mtype=EM -i delwaq_build.ini -d data_sources.yml -vv

The recommended `region options <https://deltares.github.io/hydromt/latest/user_guide/cli.html#region-options>`_ 
for a proper implementation of this model are:

- model

Alternatively, to start from a complete new region, you can start by first using HydroMT to build the linked hydrological/hydraulic 
model for your case and then build delwaq.

.. warning::

  As of now, DELWAQ models can only be built on top of existing Wflow models.

.. _config_file_EM:

Configuration file
------------------
Settings to build or update a D-Emissions model are managed in a configuration file. In this file,
every option from each :ref:`model methods <em_methods>` can be changed by the user
in its corresponding section.

In addition to the model components, the ini also contains a [global] section 
in order to specify if the model to build/update is either an EM (D-Emissions) or a WQ (D-Water Quality) model ('mtype' option).

.. code-block:: console

    [global]
    mtype = EM                     # type of Delwaq model ['EM', 'WQ']

Below is an example of :download:`.ini file <../_static/delwaq_build_EM_TN.ini>` that can be used to build a **EM model** 
for modelling **Total Nitrogen (TN)** released by **households**.

.. literalinclude:: ../_static/delwaq_build_EM_TN.ini
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

  .. toctree::
     :hidden:

     Example: Build D-Emission model <../_examples/build_EM_model.ipynb>