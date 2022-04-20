.. _model_config_build:

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

For an explanation about how to prepare the delwaq models, follow one of the links below:

.. toctree:: 
   :maxdepth: 2
   
   build_config_EM.rst
   build_config_WQ.rst   
   
.. toctree::
   :maxdepth: 2
   :caption: Available tutorials:
   
   Example: Build an emission model <../_examples/build_EM_model.ipynb>
   Example: Build a water quality model <../_examples/build_WQ_model.ipynb>
   
   
 