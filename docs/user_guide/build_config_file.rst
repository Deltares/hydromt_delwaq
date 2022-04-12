.. _model_config_file:

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
:download:`.ini file <../examples/examples/delwaq_build_EM.ini>`. Each section corresponds 
to a model component with the same name.

.. literalinclude:: ../examples/examples/delwaq_build_EM.ini
   :language: Ini

Below is an example of ini file that can be used to build a complete **WQ model**
:download:`.ini file <../examples/examples/delwaq_build_WQ.ini>`. Each section corresponds 
to a model component with the same name.

.. literalinclude:: ../examples/examples/delwaq_build_WQ.ini
   :language: Ini