.. _model_config_update_EM:

Updating a model
----------------
This plugin allows to update any components from a D-Emissions model. To do so, list the methods to update in a configuration file.
If needed edit your data library with new data sources required for the update and use the command:

.. code-block:: console

    activate hydromt-delwaq
    hydromt update demission path/to/model_to_update -o path/to/updated_model -i demission_update.yml -d data_sources.yml -vv
