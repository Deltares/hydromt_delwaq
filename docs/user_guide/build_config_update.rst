.. _model_config_update:

Updating a model
----------------
This plugin allows to update any components from a DELWAQ model. To do so, list the components to update in a configuration file,
if needed edit your data library with new data sources required for the update and use the command:

.. code-block:: console

    activate hydromt-delwaq
    hydromt update delwaq path/to/model_to_update -o path/to/updated_model -i delwaq_update.ini -d data_sources.yml -vv
	
.. toctree::
   :maxdepth: 2
   :caption: Available tutorials:
   
   Example: Update monitoring locations  <../examples/examples/update_model_monitoring.ipynb>
   Example: Update forcing data <../examples/examples/update_model_forcing.ipynb>