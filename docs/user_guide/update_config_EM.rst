.. _model_config_update_EM:

Updating a model
----------------
This plugin allows to update any components from a D-Emissions model. To do so, list the methods to update in a configuration file.
If needed edit your data library with new data sources required for the update and use the command:

.. code-block:: console

    activate hydromt-delwaq
    hydromt update delwaq path/to/model_to_update -o path/to/updated_model -i delwaq_update.ini -d data_sources.yml -vv

.. toctree::
   :hidden:
   :maxdepth: 2
   :caption: Available tutorials:

   Example: Add (global) emission data <../_examples/adding_global_emission.ipynb>
   Example: Add (local) emission data  <../_examples/adding_local_emission.ipynb>
   Example: Update forcing data <../_examples/update_model_forcing.ipynb>
