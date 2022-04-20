.. _model_config_data:

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
   :caption: Available tutorials:
  
   Example: Add (global) emission data <../_examples/adding_global_emission.ipynb>
   Example: Add (local) emission data  <../_examples/adding_local_emission.ipynb>
   Example: Plot DELWAQ data <../_examples/plot_delwaq.ipynb>