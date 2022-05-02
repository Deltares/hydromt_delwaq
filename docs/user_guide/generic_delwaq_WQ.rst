.. _generic_delwaq_WQ:

Working with the D-Water Quality model 
--------------------------------------
The D-Water Quality module simulates the far- and mid-field water and sediment quality due to a variety of transport and water quality processes. To accommodate 
these, it includes several advection diffusion solvers and an extensive library of standardized process formulations with the user-selected substances. Default 
processes allow to simulate for instance the decay of BOD and nitrification, elementary growth of algae and nutrient cycling, exchange of substances with the 
atmosphere, adsorption and desorption of contaminant substances and the deposition and resuspension of particles and adsorbed substances to and from the bed.


D-Water Quality is a module of the DELWAQ software which, in combination with rainfall-runoff and hydraulic modules (e.g. wflow), and the D-Emission module, 
may form an integral catchment model.


In this section, we will detail the steps to run a D-Water Quality model using HydroMT and inputs from D-Emission:

-  :ref:`Selection of substance(s) <generic_delwaq_WQ_substances>` of interest and their relevant processes in surface waters
-  :ref:`Preparing the WQ model with HydroMT <build_config_WQ>` (WQ model setup) 
-  :ref:`Editing the WQ run information <generic_delwaq_WQ_prepare>`
-  :ref:`Link with the outputs from D-Emissions <coupling_delwaq_single>`
-  :ref:`Running D-Water Quality <generic_delwaq_WQ_run>` (WQ-plugin)

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Table of Contents
   
   WQ_model.rst
   generic_delwaq_WQ_substances.rst
   build_config_WQ.rst
   update_config_WQ.rst
   generic_delwaq_WQ_prepare.rst
   generic_delwaq_WQ_run.rst