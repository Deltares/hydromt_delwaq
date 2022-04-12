.. _generic_delwaq_EM:

Working with the D-Emissions model
----------------------------------
The first step in solving or mitigating many water quality and health problems, is to get a quantitative overview of the emissions of substances causing such problems.
`D-Emissions <www.deltares.nl/en/software/module/D-Emissions/>`_ provides a perfect tool to assess the quantity and spatial distribution of the emitted substances in the catchment. 
Because it simulates the true cause-effect chain, it provides predictive power and allows for “what if” scenarios driven by policy projections and implemented measures.


D-Emissions is a library of the Delwaq software which, in combination with rainfall-runoff and hydraulic modules (e.g. wflow), may form an integral catchment model.

In this section, we will detail the full steps to prepare and run D-Emissions models:

-  :ref:`Selection substance(s) <generic_delwaq_EM_substances>` of interest and their relevant sources (based on available data)
-  :ref:`Preparing the EM model with HydroMT <build_config_EM>` (EM model setup and emission data preparation) 
-  :ref:`Editing the EM run information <generic_delwaq_EM_prepare>`
-  :ref:`Running D-Emissions <generic_delwaq_EM_run>` (EM-Plugin)

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Table of Contents
   
   generic_delwaq_EM_substances.rst
   generic_delwaq_EM_prepare.rst
   generic_delwaq_EM_run.rst