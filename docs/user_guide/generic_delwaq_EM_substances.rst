.. _generic_delwaq_EM_substances:

Selection of substance(s) and sources
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The very first part of preparing a water quality model is to define which substances and which emission sources you wish to include in the model. 

When a substance is chosen, the sources (or losses) the user wants to consider can be defined (e.g.: atmospheric deposition, domestic, agriculture...). 
The sources and their emission releases are then quantified or initialized based on available data. There are two types of sources in D-Emissions:

-  Type A: direct substance release rates by the source (usually in g/day, e.g. atmospheric deposition)
-  Type B: substance release rates based on a locator: for example a domestic source can be charaterized by the locator population (in cap) and its emission factor (in g/cap/day).

Then in D-Emissions, the substances can enter the catchment via different receptors that are:

-  Sewage system
-  Paved areas
-  Unpaved areas
-  Stormwater system
-  Surface waters directly
-  Soil

When substances reach the receptors, they can then undergo several processes such as burial or removal by treatment in a WWTP before reaching the final receptors of D-Emissions, 
the surface waters and the soil. The full processes included the D-Emissions model are represented in the scheme below:

.. image:: ../img/d-emission-pathways.png

As modelling of several substances with multiple sources and processes can be complex, the D-Emissions model has been designed to model
each substance one-by-one. This implies that in order to model several substances, the D-Emissions model needs to be run several times.