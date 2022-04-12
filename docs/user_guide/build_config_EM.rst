.. _build_config_EM:

Preparing the EM model with HydroMT
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
D-Emissions (Delwaq) is an unstructured model and composed of unordered segments. D-Emissions segments represent one environmental layer unit (eg open water, soil..). 
The direction of the flows in D-Emissions between each segment is defined in the pointer file. This pointer file defines not only the direction 
of lateral flows between segments (eg discharge), but also the direction of vertical flows (flows between the different layers/compartments of a single geographic cell). 
This means that external flows into / out of a segment (for example precipitation from the atmosphere) are defined in the pointer as flows between a segment and 
a boundary (for precipitation the boundary is the atmosphere).

By only modelling substances one-by-one, the generic D-Emissions version only considers one environment layer which is the **emission or substance layer**. As D-Emissions does not deal 
with fate and transport (no lateral flux between segments), there is then no need to define either a pointer or a boundary file. In essence, D-Emissions calculates the various emitted loads per wflow cell.

To build an EM model with HydroMT you can use the ** :ref:`build command line <model_config_build>` ** to prepare all the different ** :ref:`components <model_components>` ** of your model.
Below is an example :download:`.ini file <../_static/delwaq_build_EM_TN.ini>` for the **Total Nitrogen (TN)** released by **households** example.

As a reminder, for an EM model, the required hydrological data (from wflow_output in our example) is:

- precipitation
- the amount of the precipitation that infiltrates into the soil from unpaved areas
- the amount of the precipation that goes directly to surface runoff from paved areas
- the amount of the precipation that goes directly to surface runoff from unpaved areas

.. literalinclude:: ../_static/delwaq_build_EM_TN.ini
   :language: Ini

.. note::

   In this version of EM, the *the amount of the precipitation that **infiltrates** into the soil from unpaved areas* is dealt with as a sink, since the WQ version does not yet contain a Soil compartment.