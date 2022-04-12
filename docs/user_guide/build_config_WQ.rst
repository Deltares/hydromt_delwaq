.. _build_config_WQ:

Preparing the WQ model with HydroMT
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
D-Water Quality (Delwaq) is an unstructured model and composed of unordered segments. D-Water Quality segments represent one environmental layer unit (eg open water, soil..). 
The direction of the flows in Delwaq between each segment is defined in the pointer file. This pointer file defines not only the direction 
of lateral flows between segments (eg discharge), but also the direction of vertical flows (flows between the different layers/compartments of a single geographic cell). 
This means that external flows into / out of a segment (for example precipitation from the atmosphere) are defined in the pointer as flows between a segment and 
a boundary (for precipitation the boundary is the atmosphere).

Compared to D-Emission, lateral transport of the substance in the surface water between Delwaq segments is considered and a pointer file is then necessary.

To build an WQ model with HydroMT you can use the ** :ref:`build command line <model_config_build>` ** to prepare all the different ** :ref:`components <model_components>` ** of your model.
Below is an example :download:`ini file <../_static/delwaq_build_WQ_TN.ini>` for the **Total Nitrogen (TN)** released by **households** example.

As a reminder, for a WQ model, the required hydrological data (from wflow_output in our example) is:

- surface runoff
- surface water volumes
- total water entering the surface waters (e.g. precipitation, exfiltration from the soil, evaporation…) in order to close the surface water mass balance

.. literalinclude:: ../_static/delwaq_build_WQ_TN.ini
   :language: Ini

.. note::

   In some cases, you may require some additional GIS data in order to run a WQ model, for example, when you do not need to run the D-Emission model beforehand, 
   but include point sources that emit pollutants into the surface water directly. When this is the case, you can still use the various [setup_emission_*] components of hydroMT_delwaq.