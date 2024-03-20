.. _generic_delwaq_WQ_substances:

Selection of substance(s) and processes
=======================================

.. Tip::

    This page contains additional information on the generic D-Water Quality module, not direclty related to the HydroMT-delwaq plugin.

When setting-up D-Emissions, we decided to model the substance Total Nitrogen (TN) emitted by households. As TN is a  quite complex substance with many processes
and interactions with different forms of Nitrogen molecule, in the following tutorials, as a simplified example, we will model **Total Nitrogen (TN)** as a
simple constant tracer cTR.

For more information on available processes per substances, see the DELWAQ documentation:
 -  `<content.oss.deltares.nl/delft3d/manuals/D-Water_Quality_Processes_Technical_Reference_Manual.pdf>`_
 -  `<content.oss.deltares.nl/delft3d/manuals/D-Water_Quality_Processes_Library_Tables.pdf>`_
 -  `<content.oss.deltares.nl/delft3d/manuals/D-Water_Quality_Open_Proc_Lib_User_Manual.pdf>`_

D-Water Quality (DELWAQ) is an unstructured model and composed of unordered segments. D-Water Quality segments represent one environmental layer unit (eg open water, soil..).
The direction of the flows in DELWAQ between each segment is defined in the pointer file. This pointer file defines not only the direction
of lateral flows between segments (eg discharge), but also the direction of vertical flows (flows between the different layers/compartments of a single geographic cell).
This means that external flows into / out of a segment (for example precipitation from the atmosphere) are defined in the pointer as flows between a segment and
a boundary (for precipitation the boundary is the atmosphere).

Compared to D-Emissions, lateral transport of the substance in the surface water between DELWAQ segments is considered and a pointer file is then necessary.

For a D-Water Quality model, the required hydrological data is:

- surface runoff
- surface water volumes
- total water entering the surface waters (e.g. precipitation, exfiltration from the soil, evaporation…) in order to close the surface water mass balance

.. note::

   In some cases, you may require some additional GIS data in order to run a WQ model, for example, when you do not need to run the D-Emissions model beforehand,
   but include point sources that emit pollutants into the surface water directly. When this is the case, you can still use the various [setup_emission_*] components of HydroMT-delwaq.
