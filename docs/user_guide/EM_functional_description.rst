.. _EM_functional_description:

===============================
Detailed functional description
===============================

Generic Single Substance Emission Model (EM_GSS)
================================================

Introduction
------------

This process quantifies releases to various compartments and routes them to the final recipients: surface waters and soils. It is applied for a single substance.

Releases are defined for:

•	Atmospheric deposition: dry deposition as mass/time/surface area, wet deposition as mass/volume in precipitation;
•	A variable number of “type A” sources;
•	A variable number of “type B” sources.


The quantification of the releases of polutants associated to the various sources proceeds by the emission factor method. A variable collection of sources can be considered. Releases (L) of a pollutant “p” for a certain socio-economic activity “a” are calculated by multiplying an activity rate (ARa) by an emission factor for this activity and a certain pollutant (EFp,a):

Lp,a = ARa x EFp,a

Releases are distributed in space according to two alternative methods:

A.	The activity rate is known for a larger geographic area (region, country). The releases are first calculated at this aggregated level and then distributed in space using an auxiliary spatial variable called a “locator”.
B.	The activity rate is already a spatially distributed variable and can be directly used to calculate spatially variable releases.

The quantified releases are allocated to various initial receptors and routed towards the final receiving compartments: the soil system or the surface waters. Initial receptors can be: (1) the surface waters (directly), (2) the soil system (directly), (3) impermeable surfaces, (4) permeable surfaces, (5) the sewage collection system (combined or separated) and (6) the separate rainwater collection system.

The emissions related to stormwater and wastewater are simulated as follows:

•	The allocation of the wastewater releases from households to initial receptors should take into account the share of unconnected households and the share of households using septic tanks or other “individual or appropriate systems” (IAS, as defined under the Urban Wastewater Treatment Directive).
•	The substances washed off from impermeable areas find their way to a separate rainwater collection system, a combined collection system for stormwater and wastewater or to surface waters and soils in places where there is no collection system.
•	Separate collection systems discharge to surface waters, while a retention term can be defined that is partly allocated to soils (re-use or distribution of sludge).
•	Combined collection systems discharge via a WWTP, where the treatment level can be variable.
•	Combined collection systems feature overflow events (CSOs), during which the treatment is bypassed; CSOs are a fixed fraction or occur if a daily rainfall threshold is exceeded.

Implementation
--------------

The EM is currently configured for 6 compartments.

.. list-table::
   :widths: 5, 25
   :header-rows: 1

   * - Abbreviation
     - Description
   * - ``Sew`` 
     - Sewer system that receives wastewater and (optionally) stormwater
   * - ``Pav``
     - Paved or impermeable surfaces
   * - ``Unp``
     - Unpaved or permeable surfaces
   * - ``Stw``
     - Sewer system that receives only stormwater
   * - ``Sfw``
     - Final recipient surface water, all matter that ends up here is “emissions to surface waters”
   * - ``Soi``
     - Final recipient soil system, all matter that ends up here is “emissions to soils” 


In the EM software, these compartments are mathematically represented by substances. The mass of the simulated substances is expressed in grams. Transports from one compartment to another (in g s-1) are mathematically represented by “process fluxes” (transformations of one substance into another substance). This comes back in the mass balance output.

As the mathematical concept of a “substance” is used to represent compartments, EM has been set up to run for one “true” substance at a time.

Formulation
-----------


Input
-----

Output
------

Waste water management supportive process (GenWWman)
====================================================

Introduction
------------

Implementation
--------------

Formulation
-----------

Input
-----

Output
------

The supportive process produces the below output, that feeds directly into the main emission modelling process.

.. list-table::
   :widths: 5, 25, 5
   :header-rows: 1

   * - Name in model
     - Definition
     - Unit
   * - ``Eff_WWTP`` 
     - fraction of substance in WWTP influent that reaches the effluent (substance dependent)
     - (-)
   * - ``Sld_WWTP`` 
     - fraction of substance in WWTP influent that reaches the sludge (substance dependent)
     - (-)
   * - ``WWtoSew`` 
     - fraction of wastewater allocated to mixed sewers
     - (-)
   * - ``WWtoSfw`` 
     - fraction of wastewater allocated to surface waters
     - (-)
   * - ``WWtoSoi`` 
     - fraction of wastewater allocated to soils
     - (-)