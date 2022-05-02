.. _generic_delwaq_WQ_prepare:

Editing the WQ run information
------------------------------

.. Tip::

    This page contains additional information on the generic D-Water Quality module, not direclty related to the HydroMT-delwaq plugin.

In order to run the D-Water Quality module of DELWAQ, an input file is used to set up all the relevant run information such as substances and processes to include, run time information, data and model parameters 
to use etc. This is set up in a **.inp** file. A template of such a :download:`file <../_static/delwaq.inp>`, for an WQ run, is shown below.

The input file is separated into 10 input blocks each requiring different type of information. You can find more info on what is expected 
where in the `Delwaq input file documentation <content.oss.deltares.nl/delft3d/manuals/D-Water_Quality_Input_File_Description.pdf>`_. In short, the ten blocks are:

- B1: Identification, selected substances
- B2: Timers, integration, monitoring
- B3: Grid and values of the volumes
- B4: Hydrodynamic data
- B5: Open boundary conditions
- B6: Loads and withdrawals
- B7: Process steering
- B8: Initial conditions
- B9: Model output
- B10: Statistical output

With hydroMT, different files were prepared containing information that can be directly linked and included in this input file. These are:

- The different files in the *config* folder: each filename has a prefix containing the block number it should be inserted in (ex: B1_timestamp.inc).
- Emission and other grid data stored in *staticdata* folder
- Time-dependant data such as hydrological fluxes stored in *dynamicdata* folder.
   
Depending on their format (binary or ASCII), the files can be included in the input file using the following syntax.

- For ASCII files: keyword **INCLUDE** + path/to/ascii

::

    INCLUDE 'config\B1_timestamp.inc'

- For binary files: keywords **ALL BINARY_FILE** + path/to/binary

::

    ALL BINARY_FILE 'dynamicdata\hydrology.bin'

In the template below, all (mandatory) files produced by HydroMT have been added and linked. Some explanations on the specificities of the settings for 
the D-Water Quality linked to one-substance D-Emissions models:

- In block 1, you can see that we are going to model TN as a tracer for simplification (less expert processes to define). The generic name for a tracer in Delwaq is **cTR1**:

::

    ;Nr      Name
    1       cTR1

- The boundary conditions and initial concentrations of our substance **cTR1** needs to be defined in block 5 and 8 respectively and can be initialized with zero:

::

    ITEM BD_1    CONCENTRATIONS cTR1 DATA 0.0
    #5

::

    INITIALS  cTR1 DEFAULTS  0.0
    #8

- The link of our tracer **TRA** with D-Emissions output and to the relevant process is activated in block 7 of the input file (1st line), the 2nd line is added to map all of the emitted tracer *Tra* to a conservative tracer substance *cTR1* (Note: not used in this tutorial, but EM also supports to simulate a decayable tracer *dTR1*, which requires to also define decay rates):

::

    CONSTANTS Active_EM_TRA DATA 1
    CONSTANTS TratocTR1     DATA 1.0

.. literalinclude:: ../_static/delwaq.inp