.. _generic_delwaq_EM_prepare:

Editing the EM run information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In order to run the EM-Plugin of DELWAQ, an input file is used to set up all the relevant run information such as substances and processes to include, run time information, data and model parameters 
to use etc. This is set up in a **.inp** file. A template of such a :download:`file <../_static/espace.inp>`, for an EM run, is shown below.

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
the one-substance D-Emissions model:

- In block 1, you can see that in the one-substance D-Emissions model, instead of the substances (nitrogen, TSS...), we model the emission receptors compartments in a specific order. 
  With this, emission modelling can be more generic, run time faster and there is no need to update the procdef tables for each D-Emissions model case.

::

    ;Nr     Name
    1       Sew 
    2       Pav 
    3       Unp 
    4       Stw 
    5       Sfw 
    6       Soi 

- The different emission sources are defined in block 7. Currently, in the one-substance D-Emissions model, **one type A source** and **three type B sources** can be defined. All names 
  respect a certain convention (EV_A01, EV_B01). The user should also define to which receptor (from block 1) the source emission goes to. Here is an example for a type B source :

::

    ; source Type B01
    PARAMETERS EV_B01        ALL BINARY_FILE 'staticdata\GHS-POP_2015.dat' ; locator/EV variable type B    (cap)                
    CONSTANTS  EF_B01           DATA 2.2          ; emission factor                                        (Kg/d/cap)           
    CONSTANTS  B01toSew         DATA 1.0          ; released fraction receptor                             (-)                 
    CONSTANTS  B01toPav         DATA 0.0          ; released fraction receptor                             (-)                 
    CONSTANTS  B01toUnp         DATA 0.0          ; released fraction receptor                             (-)
    CONSTANTS  B01toSoi         DATA 0.0          ; released fraction receptor                             (-)
    CONSTANTS  B01toStw         DATA 0.0          ; released fraction receptor                             (-)
    CONSTANTS  B01toSfw         DATA 0.0          ; released fraction receptor                             (-)

- All unused type A and/or B sources should still be listed in the input file with zero data:

::

    ; source Types NOT USED
    CONSTANTS  EV_A01        DATA 0.    
    CONSTANTS  EV_B02        DATA 0.    
    CONSTANTS  EV_B03        DATA 0.    

.. note::

  The unit of emission factor EF is kg/d/X, where X is the emission unit provided in the Emission Variable EV, e.g. capita for Population.
  Here the EF factor was set to a constant dummy value but you can use HydroMT to prepare and distribute 'real' values from local/global data.

.. literalinclude:: ../_static/espace.inp