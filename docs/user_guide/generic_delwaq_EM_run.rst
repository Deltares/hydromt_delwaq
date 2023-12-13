.. _generic_delwaq_EM_run:

Running D-Emissions
-------------------

.. Tip::

    This page contains additional information on the generic D-Emissions module, not direclty related to the HydroMT-delwaq plugin.

In order to run DELWAQ for D-Emissions, we need to include the EM-Plugin to the run. DELWAQ is also run via command lines using two executables, one that reads and process all the information and data 
in the input file **delwaq1.exe**, and one that runs DELWAQ with this data **delwaq2.exe**. The command line also include the path to the input file (**espace.inp**) and to 
specific libraries in the EM-Plugin (or D-Emissions plugin). Here is an example:

::

    erase *.wrk
    
    REM Read input
    ..\Delwaq\bin\delwaq1.exe espace.inp -A -p..\EM-Plugin\Tables\proc_def 
    
    REM Simulate
    ..\Delwaq\bin\delwaq2.exe espace.inp -openpb..\EM-Plugin\x64\Release\D3Dwaq_OpenPL.dll
    pause