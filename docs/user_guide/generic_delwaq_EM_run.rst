.. _generic_delwaq_EM_run:

Running D-Emission
^^^^^^^^^^^^^^^^^^
In order to run delwaq for D-Emissions, we need to include the EM-Plugin to the run. Delwaq is also run via command lines using two executables, one that reads and process all the information and data 
in the input file **delwaq1.exe**, and one that runs Delwaq with this data **delwaq2.exe**. The command line also include the path to the input file (**espace.inp**) and to 
specific libraries in the EM-Plugin (or D-Emissions plugin). Here is an example:

::

    erase *.wrk
    
    REM Read input
    ..\Delwaq\bin\delwaq1.exe espace.inp -A -p..\EM-Plugin\Tables\proc_def 
    
    REM Simulate
    ..\Delwaq\bin\delwaq2.exe espace.inp -openpb..\EM-Plugin\x64\Release\D3Dwaq_OpenPL.dll
    pause