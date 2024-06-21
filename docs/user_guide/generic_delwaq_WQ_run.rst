.. _generic_delwaq_WQ_run:

Running D-Water Quality
-----------------------

.. Tip::

    This page contains additional information on the generic D-Water Quality module, not direclty related to the HydroMT-delwaq plugin.

In order to run DELWAQ for D-Water Quality, we need to include the WQ-Plugin to the run. DELWAQ is also run via command lines using two executables, one that reads and process all the information and data
in the input file **delwaq1.exe**, and one that runs DELWAQ with this data **delwaq2.exe**. The command line also include the path to the input file (**delwaq.inp**) and to specific libraries in the WQ-Plugin (or D-Water Quality plugin). You can also notice that there is no direct link to the outdata_em files either in the delwaq.inp or in the run-WQ-tra.bat files. This is why respecting the naming conventions of either xxx_em.bin and xxx_em.txt or xxx_em.def are mandatory.

::

    erase *.wrk

    REM Read input
    ..\Delwaq\bin\delwaq1.exe delwaq.inp -A -p..\WQ-Plugin\Tables-Integral\proc_def

    REM Simulate
    ..\Delwaq\bin\delwaq2.exe delwaq.inp -openpb..\WQ-Plugin\x64\Release\D3Dwaq_OpenPL.dll
    pause
