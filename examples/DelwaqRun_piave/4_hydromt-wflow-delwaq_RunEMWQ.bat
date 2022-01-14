REM Manual adjust the Delwaq .inp files to include the required data for the substance to model
copy hydroMT-setup\for_EM\espace.inp EM_piave
copy hydroMT-setup\for_WQ\delwaq.inp WQ_piave
copy hydroMT-setup\for_WQ\tra_em.def WQ_piave
pause

REM Run EM
cd EM_piave
REM Read input
..\Delwaq_bin\Delwaq\bin\delwaq1.exe espace.inp -A -p..\Delwaq_bin\EM-Plugin\Tables\proc_def 
REM Simulate
..\Delwaq_bin\Delwaq\bin\delwaq2.exe espace.inp -openpb..\Delwaq_bin\EM-Plugin\x64\Release\D3Dwaq_OpenPL.dll
erase *.wrk
cd ..
pause

REM Run WQ
cd WQ_piave
REM Read input
..\Delwaq_bin\Delwaq\bin\delwaq1.exe delwaq.inp -A -p..\Delwaq_bin\WQ-Plugin\Tables-Integral\proc_def
REM Simulate
..\Delwaq_bin\Delwaq\bin\delwaq2.exe delwaq.inp -openpb..\Delwaq_bin\WQ-Plugin\x64\Release\D3Dwaq_OpenPL.dll
erase *.wrk
cd ..
pause