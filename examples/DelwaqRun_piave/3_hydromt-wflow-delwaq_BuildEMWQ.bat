REM Build Delwaq EM and WQ case from wflow
copy hydroMT-setup\for_EM\delwaq_build_EM.ini .\
copy hydroMT-setup\for_EM\delwaq_update_EM_forcing.ini .\
copy hydroMT-setup\for_WQ\delwaq_build_WQ.ini .\
copy hydroMT-setup\for_WQ\local_sources.yml .\
pause
call activate hydromt-delwaq
hydromt build delwaq EM_piave "{'wflow': './wflow_piave'}" -i ./delwaq_build_EM.ini -vvv 
pause
hydromt update delwaq EM_piave -i ./delwaq_update_EM_forcing.ini -d ./local_sources.yml -vvv 
pause
hydromt build delwaq WQ_piave "{'wflow': './wflow_piave'}" -i ./delwaq_build_WQ.ini -vvv
pause