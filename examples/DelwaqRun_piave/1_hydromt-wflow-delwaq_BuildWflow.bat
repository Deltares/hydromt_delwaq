REM Prepare wflow model
call activate hydromt-wflow
hydromt build wflow "./wflow_piave" "{'subbasin': [12.2051, 45.8331], 'strord': 4, 'bounds': [11.70, 45.35, 12.95, 46.70]}" -r 0.00833 -i hydroMT-setup/for_wflow/wflow_build.ini -v
call conda deactivate
pause