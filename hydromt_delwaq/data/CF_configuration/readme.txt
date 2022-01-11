#Description of the skeleton config setup

sibrenloos:
Notes from working session 14 July 2021:
• For Ribasim skeleton config besides ribasim there is an environment folder
	o If more models use environment folder this is potentially not multi model proof. To discuss: e.g. share same Python version / environment for all models.
• delwaq moduleDataset setup needs to be ‘refurbished’ to meet HydroMT, current version still based on old modelbuilder
• delwaq parameters not known on forehand:
	o option1: read substances from delwaq inc file using read_model function in HydroMT
	o option2: limit list of substances to ~20 (e.g. BOD, E-coli, TN, TP, cTR1, dTR1, and for fractions: Initial, Area1, Area2, Area3) based on template list and use CF to turn these on or off; make use of html page with dropdown list
	o FEWS has functionality that supports that if a locationset is not found (or is empty?) it is not shown in spatialDisplay; when using (location)attributes and use a separate grid per substance, with the attributes one can indicate if the grid exists or not (for a complete explanation please contact Peter Gijsbers)
• connecting hydroMT to CF:
	o Use Build command, where last option is export to CF (Yes / No option); if Yes then the CF-export transfers the model to moduleDataset.zip and perform other actions required for integration in FEWS
• Fews ids such as zip file names may not exceed limit of 64 characters

gijsber:
A few more decisions from July 14:
• For generic configuration items (applicable for all models of a certain modelcode):
	o filenames should include the modelcode-reference (eg wflow, delwaq or demm)
	o files are stored in a subfolder for the model code (e.g. ModuleConfigFiles/delwaq/run_delwaq.xml)
• For modelschematisation specific files
	o all location/moduleinstance identifiers should include a modelcode reference to prevent clashes with other modelcodes
	o it is recommended to also include a modelschematisation reference in those identifiers
	o parameters that may occur in multiple codes (e.g. Q.sim) should get a modelcode reference postfix (e.g. Q.sim.wflow, Q.sim. rbs)
	o For easy management purposes it is recommendeded that each modelschematisation has its own files
	o filenames should include a modelcode-references. and are recommended to include a model schematization indentifier (e.g. LocationSets__delwaq_Rhine1x1km.xml)
	o Since we are modelling a region or basin, we introduce a rootfolder level for the region, followed by a modelcode folder and optionally followed by a model schematzation folder, i.e. RegionConfigFiles/Rhine/delwaq/Rhine1x1/grids_rhine1x1km.xml
	o Tree-structures in the GUI (Sptial Display, DisplayGroups, Filters) als are organized as // e.g,. Rhine/Water quality/BOD
	o The Task Tree will recieve its own design (work in progress)
	
Actions 15 oct 2021:
- create Delwaq templates based on structure above and c:\Users\loos_sb\OneDrive - Stichting Deltares\loos_sb\PROJECTEN\11205287-012_SO IMM Global Emissions\FEWS_Export\RIBASIM-CF\LATEST\Ribasim8.zip\ example (https://repos.deltares.nl/repos/ModelConverterRibasim/external/RIBASIM-CF/generator)
	- start with simple FEWS acc application as start (incl WFLOW)
	- create template config files
		- CF_configuration\Ribasim8\Modules\ribasim\system\cf_configurator\config_templates\
			
			- \cf_configurator\config_templates\AdapterConfigFiles\ribasim_postadapter.xml
			- \cf_configurator\config_templates\AdapterConfigFiles\ribasim_preadapter.xml
			
			- \cf_configurator\config_templates\RegionConfigFiles\Topology.xml
			- \cf_configurator\config_templates\RegionConfigFiles\ModuleInstanceSets.xml
			- \cf_configurator\config_templates\DisplayConfigFiles\SpatialDisplay.xml
			
			- \cf_configurator\config_templates\SystemConfigFiles\LocationIcons.xml
			
			- \cf_configurator\definitions\parameters\all_rib_parameters_complete.csv
			- \cf_configurator\definitions\mappings\spatialdisplay_mapping.csv
			
			- \cf_configurator\config_templates\RegionConfigFiles\ModifierTypes.xml
			- \cf_configurator\config_templates\DisplayConfigFiles\ModifiersDisplay.xml
			
		- use Tracer calculation Moselle from FEWS-Acc (c:\Users\loos_sb\REPOS\wflow\demo_peru\fews_Moselle_WQimport\) 
			- adapt GA to run EM and WQ
	- create html file that writes ini file
		- CF_configuration\Ribasim8\Modules\bin\startup\index_ribasim.html
		- CF_configuration\Ribasim8\Preparation\newbasin\ribasim_cf_configurator.ini
	- use and adapt py scripts for delwaq => delwaq_generator
		- see c:\Users\loos_sb\REPOS\hydromt_delwaq\hydromt_delwaq\data\CF_configuration\Ribasim8\Modules\bin\pyscripts\ribasim_generator\
	- test python script to generate complete config
	- integrate in hydroMT (requires wflow_generator besides delwaq generator)
	
	
[15-10-2021 10:31] Sibren Loos
Hoi Peter,
Ik wil met de FEWS export vanuit HydroMT aan de slag (richting CF dus)., testcase Delwaq. 
Zijn er nog recente ontwikkelingen binnen CF waar ik rekening mee moet houden?
[15-10-2021 10:33] Peter Gijsbers
Kijk nog even naar de laatste Ribasim versie die ik bij Mark in de brievenbus heb gezet
[15-10-2021 10:48] Sibren Loos
Ribasim8.zip\Modules\bin\pyscripts\ribasim_cf_configurator.py is denk ik het hoofd py script om de Ribasim config op te zetten obv templates? Dus deze moeten we in HydroMT zien te integreren (of ben je daarover al met Mark in gesprek?) en ook beschikbaar maken voor andere model?Moet ik verder nog iets weten? Zo niet, kijk ik of ik Delwaq op vergelijkbare manier op kan zetten en daarna de afstemming/integratie met hydroMT
[15-10-2021 10:49] Peter Gijsbers
correct, staan ook op svn https://repos.deltares.nl/repos/ModelConverterRibasim/external/RIBASIM-CF/generator
[15-10-2021 10:50] Sibren Loos
is er al iets met Mark afgestemd wbt integratie hydroMT?
[15-10-2021 10:56] Peter Gijsbers
inpassing in hydroMT is nog niet voldoende afgestemd. wel willen we daar volgend jaar mee aan de slag onder de noemer integratie RIBASIM-WFlow
[15-10-2021 12:13] Sibren Loos
klopt het dat jij start met een 'werkende' RIBASIM in Ribasim8\Preparation\newbasin\ (die bijv. mbv HydroMT gemaakt kan worden)? En integreer (append_config) je het daarna in FEWS obv ini file (aangemaakt mbv html pagina) en template config files?
[15-10-2021 12:14] Peter Gijsbers
correct
[15-10-2021 12:15] Sibren Loos
Mijn FEWS startapplic is dan aangemaakt mbv FEWS accelarator, waarna ik mbv delwaq_generator het hydroMT delwaq model op pak (en tzt ook wflow_generator). Klinkt logisch/werkbaar
de Update Ribasim Model topology node zal ik dan niet hebben (of wordt een link naar HydroMT tzt) in de eerste delwaq versie
[15-10-2021 12:17] Peter Gijsbers
ja, maar ik denk dat we volgend jaar een omhullend iets krijgen voor HydroMT zodat het start proces met toelevering van data (en dus de prepare Ribasim Basin workflow) buiten Fews wordt gedraaid
die update Ribasim Model is inderdaad alleen nodig om het Ribasim netwerk verder af te ronden
[15-10-2021 12:20] Sibren Loos
ok dat lijkt me handig buiten FEWS, maar de stap van een startConfig met appendConfig naar nieuwe Config blijft. In hydroMT geef je dan als argument je startConfig mee evt. optie createConfig mbv FEWS accelarator onderdelen.Ik kan voor nu in ieder geval even aan de slag.Het is mij helder nu, de opzet/logica
[15-10-2021 12:22] Peter Gijsbers
inderdaad je hebt een skeleton (startConfig) en het appendConfig process (op basis van je schematisatie) blijft in beide situaties bestaan. als je een omhullende buiten Fews hebt zal die appendConfig buiten Fews gedraaid worden en in combinatie met de startConfig een 'eind-applicatie'leveren


[15-10-2021 10:49] Sibren Loos
Wat doet Modules\bin\venv_ribasim precies, zet je een specifieke python environment voor ribasim op?
[15-10-2021 10:49] Peter Gijsbers
correct
[15-10-2021 10:50] Sibren Loos
waar is dat voor nodig? hebben alle modellen een eigen env nodig straks?
[15-10-2021 10:51] Peter Gijsbers
nee, dit had ook als complete python installatie meegeleverd inclusief de sitepackages 
[15-10-2021 10:52] Peter Gijsbers
die venv is een overblijfsel uit het verleden