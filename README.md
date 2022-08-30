#LBRM
contact: james.kessler@noaa.gov

This is a stand-alone version of the large basin runoff model which is a component of the Great Lakes Seasonal Hydrologic Forecast System (GLSHFS).  
This version has 3 different options available to represent Potential Evapotranspiration.

For a more comprehensive description see the README at github.com/NOAA-GLERL/GLSHFS



source:
	 contains all the FORTRAN code for LBRM

parmfiles:
	contains Parameter files for each subbasin. (As of Aug 2022, these files do not have 
	Net Radiation climatology so, these values must be generated PET must be < 3, see PETMETHOD below)

run:
	contains sample input files required for running the standalone lbrm:
	- bounds_sup01.txt   
		boundary/initial conditions
	- config.in          
		configuration file: points to each of these files, and sets the Potential Evapotranspiration scheme
	- met_sup01.txt 
		meterological forcing for run 
		this version contains a column for Net Radiation but this can be omitted if not using PETMETHOD=3
	- param_sup01.txt    
		parameters for run, this version contains Air T climatology and Net Rad Climatology
		this version contains a block for both Air T and a block for Net Rad
		(Air T can be omitted if PETMETHOD=1; Net Rad can be omitted if PETMETHOD < 3)
		The values in this file may or *MAY NOT BE REALISTIC* and are just used as an exmaple



Gettting Started:
1. compile lbrm in the source dir and cp/link executable to "run" dir
2. modify config.in as desired
3. ./lbrm config.in





Changes related to Potential Evapotranspiration:
	PETMETHOD is now a run-time setting in the config file
		PETMETHOD=1;  the original PET scheme that LBRM used upon its creation
		PETMETHOD=2;  2016 method that applies Clausius-Clayperon relationship
		PETMETHOD=3;  2021 method based on Priestly-Taylor scheme 

	PETMETHOD must be set to 1,2, or 3, otherwise desecrptive error occurs
	If PETMETHOD=2, Air T must exist in the parameter file (otherwise error occurs)
	If PETMETHOD=3, Above line is true + Net Rad must exist in the parameter file and meteorological file (otherwise error occurs)
	However, "extra" data is handled gracefully (e.g. a Net Rad column in the metfile is simply ignored if PETMETHOD < 3)

	additionally, the PET method being used is written to stdout AND to the output file (line 2)
	

	
	



