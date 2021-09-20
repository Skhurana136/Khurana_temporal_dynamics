# Spatio-temporal_heterogeneities_saturated_zone_microbial_dynamics
This work is associated with the following publication submitted for peer-review:
Khurana et al. 2021. Should we worry about surficial dynamics when assessing nutrient cycling in the groundwater?. Frontiers in Water.
Registration number: XX-XXXX-XX.

- This repository contains the following:
	- Python scripts used for processing simulation results, further analysis and generating graphs for the publication.
	- Files used for running simulations in OGSBRNS
	- Processed results (csv and xls files)

## Simulation Results
- Available on request with the authors.

## Subdirectory Simulations
- Trial324_DLL.mws: Describes the reaction network. This is processed in MAPLE. It generates a spread.m file.
- spread.m: This file generates fortran files that are used by BRNS to generate a DLL that is linked with OGS.
- brnsDLL.dll: The dynamically linked library to run OGSBRNS.
- ogs.exe: The OGS executable that is solving the groundwater flow component of the model and exchanging information with BRNS.
- Subdirectories:
	- EqualAR_1: Contains simulation files for 49 scenarios for the medium flow regime (as described in the manuscript) for Time Series T1.
	- EqualAR_2: Contains simulation files for 49 scenarios for the medium flow regime (as described in the manuscript) for Time Series T2.
	- EqualAR_5: Contains simulation files for 49 scenarios for the medium flow regime (as described in the manuscript) for Time Series T5.
	- SlowAR_1: Contains simulation files for 49 scenarios for the slow flow regime (as described in the manuscript) for Time Series T1.
	- SlowAR_2: Contains simulation files for 49 scenarios for the slow flow regime (as described in the manuscript) for Time Series T2.
	- SlowAR_5: Contains simulation files for 49 scenarios for the slow flow regime (as described in the manuscript) for Time Series T5.
	- FastAR_1: Contains simulation files for 49 scenarios for the fast flow regime (as described in the manuscript) for Time Series T1.
	- FastAR_2: Contains simulation files for 49 scenarios for the fast flow regime (as described in the manuscript) for Time Series T2.
	- FastAR_5: Contains simulation files for 49 scenarios for the fast flow regime (as described in the manuscript) for Time Series T5.
	
## Subdirectory Scripts
- All_figures.py: Generates figures used in the manuscript using *.csv files in the Results subdirectory.
- massfluxdata_generation_steadystate.py:
	- Processes the simulation results (numpy arrays) and outputs the mass flux at the inlet and outlet boundaries.
	- The results are saved for all the reactive species in a *.csv file (see Results subdirectory).
	- Access to functions in dsenvsci repository is required to run this script.
	- Output file: massflux_steadystate_BG.csv
- biomassdata_generate_steadystate.py:
	- Processes the simulation results (numpy arrays) and outputs the biomass in the domain.
	- The results are saved for all the reactive species in a *.csv file (see Results subdirectory).
	- Access to functions in dsenvsci repository is required to run this script.
	- Output file: biomass_steadystate_BG.csv
- massfluxdata_comparison_steadystate.py:
	- Compares the change in mass flux in all domains with that in the base case (homogeneous case) in each flow regime.
	- Input file required: massflux_steadystate_BG.csv.
	- Output file: massflux_comparison_steadystate_BG.csv
- biomass_comparison_steadystate.py:
	- Compares the biomass in all domains with that in the base case (homogeneous case) in each flow regime.
	- Input file required: biomass_steadystate_BG.csv.
	- Output file: biomass_comparison_steadystate_BG.csv
- Tracer_studies_saturated_flow.py:
	- Processes the tracer studies' simulation results (data set is large, and it can be made available upon request).
	- Generates the breakthrough time in each scenario in each flow regime.
	- Access to functions in dsenvsci repository is required to run this script.
	- Output file: tracer_combined_05032020.csv
- Calculating_Da_BG.ipynb:
	- Processes the simulation results to generate the Damkohler number.
	- Access to functions in dsenvsci repository is required to run this script.
	- Output file: Da_BG.csv
	
## Subdirectory Results
Information consistent throughout all the files:
- Column **Trial**: Internal identifier for scenario. 'H' refers to homogeneous domain. All other numbers refer to heterogeneous scenarios.
- Column **Variance**: Variance in the log hydraulic conductivity field to generate the spatial random field corresponding to the Trial.
- Column **Anisotropy**: Anisotropy to generate the spatial random field corresponding to the Trial.
- Column **Regime**: Flow regime as referred to in the manuscript.
- Column **Chem**: Reactive, or non-reactive, or microbial  species of concern.

The processed results with some metadata information:
- Da_BG.csv
	- Column **Conc_in**: Flux averaged concentration of species at the inlet (uM).
	- Column **Conc_out**: Flux averaged concentration of species at the outlet (uM).
	- Column **Normconc**: Ratio of Conc_in and Conc_out.
	- Column **base**: Normconc in the base case.
	- Column **k**: Pseudo first order reaction rate constant (d-1).
	- Column **tau**: Characteristic reaction time scale (d).
	- Column **Da**: Estimated Damkohler number.
	- Refer to massflux_comparison_steadystate_BG.csv for the description of the rest of the columns.
- tracer_combined_05032020.csv
	- Column **Time**: Breakthrough time (days)
	- Column **fraction**: Ratio of breakthrough time in heterogeneous scenario and that in homogeneous scenario (base case)