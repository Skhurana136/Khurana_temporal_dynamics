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
- biomassdata_generation.py
	- Processes simulation results to generate average mass contribution of each microbial subpopulation.
	- Combines this with breakthrough time for each simulation scenario.
	- Output file: biomass_Original_complete.csv.
- biomassdata_comparison.py
	- Normalizes mean mass contribution of each microbial subpopulation with that in steady state conditions.
	- Output file: biomass_comparison_Original_complete.csv.
- massfluxdata_generation.py
	- Processes simulation results to derive the relative change in mass flux for all chemical species.
	- Combines change in mass flux with tracer breakthrough time in each simulation scenario.
	- Output file: massflux_complete_28022021.csv.
- massfluxdata_comparison.py
	- Normalizes change in mass flux with that in homogeneous domain at steady state conditions.
	- Normalizes change in mass flux with that in steady state conditions.
	- Output file: massflux_comparison_Original_complete_16062021.csv
- temp_animation.py
	- Script used to generate a video animation to visualize temporal dynamics in multiple heterogeneous domains.
	- Output file: spatio-temp-dynamics.mp4
- mass_flux_timeseries_generation_analysis.py
	- Processes simulation results to derive the responsiveness of the chemical species in the system and associated backward traceability.
	- Output files
		- normalized_sensitivity_19022021.csv
		- mass_flux_sensitivity_generalized_19022021.csv
		- biomass_sensitivity_generalized_19022021.csv
		- crosschor_memory_chem.csv
- biomass_timeseries_generation_analysis.py
	- Processes simulation results to assess responsiveness, cross-correlation and backward traceability of microbes in the system.
	- Output files
		- Normalized_RMSamplitude_biomass.csv
		- crosschor_memory_biomass.csv
	
## Subdirectory Results
Information consistent throughout all the files:
- Column **Trial**: Internal identifier for scenario. 'H' refers to homogeneous domain. All other numbers refer to heterogeneous scenarios.
- Column **Variance**: Variance in the log hydraulic conductivity field to generate the spatial random field corresponding to the Trial.
- Column **Anisotropy**: Anisotropy to generate the spatial random field corresponding to the Trial.
- Column **Regime**: Flow regime as referred to in the manuscript.
- Column **Chem**: Reactive, or non-reactive, or microbial  species of concern.

The processed results with some metadata information:
- CSV files:
	- biomass_comparison_Original_complete - Aggregated results compared with base case
	- biomass_Original_complete - Aggregated results
	- biomass_sensitivity_generalized_19022021 - Responsiveness of biomass normalized by base case
	- crosschor_memory_biomass - Time lag correlation analysis
	- crosschor_memory_chem - Time lag correlation analysis
	- Da_29012021_95pcloss - Damkohler number at steady state conditions
	- Headat inlet - Variation of groundwater head at inlet
	- mass_flux_sensitivity_generalized_19022021 - Responsiveness of chemical species normalized by base case
	- massflux_comparison_Original_complete - Aggregated mass removal over 15 years compared with base case
	- massflux_comparison_Original_complete_16062021 - Aggregated mass removal over 15 years compared with base case
	- massflux_Original_complete - Aggregated mass removal of 15 years
	- Normalized_RMSamplitude_biomass - Responsiveness of biomass
	- Normalized_RMSamplitude_biomass_190022021 - Responsiveness of biomass
	- Normalized_RMSamplitude_chem - Responsiveness of chemical species
	- Normalized_RMSamplitude_chem_withloc - Responsiveness of chemical species
	- normalized_sensivity_19022021 - Responsiveness of chemical species
	- tracer_combined_05032020 - Tracer breakthrough at steady state conditions
- HDF5 files:
- Temporal_analysis_biomass - Time series of biomass of all simulation scenarios
- Temporal_analysis_biomass_ratiopops - TIme series of biomass subpopulation ratios
- Temporal_analysis_full_Dat - Time series of chemical concentration leaving the domain divided on the basis of chemical regime (Da)
- Temporal_analysis_full_data - Time series of chemical concentrations leaving the domain of all scenarios
- Temporal_analysis_regimes_Dat - Time series of chemical concentration leaving the domain divided on the basis of chemical regime (Da)