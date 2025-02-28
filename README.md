# Modelling the effects of diurnal temperature variation on malaria infection dynamics in mosquitoes.

## Authors: Isaac J. Stopard, Antoine Sanou, Eunho Suh, Lauren J. Cator, Matthew B. Thomas, W. Moussa Guelbéogo, N'Falé Sagnon, Ben Lambert & Thomas S. Churcher 

This GitHub repository provides the code for the paper. We ran the model using R version 4.2.1.

#### Last updated: 27/02/2025

#### scripts 

:one: *1_run_SMFA.R* - script to fit the model of mosquito infection dynamics during standard membrane feeding assays.

:two: *2_vector_species_comparison_temp_only.R* - script to compare the extrinsic incubation period (EIP) and human-to-mosquito transmission probability (HMTP) estimates. Generates Figure 2.

:three: *3_variance_comparison.R* - script to estimate the Erlang distribution shape parameter.

:four: *4_temp_BF.R* - script to predict the EIP and HMTP in Tiefora, Burkina Faso, using the best fitting model.

:five: *5_spz_prev.R* - script to fit the generalised additive models to the previously published human biting rate data.

:six: *6_EIR_fit.R* - script to fit the Ross-Macdonald type malaria transmission model to the sporozoite prevalence data from Tiefora.

#### helper scripts

:one: *read_libraries_data.R* - script to source commonly used data and load packages.

:two: *read_bt_data.R* - script to read and format the biting time data.

:three: *data_wrangling_functions.R* - functions used across scripts to wrangle SMFA data and extract values from the model fits.

:four: *model_functions.R* - SMFA ODE model and helper functions.

:five: *format_temp_data.R* - script to wrangle and interpolate the temperature data.

### Notes

:warning: Please note that outputs are not saved on the Github repository.

:warning: This code is released with no support and we cannot endorse any outputs from the model other than those we generate.

:warning: This model is under development, so code may change without notice.

:warning: No liability is accepted by the authors.
