# Modelling the effects of diurnal temperature variation on malaria infection dynamics in mosquitoes.

## Authors: Isaac J. Stopard, Antoine Sanou, Eunho Suh, Lauren J. Cator, Matthew B. Thomas, W. Moussa Guelbéogo, N'Falé Sagnon, Ben Lambert & Thomas S. Churcher 

This GitHub repository provides the code for the paper. We ran the model using R version 4.2.1.

#### Last updated: 15/07/2024

#### scripts 

:one: *read_libraries_data.R* - script to source commonly used data and load packages.

:two: *vector_species_comparison_temp_only.R* - script to compare the extrinsic incubation period (EIP) and human-to-mosquito transmission probability (HMTP) estimates. Generates Figure 2.

:three: *variance_comparison.R* - script to estimate the Erlang distribution shape parameter.

:four: *run_SMFA.R* - script to fit the model of mosquito infection dynamics during standard membrane feeding assays.

:five: *run_SMFA.R* - script to fit the model of mosquito infection dynamics during standard membrane feeding assays.

:six: *format_temp_data.R* - script to wrangle and interpolate the temperature data.

:seven: *microclimate.R* - script to predict the EIP and HMTP in Tiefora, Burkina Faso, using the best fitting model.

:eight: *spz_prev.R* - script to fit the generalised additive models to the previously published human biting rate data.

:nine: *EIR_fit.R* - script to fit the malaria transmission model to the sporozoite prevalence data from Tiefora.

#### helper scripts

:one: *read_bt_data.R* - script to read and format the biting time data.

:two: *data_functions.R*

:three: *model_functions.R*

### Notes

:warning: Please note that outputs are not saved on the Github repository.

:warning: This code is released with no support and we cannot endorse any outputs from the model other than those we generate.

:warning: This model is under development, so code may change without notice.

:warning: No liability is accepted by the authors.
