# HorizonExperimentsAnalysis
This repository contains all code for for analyses in the paper "Exploratory choice reflects the future value of information", to be published in the journal Decision in 2017. This repository will be updated to import the experiment data from the Decision website when it becomes avaiable.

Code for running the experiment can be found [here](https://github.com/NYUCCL/HorizonExperiments).

The files `analysis_1mushroom.r` and `analysis_4mushroom.r` contain non-Stan-related analyses and graphs of the raw data.

`fitstan_1mushroom.r`, `fitstan_4mushroom.r` and `fitstan_1mushroom_freehorizon.r` contain code to run the Stan models in
`stan_models` on the experment data.

`stan_analysis.r` contains code to analyze and graph the samples from the fitted Stan models.

`optimal_models.R` and `optimal_sim_graphs.R` contain code to calculate the
optimal choice for each information set and horizon, and to simulate and graph
optimal performance.
