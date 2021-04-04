This repository contains the scripts that were used to produce the tables in the paper: "AntMAN: Anthology of Mixture ANalysis tools". 

Specifically, we produced two sets of tables for:
- Benchmarking AntMAN against BNPmix for the univariate Gaussian model. Results of samplers implemented in both AntMAN and BNPmix are returned.
- Illustrating the results of AntMAN for the multivariate Bernoulli sampler.  

### How to reproduce the results

The main scripts to be run are:
- uvn_gaussian_AntMAN.R
- uvn_gaussian_BNPmix.R
- mvn_bernoulli_AntMAN.R

run_experiment_AntMAN.R and run_experiment_BNPmix.R contain helper functions that allow us to reproduce the documented experimental results. Ensure that these two scripts are placed in the same directory as the main scripts to be run. 
