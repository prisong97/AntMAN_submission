########################################################
library("devtools")
load_all('/Users/Priscilla/Desktop/antman/AntMAN/AntMAN')

# load external file 
source("run_experiment_AntMAN.R")

########################################################

data(galaxy, package = "AntMAN")
y_uvn <- galaxy

#initialise the distributions with hyperparameters specified in the infinity paper
mixture_uvn_params <- AM_mix_hyperparams_uninorm (m0=20.8315, k0=0.01,
                                                  nu0=4/2, sig02= 4/(2*4))
mcmc_params <- AM_mcmc_parameters(niter=55000, burnin=5000, thin=10)

#specify the variables to experiment with
lambda_gamma_pair = c("100, 2e-4", "10, 2e-3", "1, 10e-2",
                      "100, 1e-2", "10, 0.143", "5, 0.5",
                      "1000, 2.8e-3", "100, 3.2e-2", "10, 1.8")

#specifying the group list for plotting purposes
group_list = c("A1"  , "A2" , "A3",
               "B1", "B2", "B3",
               "C1", "C2", "C3")


#initialise list to store all the experimental results
univariate_antman_gaussian_results <- list()


for (i in (1: length(lambda_gamma_pair))){
  
  #retrieve group list, and corresponding hyperparameters
  group = group_list[i]
  pair = lambda_gamma_pair[i]
  pair_df = sapply(strsplit(pair, ","), function(x) {x <- as.numeric(x)})
  lambda = pair_df[1]
  gamma = pair_df[2]
  
  # fix Lambda and gamma
  components_prior <- AM_mix_components_prior_pois (Lambda=lambda)
  weights_prior <- AM_mix_weights_prior_gamma(gamma=gamma)
  
  experiment_results <- run_benchmarks(group, y_uvn, initial_clustering_ = NULL, mixture_uvn_params, components_prior, weights_prior, mcmc_params, thinning_interval = 10, no_of_iterations = 10, to_plot = FALSE, is_lca = FALSE)
  experimental_results_agg <- sapply(experiment_results[[1]], function(x) c("Mean" = mean(x), "std dev" = sd(x)))
  
  univariate_antman_gaussian_results[[i]] <- experimental_results_agg
}

