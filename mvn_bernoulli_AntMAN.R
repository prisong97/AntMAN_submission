########################################################
library("devtools")
load_all('/Users/Priscilla/Desktop/antman/AntMAN/AntMAN')

# load external file 
source("run_experiment_AntMAN.R")

########################################################

#load the data
data(carcinoma, package="AntMAN")
y_mvb <- carcinoma
n <- dim(y_mvb)[1]
d <- dim(y_mvb)[2]

# specify the kernel and priors
mixture_mvb_params <- AM_mix_hyperparams_multiber(a0=rep(1,d),b0= rep(1,d))
weights_prior <- AM_mix_weights_prior_gamma(init=5, a=1, b=1)
components_prior <- AM_mix_components_prior_negbin(R=1, init_P=0.1,a_P=1,b_P =1)

# specify the MCMC parameters
mcmc_params <- AM_mcmc_parameters(niter=255000, burnin=5000, thin=50,
                                  verbose=1)

# run multivariate bernoulli model
fit  <- run_benchmarks("group", y_mvb, initial_clustering_ = NULL, mixture_mvb_params, components_prior, weights_prior, mcmc_params, thinning_interval = 50, no_of_iterations = 10, to_plot = FALSE, is_lca = TRUE)