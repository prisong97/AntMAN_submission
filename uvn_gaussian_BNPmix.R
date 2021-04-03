########################################
source("run_experiment_BNPmix.R")
########################################

bnp_univariate_gaussian_results <- c()
alpha_list <- c(0.160267,1.17344, 2.98651)
summary_stats_val = FALSE


for (i in (1: length(alpha_list))){
  
  alpha <- alpha_list[i]
  
  prior_list <- list(model = 'LS', strength= alpha,
                     m0 = 20.8315, k0= 0.01, a0 = 4/2, b0 = 4/(2*4))
  mcmc_list <- list(niter = 55000, nburn = 5000, method = 'ICS')
  
  experiment_results <- run_benchmarks_bnpmix(y_uvn, prior_list, mcmc_list, thinning_interval=10, no_of_iterations=10, summary_stats = summary_stats_val)
  
  if (summary_stats_val == TRUE){
    experimental_results_agg <- sapply(experiment_results, function(x) c("Mean" = mean(x), "std dev" = sd(x)))
    
    bnp_univariate_gaussian_results[[i]] <- experimental_results_agg}
  
  else {bnp_univariate_gaussian_results[[i]] <- experiment_results}
}
