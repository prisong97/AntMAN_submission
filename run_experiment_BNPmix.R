#load the libraries required for benchmarking

library("BNPmix")
library("coda")

#' run the benchmarks using BNPmix
#'
#' @param y_arg The y_data vector to be clustered
#' @param prior_list A list detailing the type and hyperparameters of the model (refer to BNPmix documentation for more information)
#' @param mcmc_list Parameters specifying the MCMC experimental set-up
#' @param thinning_interval Integer specifying the thinning interval of the MCMC object
#' @param no_of_iterations Number of experiments to average results over. Default is 10.
#' @param summary_stats Boolean indicating whether summary information of K should be returned. Default is TRUE
#'  
#' @return the averaged results (mean, sd) of the MCMC experiments for AntMAN
#' 
#'  
run_benchmarks_bnpmix <- function(y_arg, prior_list, mcmc_list, thinning_interval, no_of_iterations, summary_stats = TRUE){
  
  #initialise empty dataframe to store results
  
  if (summary_stats == TRUE){
    bnpmix_results <- data.frame(mean_K = numeric(0),
                                 ESS_K = numeric(0),
                                 IAC_K = numeric(0),
                                 Time_taken = numeric(0))}
  else (bnpmix_results <- list())
  
  for (i in (1: no_of_iterations)){
    
    set.seed(i)
    
    # start the timer
    start_time <- proc.time()
    
    est_model <- PYdensity(y_arg, mcmc = mcmc_list,
                           prior = prior_list,
                           output = list(out_param = TRUE, mcmc_dens = TRUE))
    
    # end the timer
    end_time <- proc.time()
    time_taken <- getElement(end_time - start_time, "elapsed")
    
    #manually thin the MCMC chain 
    numClusters <- vapply(est_model$mean, function(x) {length(x)}, numeric(1))
    numClusters <- numClusters[seq(1, length(numClusters), by = thinning_interval)]
   
     chain <- mcmc.list(mcmc(cbind(NumClusters = numClusters), thin = thinning_interval))
    
    if (summary_stats == TRUE){
      #obtain the results of the Gibbs Sampler
      ESS_K <- effectiveSize(chain)
      IAC_K <- IAT(numClusters)
      mean_K <- mean(numClusters)
      
      bnpmix_results[i, ] <- c(mean_K, ESS_K, IAC_K, time_taken)
    }
    else {bnpmix_results[[i]] <- numClusters}
  }
    return(bnpmix_results) 
}




