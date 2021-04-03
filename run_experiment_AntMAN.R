# load the libraries required for benchmarking

library("coda")
library("LaplacesDemon")
library("dplyr")
library("AntMAN")
library("functional")


#' plotting function to store plots of M and K 
#'
#'@param count_list The vector storing the count data to be plotted
#'@param x_lab String specifying the x_label of the plot
#'@param header String specifying the title of the plot

plotit <- function(count_list, x_lab, header){
  normalising_const <- sum(table(count_list))
  function() plot(table(count_list)/normalising_const, xlab=x_lab, ylab='post prob', main=header)}

#' run the benchmarks using AntMAN
#'
#' @param group String specifying the group list of the experiment e.g. A1, A2, A3, B1, B2, B3, C1, C2, C3
#' @param y_arg The y_data vector to be clustered
#' @param initial_clustering_ A vector specifying any initial cluster assignment of the observations. Default setting is NULL
#' @param kernel_hyperparams_ The kernel hyperparameters of the mixture model
#' @param components_prior_ The distribution of the prior on the number of components
#' @param weights_prior_ The distribution of the prior on the weights 
#' @param mcmc_params_ Parameters specifying the MCMC experimental set-up
#' @param thinning_interval Integer specifying the thinning interval of the MCMC object
#' @param no_of_iterations Number of experiments to average results over. Default is 10.
#' @param to_plot Boolean indicating whether plots of the posterior M and K should be generated. If TRUE, list of M and K are returned to facilitate
#'        plotting. Default is FALSE.
#' @param is_lca Boolean indicating whether we want the results in the format compatible with the lca experiments. Default is FALSE.
#'  
#' @return the averaged results (mean, sd) of the MCMC experiments for AntMAN
#' 
#'  



run_benchmarks <- function(group, y_arg, initial_clustering_ = NULL,
                           kernel_hyperparams_, components_prior_, weights_prior_,
                           mcmc_params_,  thinning_interval, no_of_iterations = 10, 
                           to_plot = FALSE, is_lca = FALSE 
){
  
  #initialise empty dataframe to store results
  
  antman_results <- data.frame(mean_K = numeric(0),
                               mean_M = numeric(0),
                               ESS_K = numeric(0),
                               IAC_K = numeric(0),
                               ESS_M = numeric(0),
                               IAC_M = numeric(0),
                               Time_taken = numeric(0))
  
  if (to_plot){
    
    #initialise the lists to store the plots
    M_plot_list = list()
    K_plot_list = list()
    
  }
  
  if (is_lca){
    
    # initialise dataframe to store the prob associated with the number of clusters
    # assume max number of clusters is 10
    prob_no_of_clusters <- data.frame(matrix(ncol=10, nrow=0))
    colnames(prob_no_of_clusters) <- c(1:10)
    
    
    # initialise dataframe to store the prob associated with each cluster, given the most
    # probable number of clusters
    prob_per_cluster <- data.frame(matrix(ncol=10, nrow=0))
    colnames(prob_per_cluster) <- c(1:10)
  }
  
  
  for (i in (1:no_of_iterations)){
    
    set.seed(i)
    
    # start the timer
    start_time <- proc.time()
    
    fit <- AM_mcmc_fit(
      y = y_arg, initial_clustering = initial_clustering_,
      mix_kernel_hyperparams = kernel_hyperparams_,
      mix_components_prior = components_prior_,
      mix_weight_prior = weights_prior_,
      mcmc_parameters = mcmc_params_)
    
    #end the timer
    end_time <- proc.time()
    time_taken <- getElement(end_time - start_time, "elapsed")
    
    
    #extract results
    component_cluster_no <- as.data.frame(AM_extract(fit, c('K','M')))
    
    mean_results <- colMeans(component_cluster_no)
    mean_K <- getElement(mean_results, "K")
    mean_M <- getElement(mean_results, "M")
    
    #compute ESS
    chain <- mcmc.list(mcmc(cbind(component_cluster_no), thin = thinning_interval))
    ESS_total <- effectiveSize(chain)
    
    ESS_K <- getElement(ESS_total, "K")
    ESS_M <- getElement(ESS_total, "M")
    
    #Compute IAC
    
    IAC_total <- apply(component_cluster_no, MARGIN = 2, FUN = IAT)
    IAC_K <- getElement(IAC_total, "K")
    IAC_M <- getElement(IAC_total, "M")
    
    #append results to dataframe
    
    
    antman_results[i,] <- c(mean_K, mean_M, ESS_K, IAC_K,
                            ESS_M, IAC_M, time_taken)
    
    
    if (to_plot){
      group_iter <- paste(group, '.%d', sep="") 
      header_ = sprintf(group_iter, i)
      
      M_plot_list[[i]] <- component_cluster_no$M
      K_plot_list[[i]] <- component_cluster_no$K
    }
    
    
    
    if (is_lca){
      
      # distribution of no of clusters
      G <- length(fit$K)
      prob_no_of_clusters_ <- as.data.frame.matrix(rbind(table(fit$K)/G))
      prob_no_of_clusters <- bind_rows(prob_no_of_clusters, prob_no_of_clusters_)
      
      # to obtain prob per cluster, given most probable no. of clusters
      result = AM_binder(fit)
      hatc <- result$Labels
      hatk <- length(unique(hatc))
      prob_per_cluster_ <- as.data.frame.matrix(rbind(table(hatc)/n)) 
      prob_per_cluster <- bind_rows(prob_per_cluster, prob_per_cluster_)
    }
  }
  
  if (is_lca){
    
    #replace NA values before returning
    prob_no_of_clusters[is.na(prob_no_of_clusters)] <- 0
    prob_per_cluster[is.na(prob_per_cluster)] <- 0
    
    return(list(antman_results, prob_no_of_clusters, prob_per_cluster))
  }
  
  if (to_plot){
    return(list(antman_results, M_plot_list, K_plot_list))
  }
  
  return(antman_results)
}



