#' Poisson tests
#' 
#' @description
#' Make inference on one or two populations using Poisson distributed count data
#' 
#' @param x Number of events.  A vector of length one or two.
#' @param offset Time, area, etc. measured in the Poisson process.  NOTE: Do not 
#' take the log!
#' @param r optional.  If provided and inference is being made for 
#' a single population, \code{poisson_test_b} will return the posterior 
#' probability that the population rate is less than this value.
#' @param ROPE ROPE for rate ratio if inference is being made for two populations. 
#' Provide either a single value or a vector of length two.  If the former, 
#' the ROPE will be taken as (1/ROPE,ROPE).  If the latter, these will be 
#' the bounds of the ROPE.
#' @param prior Either "jeffreys" (Gamma(1/2,0)) or "flat" (Gamma(0.001,0.001)).  
#' This is ignored if prior_shape_rate is provided.
#' @param prior_shape_rate Vector of length two, giving the shape and rate parameters 
#' for the gamma distribution that will act as the prior on the population 
#' rates.
#' @param CI_level The posterior probability to be contained in the 
#' credible intervals.
#' @param plot logical.  Should a plot be shown?
#' @param seed Always set your seed!  (Unused for a single population rate)
#' @param mc_error The number of posterior draws will ensure that with 99% 
#' probability the bounds of the credible intervals of \eqn{\lambda_1/\lambda_2} 
#' will be within \eqn{\pm} \code{mc_error}. (Ignored for a single population 
#' rate.)
#' 
#' @details
#' 
#' The likelihood is
#' \deqn{
#'  y \sim Poi(\lambda t),
#' }
#' where \eqn{\lambda} is the rate, and \eqn{t} is the time or area observed 
#' and is given by the argument \code{offset}.  
#' 
#' The prior is given by 
#' \deqn{
#'  \lambda \sim \Gamma(a,b),
#' }
#' where \eqn{a} and \eqn{b} are given by the argument \code{prior_shape_rate}. 
#' If \code{prior_shape_rate} is missing and \code{prior = "jeffreys"}, 
#' then a Jeffrey's prior will be used, i.e., \eqn{\Gamma(0.5,0)} (improper), 
#' while if \code{prior = "flat"}, \eqn{\Gamma(0.001,0.001)} will be used.  
#' 
#' 
#' @examples
#' \donttest{
#' # One sample
#' poisson_test_b(x = 12)
#' ## You can compute the posterior probability that the rate is less than r
#' poisson_test_b(x = 12,
#'                r = 8)
#' 
#' # Two samples
#' poisson_test_b(x = c(12,20))
#' 
#' # Offsets can be included:
#' poisson_test_b(x = c(12,20),
#'                offset = c(10,9))
#' 
#' # Different priors can be used
#' poisson_test_b(x = c(12,20),
#'                prior = "flat")
#' poisson_test_b(x = c(12,20),
#'                prior_shape_rate = c(20,1.5))
#' }
#' 
#' 
#' @export

poisson_test_b = function(x,
                          offset,
                          r,
                          ROPE,
                          prior = "jeffreys",
                          prior_shape_rate,
                          CI_level = 0.95,
                          plot = TRUE,
                          seed = 1,
                          mc_error = 0.002){
  
  set.seed(seed)
  
  alpha_ci = 1.0 - CI_level
  
  if(missing(offset)) offset = rep(1.0,length(x))
  
  # Preliminary checks
  if(!(length(x) %in% 1:2)) 
    stop("x must be of length one or two.")
  if(length(x) != length(offset)) 
    stop("Offset must be of the same length as x.")
  
  
  # Prior distribution
  if(missing(prior_shape_rate)){
    prior = c("flat",
              "jeffreys")[pmatch(tolower(prior),
                                 c("flat",
                                   "jeffreys"))]
    
    if(prior == "flat"){
      message("Prior shape parameters were not supplied.\nA flat Gamma(0.001,0.001) prior will be used.")
      prior_shape_rate = rep(0.001,2)
    }
    if(prior == "jeffreys"){
      message("Prior shape parameters were not supplied.\nJeffrey's prior will be used.")
      prior_shape_rate = c(0.5,0.0)
    }
  }else{
    if(length(prior_shape_rate) != 2)
      stop("Length of prior_shape_rate must equal two, relating to the shape and rate of the gamma distribution.")
    if(any(prior_shape_rate <= 0))
      stop("Prior shape parameters must be positive.")
  }
  
  
  # One sample inference ----------------------------------------------------
  
  if(length(x) == 1){
    
    # Get posterior parameters
    post_shape_rate = 
      prior_shape_rate + 
      c(x,
        offset)
    
    # Compute results
    results = 
      list(x = 
             x,
           offset = 
             offset,
           posterior_mean = 
             post_shape_rate[1] / post_shape_rate[2],
           CI = 
             c(
               qgamma(0.5 * alpha_ci,
                      shape = post_shape_rate[1],
                      rate = post_shape_rate[2]),
               qgamma(1.0 - 0.5 * alpha_ci,
                      shape = post_shape_rate[1],
                      rate = post_shape_rate[2])
             )
      )
    
    # Print results
    cat("\n----------\n\nAnalysis of a single population rate using Bayesian techniques\n")
    cat("\n----------\n\n")
    cat(paste0("Number of events: ", x,"\n\n"))
    cat(paste0("Time/area base for event counts: ", offset,"\n\n"))
    cat(paste0("Prior used: Gamma(", 
               format(signif(prior_shape_rate[1], 3), 
                      scientific = FALSE),
               ",",
               format(signif(prior_shape_rate[2], 3), 
                      scientific = FALSE),
               ")\n\n"))
    cat(paste0("Posterior mean of the rate: ", 
               format(signif(results$posterior_mean, 3), 
                      scientific = FALSE),
               "\n\n"))
    cat(paste0(100 * CI_level,
               "% credible interval: (", 
               format(signif(results$CI[1], 3), 
                      scientific = FALSE),
               ", ",
               format(signif(results$CI[2], 3), 
                      scientific = FALSE),
               ")\n\n"))
    if(!missing(r)){
      results$Pr_less_than_r = 
        pgamma(r,
               shape = post_shape_rate[1],
               rate = post_shape_rate[2])
      cat(paste0("Probability that rate < ",
                 format(signif(r, 3), 
                        scientific = FALSE),
                 ": ",
                 format(signif(results$Pr_less_than_r, 3), 
                        scientific = FALSE),
                 "\n\n"))
    }
    cat("\n----------\n\n")
    
    
    # Plot (if requested)
    if(plot){
      
      results$rate_plot = 
        tibble::tibble(x = seq(qgamma(0.005,
                                      shape = post_shape_rate[1],
                                      rate = post_shape_rate[2]),
                               qgamma(0.995,
                                      shape = post_shape_rate[1],
                                      rate = post_shape_rate[2]),
                               l = 50)) |> 
        ggplot(aes(x=x)) +
        stat_function(fun = 
                        function(x){
                          dgamma(x,
                                 shape = prior_shape_rate[1],
                                 rate = prior_shape_rate[2])
                        },
                      aes(color = "Prior"),
                      linewidth = 2) + 
        
        stat_function(fun = 
                        function(x){
                          dgamma(x,
                                 shape = post_shape_rate[1],
                                 rate = post_shape_rate[2])
                        },
                      aes(color = "Posterior"),
                      linewidth = 2) + 
        
        scale_color_manual(values = c("Prior" = "#440154FF", 
                                      "Posterior" = "#FDE725FF")) +
        theme_classic(base_size = 15) +
        xlab("") + 
        ylab("") + 
        labs(color = "Distribution") + 
        ggtitle("Population rate")
      
      print(results$rate_plot)
      
    }
    
    # Add posterior parameters to returned object
    results$posterior_parameters = 
      c(shape_1 = post_shape_rate[1],
        shape_2 = post_shape_rate[2])
    
    invisible(results)
  }else{#End: One sample inference
    
    # Two sample inference ----------------------------------------------------
    
    # Get ROPE
    if(missing(ROPE)){
      ROPE = c(1.0 / 1.125, 1.125)
      # From Kruchke (2018) on rate ratios from FDA <1.25. (Use half of small effect size for ROPE, hence 0.25/2) 
    }else{
      if(length(ROPE) > 2) stop("ROPE must be given as an upper bound, or given as both lower and upper bounds.")
      if((length(ROPE) > 1) & (ROPE[1] >= ROPE[2])) stop("ROPE lower bound must be smaller than ROPE upper bound")
      if(length(ROPE) == 1) ROPE = c(1.0 / ROPE, ROPE)
    }
    
    
    # Get posterior parameters
    post_shape_rate = 
      rbind(prior_shape_rate,
            prior_shape_rate) + 
      cbind(x,
            offset)
    
    # Get posterior draws
    ## Get preliminary draws
    lambda1_draws = 
      rgamma(500,
             shape = post_shape_rate[1,1],
             rate = post_shape_rate[1,2])
    lambda2_draws = 
      rgamma(500,
             shape = post_shape_rate[2,1],
             rate = post_shape_rate[2,2])
    fhat = 
      density(lambda1_draws / lambda2_draws,
              from = .Machine$double.eps)
    n_draws = 
      0.5 * alpha_ci * (1.0 - 0.5 * alpha_ci) *
      (
        qnorm(0.5 * (1.0 - 0.99)) / 
          mc_error /
          fhat$y[which.min(abs(fhat$x - 
                                 quantile(lambda1_draws / lambda2_draws, 0.5 * alpha_ci)))]
      )^2 |> 
      round()
    
    ## Finish posterior draws
    lambda1_draws = 
      c(lambda1_draws,
        rgamma(n_draws - length(lambda1_draws),
               shape = post_shape_rate[1,1],
               rate = post_shape_rate[1,2]))
    lambda2_draws = 
      c(lambda2_draws,
        rgamma(n_draws - length(lambda2_draws),
               shape = post_shape_rate[2,1],
               rate = post_shape_rate[2,2]))
    
    
    rate_ratios = 
      lambda1_draws / lambda2_draws
    
    
    # Find CI for rate ratios
    CI_bounds = 
      quantile(lambda1_draws / lambda2_draws,
               c(0.5 * alpha_ci,
                 1.0 - 0.5 * alpha_ci))
    
    # Compute results
    results = 
      list(x = x,
           offset = offset,
           posterior_mean_pop1 = 
             post_shape_rate[1,1] / post_shape_rate[1,2],
           posterior_mean_pop2 = 
             post_shape_rate[2,1] / post_shape_rate[2,2],
           CI_pop1 = 
             c(
               qgamma(0.5 * alpha_ci,
                      shape = post_shape_rate[1,1],
                      rate = post_shape_rate[1,2]),
               qgamma(1.0 - 0.5 * alpha_ci,
                      shape = post_shape_rate[1,1],
                      rate = post_shape_rate[1,2])
             ),
           CI_pop2 = 
             c(
               qgamma(0.5 * alpha_ci,
                      shape = post_shape_rate[2,1],
                      rate = post_shape_rate[2,2]),
               qgamma(1.0 - 0.5 * alpha_ci,
                      shape = post_shape_rate[2,1],
                      rate = post_shape_rate[2,2])
             ),
           CI_lambda1_over_lambda2 = 
             c(CI_bounds[1],CI_bounds[2]),
           Pr_rateratio_in_ROPE = 
             mean( (rate_ratios > ROPE[1]) & 
                     (rate_ratios < ROPE[2]) )
      )
    
    # Print results
    cat("\n----------\n\nAnalysis of two population rates using Bayesian techniques\n")
    cat("\n----------\n\n")
    cat(paste0("Number of events: Population 1 = ", 
               x[1],
               "; Population 2 = ",
               x[2],
               "\n\n"))
    cat(paste0("Time/area base for event counts: Population 1 = ",
               offset[1],
               "; Population 2 = ",
               offset[2],
               "\n\n"))
    cat(paste0("Prior used: Gamma(", 
               format(signif(prior_shape_rate[1], 3), 
                      scientific = FALSE),
               ",",
               format(signif(prior_shape_rate[2], 3), 
                      scientific = FALSE),
               ")\n\n"))
    cat(paste0("Posterior mean: Population 1 = ", 
               format(signif(results$posterior_mean_pop1, 3), 
                      scientific = FALSE),
               "; Population 2 = ",
               format(signif(results$posterior_mean_pop2, 3), 
                      scientific = FALSE),
               "\n\n"))
    cat(paste0(100 * CI_level,
               "% credible interval: Population 1 = (", 
               format(signif(results$CI_pop1[1], 3), 
                      scientific = FALSE),
               ", ",
               format(signif(results$CI_pop1[2], 3), 
                      scientific = FALSE),
               "); Population 2 = (",
               format(signif(results$CI_pop2[1], 3), 
                      scientific = FALSE),
               ", ",
               format(signif(results$CI_pop2[2], 3), 
                      scientific = FALSE),
               ")\n\n"))
    cat(paste0(100 * CI_level,
               "% credible interval: (Population 1) / (Population 2) = (", 
               format(signif(results$CI_lambda1_over_lambda2[1], 3), 
                      scientific = FALSE),
               ", ",
               format(signif(results$CI_lambda1_over_lambda2[2], 3), 
                      scientific = FALSE),
               ")\n\n"))
    cat(paste0("Probability that the rate ratio (pop 1 vs. pop 2) is in the ROPE, defined to be (",
               format(signif(ROPE[1], 3), 
                      scientific = FALSE),
               ",",
               format(signif(ROPE[2], 3), 
                      scientific = FALSE),
               ") = ",
               format(signif(results$Pr_rateratio_in_ROPE, 3), 
                      scientific = FALSE),
               "\n\n")) 
    cat("\n----------\n\n")
    
    
    # Plot (if requested)
    if(plot){
      results$rate_plot = 
        tibble::tibble(x = seq(min(qgamma(0.005,
                                          shape = post_shape_rate[,1],
                                          rate = post_shape_rate[,2])),
                               max(qgamma(0.995,
                                          shape = post_shape_rate[,1],
                                          rate = post_shape_rate[,2])),
                               l = 50)) |> 
        ggplot(aes(x=x)) +
        stat_function(fun = 
                        function(x){
                          dgamma(x,
                                 shape = prior_shape_rate[1],
                                 rate = prior_shape_rate[2])
                        },
                      aes(color = "Prior"),
                      linewidth = 2) + 
        stat_function(fun = 
                        function(x){
                          dgamma(x,
                                 shape = post_shape_rate[1,1],
                                 rate = post_shape_rate[1,2])
                        },
                      aes(color = "Posterior (Pop1)"),
                      linewidth = 2) + 
        stat_function(fun = 
                        function(x){
                          dgamma(x,
                                 shape = post_shape_rate[2,1],
                                 rate = post_shape_rate[2,2])
                        },
                      aes(color = "Posterior (Pop2)"),
                      linewidth = 2) + 
        scale_color_manual(values = c("Prior" = "#440154FF", 
                                      "Posterior (Pop1)" = "#21908CFF", 
                                      "Posterior (Pop2)" = "#FDE725FF")) +
        theme_classic(base_size = 15) +
        xlab("") + 
        ylab("") + 
        labs(color = "Distribution") + 
        ggtitle("Population rate")
      
      print(results$rate_plot)
      
    }
    
    # Add posterior parameters to returned object
    results$posterior_parameters = list()
    results$posterior_parameters$population_1 = 
      c(shape_1 = post_shape_rate[1,1],
        shape_2 = post_shape_rate[1,2])
    results$posterior_parameters$population_2 = 
      c(shape_1 = post_shape_rate[2,1],
        shape_2 = post_shape_rate[2,2])
    
    
    invisible(results)
  }#End: 2 sample inference
}
