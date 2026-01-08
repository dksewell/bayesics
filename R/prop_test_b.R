#' Bayesian test of Equal or Given Proportions
#' @aliases binom_test_b
#' @aliases prop_test_b
#' 
#' @description
#' \code{prop_test_b} either makes inference on a single population 
#' proportion, or else compares two population proportions.  
#' \code{binom_test_b} is the same as \code{prop_test_b}.
#' 
#' 
#' @details
#' The likelihood is given by
#' \deqn{
#'  y \sim \text{Binom}(n,p),
#' }
#' and the prior on \eqn{p} is
#' \deqn{
#'  p \sim Beta(a,b),
#' }
#' where \eqn{a} and \eqn{b} are given by the argument \code{prior_shapes}.  If 
#' \code{prior_shapes} is missing and \code{prior = "jeffreys"}, then a 
#' Jeffreys prior will be used (\eqn{Beta(1/2,1/2)}), and if 
#' \code{prior = "uniform"}, then a uniform prior will be used (\eqn{Beta(1,1)}).
#' 
#' 
#' @param n_successes integer/numeric vector of length 1 (for 1 population) or 
#' 2 (for 2 populations) providing the number of "successes"
#' @param n_failures Similar to n_successes, but for failures.  Only provide this 
#' OR n_total.
#' @param n_total Similar to n_successes, but for total number of trials.  Only provide this 
#' OR n_failures.
#' @param p optional.  If provided and inference is being made for 
#' a single population, \code{prop_test_b} will return the posterior 
#' probability that the population proportion is less than this value.
#' @param predict_for_n Number in a future trial.  If missing, \code{prop_test_b} 
#' will use the observed number of trials.
#' @param ROPE ROPE for odds ratio if inference is being made for two populations. 
#' Provide either a single value or a vector of length two.  If the former, 
#' the ROPE will be taken as (1/ROPE,ROPE).  If the latter, these will be 
#' the bounds of the ROPE.
#' @param prior Either "jeffreys" (Beta(1/2,1/2)) or "uniform" (Beta(1,1)).  
#' This is ignored if prior_shapes is provided.
#' @param prior_shapes Vector of length two, giving the shape parameters 
#' for the beta distribution that will act as the prior on the population 
#' proportions.
#' @param CI_level,PI_level The posterior probability to be contained in the 
#' credible and prediction intervals respectively.
#' @param plot logical.  Should a plot be shown?
#' @param seed Always set your seed!  (Unused for a single population proportion.)
#' @param mc_error The number of posterior draws will ensure that with 99% 
#' probability the bounds of the credible intervals of \eqn{p_1 - p_2} will be 
#' within \eqn{\pm} \code{mc_error}. (Ignored for a single population proportion.)
#' 
#' 
#' @examples
#' # Single population
#' prop_test_b(14,
#'             19)
#' # or another way of the same thing;
#' prop_test_b(14,
#'             n_total = 14 + 19)
#' 
#' # A null value compared against can be added:
#' prop_test_b(14,
#'             19,
#'             p = 0.5)
#' 
#' # Two populations
#' prop_test_b(c(14,22),
#'             c(19,45))
#' # or equivalently
#' prop_test_b(c(14,22),
#'             n_total = c(14,22) + c(19,45))
#' 
#' @export

prop_test_b = function(n_successes,
                       n_failures,
                       n_total,
                       p,
                       predict_for_n,
                       ROPE,
                       prior = "jeffreys",
                       prior_shapes,
                       CI_level = 0.95,
                       PI_level = 0.95,
                       plot = TRUE,
                       seed = 1,
                       mc_error = 0.002){
  
  set.seed(seed)
  
  alpha_ci = 1.0 - CI_level
  alpha_pi = 1.0 - PI_level
  
  # Preliminary checks
  if(missing(n_successes)) 
    stop("Must provide the number of successes")
  if(missing(n_failures) & missing(n_total)) 
    stop("Must provide either number of failures or the total number of trials.")
  if( !missing(n_failures) & !missing(n_total) ) 
    if(!isTRUE(all.equal(n_failures,n_total - n_successes)))
      stop("Number of failures provided is inconsistent with total number of trials")
  if(missing(n_failures))
    n_failures = n_total - n_successes
  if(missing(n_total))
    n_total = n_failures + n_successes
  if(missing(predict_for_n)){
    predict_for_n = n_total
  }
  if(length(n_successes) != length(n_failures))
    stop("Length of n_successes must be the same length as n_successes/n_failures")
  if(length(n_successes) > 2)
    stop("The number of populations must be either 1 or 2.")
  
  
  # Prior distribution
  if(missing(prior_shapes)){
    prior = c("uniform",
              "jeffreys")[pmatch(tolower(prior),
                                 c("uniform",
                                   "jeffreys"))]
    
    if(prior == "uniform"){
      message("Prior shape parameters were not supplied.\nA uniform prior will be used.")
      prior_shapes = rep(1.0,2)
    }
    if(prior == "jeffreys"){
      message("Prior shape parameters were not supplied.\nJeffrey's prior will be used.")
      prior_shapes = rep(0.5,2)
    }
  }else{
    if(any(prior_shapes <= 0))
      stop("Prior shape parameters must be positive.")
  }
  
  
  # One sample inference ----------------------------------------------------
  
  if(length(n_successes) == 1){
    
    # Get posterior parameters
    post_shapes = 
      prior_shapes + 
      c(n_successes,
        n_failures)
    
    # Get CDF of predictive posterior (beta-binomial)
    cdf_probs = 
      extraDistr::pbbinom(0:predict_for_n,
                          predict_for_n,
                          post_shapes[1],
                          post_shapes[2])
    
    # Compute results
    results = 
      list(successes = 
             n_successes,
           failures = 
             n_failures,
           posterior_mean = 
             post_shapes[1] / sum(post_shapes),
           CI = 
             c(
               qbeta(0.5 * alpha_ci,
                     post_shapes[1],
                     post_shapes[2]),
               qbeta(1.0 - 0.5 * alpha_ci,
                     post_shapes[1],
                     post_shapes[2])
             ),
           PI = 
             c(
               max(c(0,which(cdf_probs <= 0.5 * alpha_pi))),
               min(c(predict_for_n,
                     1 + which(cdf_probs > 1.0 - 0.5 * alpha_pi)))
             )
      )
    
    # Print results
    cat("\n----------\n\nAnalysis of a single population proportion using Bayesian techniques\n")
    cat("\n----------\n\n")
    cat(paste0("Number of successes: ", n_successes,"\n\n"))
    cat(paste0("Number of failures: ", n_failures,"\n\n"))
    cat(paste0("Prior used: Beta(", 
               format(signif(prior_shapes[1], 3), 
                      scientific = FALSE),
               ",",
               format(signif(prior_shapes[2], 3), 
                      scientific = FALSE),
               ")\n\n"))
    cat(paste0("Posterior mean: ", 
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
    if(!missing(p)){
      results$Pr_less_than_p = 
        pbeta(p,post_shapes[1],post_shapes[2])
      cat(paste0("Probability that p < ",
                 format(signif(p, 3), 
                        scientific = FALSE),
                 ": ",
                 format(signif(results$Pr_less_than_p, 3), 
                        scientific = FALSE),
                 "\n\n"))
    }
    cat(paste0(100 * PI_level,
               "% prediction interval for ",
               ifelse(predict_for_n == n_total,
                      "another ",
                      ""),
               predict_for_n,
               " trials: (", 
               results$PI[1],
               ", ",
               results$PI[2],
               ")\n\n"))
    cat("\n----------\n\n")
    
    
    # Plot (if requested)
    if(plot){
      
      results$prop_plot = 
        tibble::tibble(x = seq(0.001,0.999,#seq(.Machine$double.eps,1.0 - .Machine$double.eps,
                               l = 50)) |> 
        ggplot(aes(x=.data$x)) +
        stat_function(fun = 
                        function(x){
                          dbeta(x,
                                prior_shapes[1],
                                prior_shapes[2])
                        },
                      aes(color = "Prior"),
                      linewidth = 2) + 
        
        stat_function(fun = 
                        function(x){
                          dbeta(x,
                                post_shapes[1],
                                post_shapes[2])
                        },
                      aes(color = "Posterior"),
                      linewidth = 2) + 
        
        scale_color_manual(values = c("Prior" = "#440154FF", 
                                      "Posterior" = "#FDE725FF")) +
        theme_classic(base_size = 15) +
        xlab("") + 
        ylab("") + 
        labs(color = "Distribution") + 
        ggtitle("Population proportion")
      
      print(results$prop_plot)
      
    }
    
    # Add posterior parameters to returned object
    results$posterior_parameters = 
      c(shape_1 = post_shapes[1],
        shape_2 = post_shapes[2])
    
    invisible(results)
  }else{#End: One sample inference
    
    # Two sample inference ----------------------------------------------------
    
    # Get ROPE
    if(missing(ROPE)){
      ROPE = c(1.0 / 1.125, 1.125)
      # From Kruchke (2018) on rate ratios from FDA <1.25. (Use half of small effect size for ROPE, hence 0.25/2) 
      #   Use the same thing for odds ratios.
    }else{
      if(length(ROPE) > 2) stop("ROPE must be given as an upper bound, or given as both lower and upper bounds.")
      if((length(ROPE) > 1) & (ROPE[1] >= ROPE[2])) stop("ROPE lower bound must be smaller than ROPE upper bound")
      if(length(ROPE) == 1) ROPE = c(1.0 / ROPE, ROPE)
    }
    
    
    # Get posterior parameters
    post_shapes = 
      rbind(prior_shapes,
            prior_shapes) + 
      cbind(n_successes,
            n_failures)
    
    # Get posterior draws
    ## Get preliminary draws
    p1_draws = 
      rbeta(500,
            post_shapes[1,1],
            post_shapes[1,2])
    p2_draws = 
      rbeta(500,
            post_shapes[2,1],
            post_shapes[2,2])
    fhat = 
      density(p1_draws - p2_draws,
              from = -1.0 + .Machine$double.eps,
              to = 1.0 - .Machine$double.eps)
    n_draws = 
      0.5 * alpha_ci * (1.0 - 0.5 * alpha_ci) *
      (
        qnorm(0.5 * (1.0 - 0.99)) / 
          mc_error /
          fhat$y[which.min(abs(fhat$x - 
                                 quantile(p1_draws - p2_draws, 0.5 * alpha_ci)))]
      )^2 |> 
      round()
    
    ## Finish posterior draws
    p1_draws = 
      c(p1_draws,
        rbeta(n_draws - length(p1_draws),
              post_shapes[1,1],
              post_shapes[1,2]))
    p2_draws = 
      c(p2_draws,
        rbeta(n_draws - length(p2_draws),
              post_shapes[2,1],
              post_shapes[2,2]))
    
    
    odds_ratios = 
      p1_draws / (1.0 - p1_draws) * (1.0 - p2_draws) / p2_draws
    
    
    # Get CDF of predictive posterior (beta-binomial)
    cdf_probs = 
      list(pop1 = 
             extraDistr::pbbinom(0:predict_for_n[1],
                                 predict_for_n[1],
                                 post_shapes[1,1],
                                 post_shapes[1,2]),
           pop2 = 
             extraDistr::pbbinom(0:predict_for_n[2],
                                 predict_for_n[2],
                                 post_shapes[2,1],
                                 post_shapes[2,2])
      )
    
    # Find CI for difference in prob
    CI_bounds = 
      quantile(p1_draws - p2_draws,
               c(0.5 * alpha_ci,
                 1.0 - 0.5 * alpha_ci))
    
    # Compute results
    results = 
      list(successes = 
             n_successes,
           failures = 
             n_failures,
           posterior_mean_pop1 = 
             post_shapes[1,1] / sum(post_shapes[1,]),
           posterior_mean_pop2 = 
             post_shapes[2,1] / sum(post_shapes[2,]),
           CI_pop1 = 
             c(
               qbeta(0.5 * alpha_ci,
                     post_shapes[1,1],
                     post_shapes[1,2]),
               qbeta(1.0 - 0.5 * alpha_ci,
                     post_shapes[1,1],
                     post_shapes[1,2])
             ),
           CI_pop2 = 
             c(
               qbeta(0.5 * alpha_ci,
                     post_shapes[2,1],
                     post_shapes[2,2]),
               qbeta(1.0 - 0.5 * alpha_ci,
                     post_shapes[2,1],
                     post_shapes[2,2])
             ),
           CI_p1_minus_p2 = 
             c(CI_bounds[1],CI_bounds[2]),
           Pr_oddsratio_in_ROPE = 
             mean( (odds_ratios > ROPE[1]) & 
                     (odds_ratios < ROPE[2]) ),
           PI_pop1 = 
             c(
               max(c(0,which(cdf_probs$pop1 <= 0.5 * alpha_pi))),
               min(c(predict_for_n[1],
                     1 + which(cdf_probs$pop1 > 1.0 - 0.5 * alpha_pi)))
             ),
           PI_pop2 = 
             c(
               max(c(0,which(cdf_probs$pop2 <= 0.5 * alpha_pi))),
               min(c(predict_for_n[length(predict_for_n)],
                     1 + which(cdf_probs$pop2 > 1.0 - 0.5 * alpha_pi)))
             )
      )
    
    # Print results
    cat("\n----------\n\nAnalysis of two population proportions using Bayesian techniques\n")
    cat("\n----------\n\n")
    cat(paste0("Number of successes: Population 1 = ", 
               n_successes[1],
               "; Population 2 = ",
               n_successes[2],
               "\n\n"))
    cat(paste0("Number of failures: Population 1 = ",
               n_failures[1],
               "; Population 2 = ",
               n_failures[2],
               "\n\n"))
    cat(paste0("Prior used: Beta(", 
               format(signif(prior_shapes[1], 3), 
                      scientific = FALSE),
               ",",
               format(signif(prior_shapes[2], 3), 
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
               "% credible interval: (Population 1) - (Population 2) = (", 
               format(signif(results$CI_p1_minus_p2[1], 3), 
                      scientific = FALSE),
               ", ",
               format(signif(results$CI_p1_minus_p2[2], 3), 
                      scientific = FALSE),
               ")\n\n"))
    cat(paste0("Probability that the odds ratio (pop 1 vs. pop 2) is in the ROPE, defined to be (",
               format(signif(ROPE[1], 3), 
                      scientific = FALSE),
               ",",
               format(signif(ROPE[2], 3), 
                      scientific = FALSE),
               ") = ",
               format(signif(results$Pr_oddsratio_in_ROPE, 3), 
                      scientific = FALSE),
               "\n\n")) 
    cat(paste0(100 * PI_level,
               "% prediction interval for ",
               ifelse(predict_for_n[1] == n_total[1],
                      "another ",
                      ""),
               predict_for_n[1],
               " trials for population 1: (", 
               results$PI_pop1[1],
               ", ",
               results$PI_pop1[2],
               ")\n\n"))
    cat(paste0(100 * PI_level,
               "% prediction interval for ",
               ifelse(predict_for_n[length(predict_for_n)] == n_total[2],
                      "another ",
                      ""),
               predict_for_n[length(predict_for_n)],
               " trials for population 2: (", 
               results$PI_pop2[1],
               ", ",
               results$PI_pop2[2],
               ")\n\n"))
    cat("\n----------\n\n")
    
    
    # Plot (if requested)
    if(plot){
      results$prop_plot = 
        tibble::tibble(x = seq(0.001,0.999,
                               l = 50)) |> 
        ggplot(aes(x=.data$x)) +
        stat_function(fun = 
                        function(x){
                          dbeta(x,
                                prior_shapes[1],
                                prior_shapes[2])
                        },
                      aes(color = "Prior"),
                      linewidth = 2) + 
        stat_function(fun = 
                        function(x){
                          dbeta(x,
                                post_shapes[1,1],
                                post_shapes[1,2])
                        },
                      aes(color = "Posterior (Pop1)"),
                      linewidth = 2) + 
        stat_function(fun = 
                        function(x){
                          dbeta(x,
                                post_shapes[2,1],
                                post_shapes[2,2])
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
        ggtitle("Population proportion")
      
      print(results$prop_plot)
      
    }
    
    # Add posterior parameters to returned object
    results$posterior_parameters = list()
    results$posterior_parameters$population_1 = 
      c(shape_1 = post_shapes[1,1],
        shape_2 = post_shapes[1,2])
    results$posterior_parameters$population_2 = 
      c(shape_1 = post_shapes[2,1],
        shape_2 = post_shapes[2,2])
    
    
    invisible(results)
  }#End: 2 sample inference
}

#' @export
binom_test_b = function(n_successes,
                        n_failures,
                        n_total,
                        p,
                        predict_for_n,
                        ROPE,
                        prior = "jeffreys",
                        prior_shapes,
                        CI_level = 0.95,
                        PI_level = 0.95,
                        plot = TRUE,
                        seed = 1,
                        mc_error = 0.002){
  
  set.seed(seed)
  
  alpha_ci = 1.0 - CI_level
  alpha_pi = 1.0 - PI_level
  
  # Preliminary checks
  if(missing(n_successes)) 
    stop("Must provide the number of successes")
  if(missing(n_failures) & missing(n_total)) 
    stop("Must provide either number of failures or the total number of trials.")
  if( !missing(n_failures) & !missing(n_total) ) 
    if(!isTRUE(all.equal(n_failures,n_total - n_successes)))
      stop("Number of failures provided is inconsistent with total number of trials")
  if(missing(n_failures))
    n_failures = n_total - n_successes
  if(missing(n_total))
    n_total = n_failures + n_successes
  if(missing(predict_for_n)){
    predict_for_n = n_total
  }
  if(length(n_successes) != length(n_failures))
    stop("Length of n_successes must be the same length as n_successes/n_failures")
  if(length(n_successes) > 2)
    stop("The number of populations must be either 1 or 2.")
  
  
  # Prior distribution
  if(missing(prior_shapes)){
    prior = c("uniform",
              "jeffreys")[pmatch(tolower(prior),
                                 c("uniform",
                                   "jeffreys"))]
    
    if(prior == "uniform"){
      message("Prior shape parameters were not supplied.\nA uniform prior will be used.")
      prior_shapes = rep(1.0,2)
    }
    if(prior == "jeffreys"){
      message("Prior shape parameters were not supplied.\nJeffrey's prior will be used.")
      prior_shapes = rep(0.5,2)
    }
  }else{
    if(any(prior_shapes <= 0))
      stop("Prior shape parameters must be positive.")
  }
  
  
  # One sample inference ----------------------------------------------------
  
  if(length(n_successes) == 1){
    
    # Get posterior parameters
    post_shapes = 
      prior_shapes + 
      c(n_successes,
        n_failures)
    
    # Get CDF of predictive posterior (beta-binomial)
    cdf_probs = 
      extraDistr::pbbinom(0:predict_for_n,
                          predict_for_n,
                          post_shapes[1],
                          post_shapes[2])
    
    # Compute results
    results = 
      list(successes = 
             n_successes,
           failures = 
             n_failures,
           posterior_mean = 
             post_shapes[1] / sum(post_shapes),
           CI = 
             c(
               qbeta(0.5 * alpha_ci,
                     post_shapes[1],
                     post_shapes[2]),
               qbeta(1.0 - 0.5 * alpha_ci,
                     post_shapes[1],
                     post_shapes[2])
             ),
           PI = 
             c(
               max(c(0,which(cdf_probs <= 0.5 * alpha_pi))),
               min(c(predict_for_n,
                     1 + which(cdf_probs > 1.0 - 0.5 * alpha_pi)))
             )
      )
    
    # Print results
    cat("\n----------\n\nAnalysis of a single population proportion using Bayesian techniques\n")
    cat("\n----------\n\n")
    cat(paste0("Number of successes: ", n_successes,"\n\n"))
    cat(paste0("Number of failures: ", n_failures,"\n\n"))
    cat(paste0("Prior used: Beta(", 
               format(signif(prior_shapes[1], 3), 
                      scientific = FALSE),
               ",",
               format(signif(prior_shapes[2], 3), 
                      scientific = FALSE),
               ")\n\n"))
    cat(paste0("Posterior mean: ", 
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
    if(!missing(p)){
      results$Pr_less_than_p = 
        pbeta(p,post_shapes[1],post_shapes[2])
      cat(paste0("Probability that p < ",
                 format(signif(p, 3), 
                        scientific = FALSE),
                 ": ",
                 format(signif(results$Pr_less_than_p, 3), 
                        scientific = FALSE),
                 "\n\n"))
    }
    cat(paste0(100 * PI_level,
               "% prediction interval for ",
               ifelse(predict_for_n == n_total,
                      "another ",
                      ""),
               predict_for_n,
               " trials: (", 
               results$PI[1],
               ", ",
               results$PI[2],
               ")\n\n"))
    cat("\n----------\n\n")
    
    
    # Plot (if requested)
    if(plot){
      
      results$prop_plot = 
        tibble::tibble(x = seq(0.001,0.999,#seq(.Machine$double.eps,1.0 - .Machine$double.eps,
                               l = 50)) |> 
        ggplot(aes(x=.data$x)) +
        stat_function(fun = 
                        function(x){
                          dbeta(x,
                                prior_shapes[1],
                                prior_shapes[2])
                        },
                      aes(color = "Prior"),
                      linewidth = 2) + 
        
        stat_function(fun = 
                        function(x){
                          dbeta(x,
                                post_shapes[1],
                                post_shapes[2])
                        },
                      aes(color = "Posterior"),
                      linewidth = 2) + 
        
        scale_color_manual(values = c("Prior" = "#440154FF", 
                                      "Posterior" = "#FDE725FF")) +
        theme_classic(base_size = 15) +
        xlab("") + 
        ylab("") + 
        labs(color = "Distribution") + 
        ggtitle("Population proportion")
      
      print(results$prop_plot)
      
    }
    
    # Add posterior parameters to returned object
    results$posterior_parameters = 
      c(shape_1 = post_shapes[1],
        shape_2 = post_shapes[2])
    
    invisible(results)
  }else{#End: One sample inference
    
    # Two sample inference ----------------------------------------------------
    
    # Get ROPE
    if(missing(ROPE)){
      ROPE = c(1.0 / 1.125, 1.125)
      # From Kruchke (2018) on rate ratios from FDA <1.25. (Use half of small effect size for ROPE, hence 0.25/2) 
      #   Use the same thing for odds ratios.
    }else{
      if(length(ROPE) > 2) stop("ROPE must be given as an upper bound, or given as both lower and upper bounds.")
      if((length(ROPE) > 1) & (ROPE[1] >= ROPE[2])) stop("ROPE lower bound must be smaller than ROPE upper bound")
      if(length(ROPE) == 1) ROPE = c(1.0 / ROPE, ROPE)
    }
    
    
    # Get posterior parameters
    post_shapes = 
      rbind(prior_shapes,
            prior_shapes) + 
      cbind(n_successes,
            n_failures)
    
    # Get posterior draws
    ## Get preliminary draws
    p1_draws = 
      rbeta(500,
            post_shapes[1,1],
            post_shapes[1,2])
    p2_draws = 
      rbeta(500,
            post_shapes[2,1],
            post_shapes[2,2])
    fhat = 
      density(p1_draws - p2_draws,
              from = -1.0 + .Machine$double.eps,
              to = 1.0 - .Machine$double.eps)
    n_draws = 
      0.5 * alpha_ci * (1.0 - 0.5 * alpha_ci) *
      (
        qnorm(0.5 * (1.0 - 0.99)) / 
          mc_error /
          fhat$y[which.min(abs(fhat$x - 
                                 quantile(p1_draws - p2_draws, 0.5 * alpha_ci)))]
      )^2 |> 
      round()
    
    ## Finish posterior draws
    p1_draws = 
      c(p1_draws,
        rbeta(n_draws - length(p1_draws),
              post_shapes[1,1],
              post_shapes[1,2]))
    p2_draws = 
      c(p2_draws,
        rbeta(n_draws - length(p2_draws),
              post_shapes[2,1],
              post_shapes[2,2]))
    
    
    odds_ratios = 
      p1_draws / (1.0 - p1_draws) * (1.0 - p2_draws) / p2_draws
    
    
    # Get CDF of predictive posterior (beta-binomial)
    cdf_probs = 
      list(pop1 = 
             extraDistr::pbbinom(0:predict_for_n[1],
                                 predict_for_n[1],
                                 post_shapes[1,1],
                                 post_shapes[1,2]),
           pop2 = 
             extraDistr::pbbinom(0:predict_for_n[2],
                                 predict_for_n[2],
                                 post_shapes[2,1],
                                 post_shapes[2,2])
      )
    
    # Find CI for difference in prob
    CI_bounds = 
      quantile(p1_draws - p2_draws,
               c(0.5 * alpha_ci,
                 1.0 - 0.5 * alpha_ci))
    
    # Compute results
    results = 
      list(successes = 
             n_successes,
           failures = 
             n_failures,
           posterior_mean_pop1 = 
             post_shapes[1,1] / sum(post_shapes[1,]),
           posterior_mean_pop2 = 
             post_shapes[2,1] / sum(post_shapes[2,]),
           CI_pop1 = 
             c(
               qbeta(0.5 * alpha_ci,
                     post_shapes[1,1],
                     post_shapes[1,2]),
               qbeta(1.0 - 0.5 * alpha_ci,
                     post_shapes[1,1],
                     post_shapes[1,2])
             ),
           CI_pop2 = 
             c(
               qbeta(0.5 * alpha_ci,
                     post_shapes[2,1],
                     post_shapes[2,2]),
               qbeta(1.0 - 0.5 * alpha_ci,
                     post_shapes[2,1],
                     post_shapes[2,2])
             ),
           CI_p1_minus_p2 = 
             c(CI_bounds[1],CI_bounds[2]),
           Pr_oddsratio_in_ROPE = 
             mean( (odds_ratios > ROPE[1]) & 
                     (odds_ratios < ROPE[2]) ),
           PI_pop1 = 
             c(
               max(c(0,which(cdf_probs$pop1 <= 0.5 * alpha_pi))),
               min(c(predict_for_n[1],
                     1 + which(cdf_probs$pop1 > 1.0 - 0.5 * alpha_pi)))
             ),
           PI_pop2 = 
             c(
               max(c(0,which(cdf_probs$pop2 <= 0.5 * alpha_pi))),
               min(c(predict_for_n[length(predict_for_n)],
                     1 + which(cdf_probs$pop2 > 1.0 - 0.5 * alpha_pi)))
             )
      )
    
    # Print results
    cat("\n----------\n\nAnalysis of two population proportions using Bayesian techniques\n")
    cat("\n----------\n\n")
    cat(paste0("Number of successes: Population 1 = ", 
               n_successes[1],
               "; Population 2 = ",
               n_successes[2],
               "\n\n"))
    cat(paste0("Number of failures: Population 1 = ",
               n_failures[1],
               "; Population 2 = ",
               n_failures[2],
               "\n\n"))
    cat(paste0("Prior used: Beta(", 
               format(signif(prior_shapes[1], 3), 
                      scientific = FALSE),
               ",",
               format(signif(prior_shapes[2], 3), 
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
               "% credible interval: (Population 1) - (Population 2) = (", 
               format(signif(results$CI_p1_minus_p2[1], 3), 
                      scientific = FALSE),
               ", ",
               format(signif(results$CI_p1_minus_p2[2], 3), 
                      scientific = FALSE),
               ")\n\n"))
    cat(paste0("Probability that the odds ratio (pop 1 vs. pop 2) is in the ROPE, defined to be (",
               format(signif(ROPE[1], 3), 
                      scientific = FALSE),
               ",",
               format(signif(ROPE[2], 3), 
                      scientific = FALSE),
               ") = ",
               format(signif(results$Pr_oddsratio_in_ROPE, 3), 
                      scientific = FALSE),
               "\n\n")) 
    cat(paste0(100 * PI_level,
               "% prediction interval for ",
               ifelse(predict_for_n[1] == n_total[1],
                      "another ",
                      ""),
               predict_for_n[1],
               " trials for population 1: (", 
               results$PI_pop1[1],
               ", ",
               results$PI_pop1[2],
               ")\n\n"))
    cat(paste0(100 * PI_level,
               "% prediction interval for ",
               ifelse(predict_for_n[length(predict_for_n)] == n_total[2],
                      "another ",
                      ""),
               predict_for_n[length(predict_for_n)],
               " trials for population 2: (", 
               results$PI_pop2[1],
               ", ",
               results$PI_pop2[2],
               ")\n\n"))
    cat("\n----------\n\n")
    
    
    # Plot (if requested)
    if(plot){
      results$prop_plot = 
        tibble::tibble(x = seq(0.001,0.999,
                               l = 50)) |> 
        ggplot(aes(x=.data$x)) +
        stat_function(fun = 
                        function(x){
                          dbeta(x,
                                prior_shapes[1],
                                prior_shapes[2])
                        },
                      aes(color = "Prior"),
                      linewidth = 2) + 
        stat_function(fun = 
                        function(x){
                          dbeta(x,
                                post_shapes[1,1],
                                post_shapes[1,2])
                        },
                      aes(color = "Posterior (Pop1)"),
                      linewidth = 2) + 
        stat_function(fun = 
                        function(x){
                          dbeta(x,
                                post_shapes[2,1],
                                post_shapes[2,2])
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
        ggtitle("Population proportion")
      
      print(results$prop_plot)
      
    }
    
    # Add posterior parameters to returned object
    results$posterior_parameters = list()
    results$posterior_parameters$population_1 = 
      c(shape_1 = post_shapes[1,1],
        shape_2 = post_shapes[1,2])
    results$posterior_parameters$population_2 = 
      c(shape_1 = post_shapes[2,1],
        shape_2 = post_shapes[2,2])
    
    
    invisible(results)
  }#End: 2 sample inference
}

