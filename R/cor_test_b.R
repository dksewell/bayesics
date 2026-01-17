#' Test for Association/Correlation Between Paired Samples
#' @name cor_test_b
#' @rdname cor_test_b
#' @export

#' @export
cor_test_b = function(x, ...){
  UseMethod("cor_test_b")
}

#' @details
#' cor_test_b relies on the robust Kendall's tau, defined to be
#' \deqn{
#'  \tau := \frac{(\# \text{concordant pairs}) - (\# \text{discordant pairs})}{(\# \text{concordant pairs}) - (\# \text{discordant pairs})},
#' }
#' where a concordant pair is a pair of points such that if the rank of the x 
#' values is higher for the first (second) point of the pair, so too the rank 
#' of the y value is higher for the first (second) point of the pair.
#' 
#' The Bayesian approach of Chechile (2020) puts a Beta prior on \eqn{phi}, the 
#' proportion of concordance, i.e., 
#' \deqn{
#'  \phi := \frac{(\# \text{concordant pairs})}{(\# \text{concordant pairs}) - (\# \text{discordant pairs})}.
#' }
#' The relationship between the two, then, is \eqn{\tau = 2\phi - 1}, or 
#' equivalently \eqn{\phi = (\tau + 1)/2}. 
#' 
#' For more information, see \link[DFBA]{dfba_bivariate_concordance} 
#' and vignette("dfba_bivariate_concordance",package = "DFBA").
#' 
#' @param x,y numeric vectors of data values. x and y must have the same length.
#' @param tau  If provided, cor_test_b will return the posterior probability that 
#' Kendall's tau is less than this value.
#' @param ROPE  If a single number, ROPE will be \eqn{\tau\pm} ROPE. If a vector 
#' of length 2, these will serve as the ROPE bounds. Defaults to \eqn{\pm}0.05.
#' @param prior Beta prior used on \eqn{phi} (see details).  Either "uniform' 
#' (Beta(1,1)), "centered' (Beta(2,2)), "positive" (Beta(3.9,2), putting 80% of 
#' the prior mass above 0.5), or "negative" (Beta(2,3.9), putting 80% of the 
#' prior mass below 0.5).
#' @param prior_shapes Vector of length two, giving the shape parameters for 
#' the beta distribution that will act as the prior on \eqn{\phi} (see details).
#' @param CI_level The posterior probability to be contained in the credible 
#' interval.
#' @param plot logical. Should a plot be shown?
#' @param ... optional arguments.
#' 
#' @returns (returned invisible) A list with the following:
#' \itemize{
#'  \item posterior mean
#'  \item lower and upper credible interval bounds
#'  \item Posterior probability that Kendall's tau is less than provided reference tau
#'  \item Posterior probability that Kendall's tau is in the ROPE
#'  \item The underlying object from the DFBA package
#'  \item Posterior and prior plot
#'  \item posterior parameters for Kendall's tau correlations
#' }
#' 
#' @references 
#' Chechile, R.A. (2020). Bayesian Statistics for Experimental Scientists: A General Introduction Using Distribution_Free Statistics. Cambridge: MIT Press.
#' 
#' Chechile, R.A., & Barch, D.H. (2021). A distribution-free, Bayesian goodness-of-fit method for assessing similar scientific prediction equations. Journal of Mathematical Psychology. https://doi.org/10.1016/j.jmp.2021.102638
#' 
#' Lindley, D. V., & Phillips, L. D. (1976). Inference for a Bernoulli process (a Bayesian view). The American Statistician, 30, 112-119.
#' 
#' 
#' @examples
#' # Generate data
#' set.seed(2025)
#' N = 50
#' x = rnorm(N)
#' y = x + 4 * rnorm(N)
#' 
#' # Test for non-zero correlation
#' cor_test_b(x,y)
#' 
#' # Input can be in the form of formula and data
#' cor_test_b(~ asdf + qwer,
#'            data = data.frame(asdf = x,
#'                              qwer = y))
#' 
#' # Other priors can be used, also.  See help for details.
#' cor_test_b(x,y,
#'            prior = "uniform")
#' cor_test_b(x,y,
#'            prior = "negative")
#' cor_test_b(x,y,
#'            prior = "positive")
#' cor_test_b(x,y,
#'            prior_shapes = c(10,10))
#' 
#' 
#' 
#' @export
#' @rdname cor_test_b
cor_test_b.default = function(x,
                              y,
                              tau = 0.0,
                              ROPE,
                              prior = "centered",
                              prior_shapes,
                              CI_level = 0.95,
                              plot = TRUE,
                              ...){
  
  if(length(x) != length(y))
    stop("x and y must be of the same length")
  if( !("numeric" %in% class(x)) | 
      !("numeric" %in% class(y)) )
    stop("x and y must be numeric vectors")
  
  
  ## Get ROPE
  if(missing(ROPE)){
    ROPE = 
      min(
        c(0.05,
          abs(1.0 - tau) * 0.5,
          abs(1.0 + tau) * 0.5)
      )
    
    ROPE_bounds = tau + c(-1,1) * ROPE
  }else{
    if(!(length(ROPE) %in% 1:2))
      stop("ROPE must be of length 1 or 2")
    if(length(ROPE) == 1){
      if(ROPE < 0)
        stop("ROPE must be positive")
      if( (tau - ROPE < -1) | (tau + ROPE > 1) )
        stop("ROPE extends from p beyond 0 or 1")
      ROPE_bounds = tau + c(-1,1) * ROPE
    }else{
      if(diff(ROPE) <= 0)
        stop("ROPE lower bound must be less than the upper bound")
      if(ROPE[1] < -1)
        stop("Lower bound of ROPE cannot be less than -1")
      if(ROPE[2] > 1)
        stop("Upper bound of ROPE cannot be greater than 1")
      ROPE_bounds = ROPE
    }
  }
  ROPE_bounds_transformed = 
    0.5 * (ROPE_bounds + 1.0)
  
  ## Prior distribution
  if(missing(prior_shapes)){
    prior = c("uniform",
              "centered",
              "positive",
              "negative")[pmatch(tolower(prior),
                                 c("uniform",
                                   "centered",
                                   "positive",
                                   "negative"))]
    
    if(prior == "uniform"){
      message("Prior shape parameters were not supplied.\nA uniform prior will be used.")
      prior_shapes = rep(1.0,2)
    }
    if(prior == "centered"){
      message("Prior shape parameters were not supplied.\nBeta(2,2) prior will be used.")
      prior_shapes = rep(2.0,2)
    }
    if(prior == "negative"){
      message("Prior shape parameters were not supplied.\nBeta(2,3.9) prior will be used.")
      prior_shapes = c(2.0,3.9)
    }
    if(prior == "positive"){
      message("Prior shape parameters were not supplied.\nBeta(3.9,2) prior will be used.")
      prior_shapes = c(3.9,2.0)
    }
  }else{
    if(any(prior_shapes <= 0))
      stop("Prior shape parameters must be positive.")
  }
  
  
  # Use DFBA package to do analysis
  dfba_object = 
    DFBA::dfba_bivariate_concordance(x,y,
                                     a0 = prior_shapes[1],
                                     b0 = prior_shapes[2],
                                     prob_interval = CI_level)
  
  # Collate results
  results = 
    list(posterior_mean = 
           2.0 * 
           dfba_object$a_post / 
           sum( dfba_object$a_post + dfba_object$b_post) - 
           1.0,
         CI = 
           2.0 * c(dfba_object$eti_lower ,
                   dfba_object$eti_upper) - 
           1.0,
         Pr_less_than_tau = 
           pbeta(0.5 * (tau + 1.0),
                 dfba_object$a_post,
                 dfba_object$b_post),
         Pr_in_ROPE = 
           pbeta(ROPE_bounds_transformed[2],
                 dfba_object$a_post,
                 dfba_object$b_post) - 
           pbeta(ROPE_bounds_transformed[1],
                 dfba_object$a_post,
                 dfba_object$b_post),
         dfba_bivariate_concordance_object = dfba_object
    )
  
  
  # Plot (if requested)
  if(plot){
    
    results$prob_plot = 
      tibble::tibble(x = seq(-0.999,0.999,
                             l = 50)) |> 
      ggplot(aes(x=x)) +
      stat_function(fun = 
                      function(x){
                        dbeta(0.5 * (x + 1.0),
                              prior_shapes[1],
                              prior_shapes[2])
                      },
                    aes(color = "Prior"),
                    linewidth = 2) + 
      
      stat_function(fun = 
                      function(x){
                        dbeta(0.5 * (x + 1.0),
                              dfba_object$a_post,
                              dfba_object$b_post)
                      },
                    aes(color = "Posterior"),
                    linewidth = 2) + 
      
      scale_color_manual(values = c("Prior" = "#440154FF", 
                                    "Posterior" = "#FDE725FF")) +
      theme_classic(base_size = 15) +
      xlab("") + 
      ylab("") + 
      labs(color = "Distribution") + 
      ggtitle("Kendall's tau")
    
    
    print(results$prob_plot)
    
  }
  
  # Add posterior parameters to returned object
  results$posterior_parameters = 
    c(shape_1 = dfba_object$a_post,
      shape_2 = dfba_object$b_post)
  
  
  # Print results
  cat("\n----------\n\nAnalysis of Kendall's tau using Bayesian techniques\n")
  cat("\n----------\n\n")
  cat(paste0("Prior used (phi): Beta(", 
             format(signif(prior_shapes[1], 3), 
                    scientific = FALSE),
             ",",
             format(signif(prior_shapes[2], 3), 
                    scientific = FALSE),
             ")\n\n"))
  cat(paste0("Posterior mean (tau): ", 
             format(signif(results$posterior_mean, 3), 
                    scientific = FALSE),
             "\n\n"))
  cat(paste0(100 * CI_level,
             "% credible interval (tau): (", 
             format(signif(results$CI[1], 3), 
                    scientific = FALSE),
             ", ",
             format(signif(results$CI[2], 3), 
                    scientific = FALSE),
             ")\n\n"))
  cat(paste0("Probability that tau < ",
             format(signif(tau, 3), 
                    scientific = FALSE),
             " = ",
             format(signif(results$Pr_less_than_tau, 3), 
                    scientific = FALSE),
             "\n\n"))
  cat(paste0("Probability that tau is in the ROPE, defined to be (",
             format(signif(ROPE_bounds[1], 3), 
                    scientific = FALSE),
             ",",
             format(signif(ROPE_bounds[2], 3), 
                    scientific = FALSE),
             ") = ",
             format(signif(results$Pr_in_ROPE, 3), 
                    scientific = FALSE),
             "\n\n")) 
  cat("\n----------\n\n")
  
  
  
  
  invisible(results)
  
}


#' @param formula ADD description!
#' @param data ADD description!
#' @rdname cor_test_b
#' @export
cor_test_b.formula = function(formula,
                              data,
                              tau = 0.0,
                              ROPE,
                              prior = "centered",
                              prior_shapes,
                              CI_level = 0.95,
                              plot = TRUE,
                              ...){
  
  v_names = all.vars(formula)
  x = data[[v_names[1]]]
  y = data[[v_names[2]]]
  
  if(length(x) != length(y))
    stop("x and y must be of the same length")
  if( !("numeric" %in% class(x)) | 
      !("numeric" %in% class(y)) )
    stop("x and y must be numeric vectors")
  
  
  ## Get ROPE
  if(missing(ROPE)){
    ROPE = 
      min(
        c(0.05,
          abs(1.0 - tau) * 0.5,
          abs(1.0 + tau) * 0.5)
      )
    
    ROPE_bounds = tau + c(-1,1) * ROPE
  }else{
    if(!(length(ROPE) %in% 1:2))
      stop("ROPE must be of length 1 or 2")
    if(length(ROPE) == 1){
      if(ROPE < 0)
        stop("ROPE must be positive")
      if( (tau - ROPE < -1) | (tau + ROPE > 1) )
        stop("ROPE extends from p beyond 0 or 1")
      ROPE_bounds = tau + c(-1,1) * ROPE
    }else{
      if(diff(ROPE) <= 0)
        stop("ROPE lower bound must be less than the upper bound")
      if(ROPE[1] < -1)
        stop("Lower bound of ROPE cannot be less than -1")
      if(ROPE[2] > 1)
        stop("Upper bound of ROPE cannot be greater than 1")
      ROPE_bounds = ROPE
    }
  }
  ROPE_bounds_transformed = 
    0.5 * (ROPE_bounds + 1.0)
  
  ## Prior distribution
  if(missing(prior_shapes)){
    prior = c("uniform",
              "centered",
              "positive",
              "negative")[pmatch(tolower(prior),
                                 c("uniform",
                                   "centered",
                                   "positive",
                                   "negative"))]
    
    if(prior == "uniform"){
      message("Prior shape parameters were not supplied.\nA uniform prior will be used.")
      prior_shapes = rep(1.0,2)
    }
    if(prior == "centered"){
      message("Prior shape parameters were not supplied.\nBeta(2,2) prior will be used.")
      prior_shapes = rep(2.0,2)
    }
    if(prior == "negative"){
      message("Prior shape parameters were not supplied.\nBeta(2,3.9) prior will be used.")
      prior_shapes = c(2.0,3.9)
    }
    if(prior == "positive"){
      message("Prior shape parameters were not supplied.\nBeta(3.9,2) prior will be used.")
      prior_shapes = c(3.9,2.0)
    }
  }else{
    if(any(prior_shapes <= 0))
      stop("Prior shape parameters must be positive.")
  }
  
  
  # Use DFBA package to do analysis
  dfba_object = 
    DFBA::dfba_bivariate_concordance(x,y,
                                     a0 = prior_shapes[1],
                                     b0 = prior_shapes[2],
                                     prob_interval = CI_level)
  
  # Collate results
  results = 
    list(posterior_mean = 
           2.0 * 
           dfba_object$a_post / 
           sum( dfba_object$a_post + dfba_object$b_post) - 
           1.0,
         CI = 
           2.0 * c(dfba_object$eti_lower ,
                   dfba_object$eti_upper) - 
           1.0,
         Pr_less_than_tau = 
           pbeta(0.5 * (tau + 1.0),
                 dfba_object$a_post,
                 dfba_object$b_post),
         Pr_in_ROPE = 
           pbeta(ROPE_bounds_transformed[2],
                 dfba_object$a_post,
                 dfba_object$b_post) - 
           pbeta(ROPE_bounds_transformed[1],
                 dfba_object$a_post,
                 dfba_object$b_post),
         dfba_bivariate_concordance_object = dfba_object
    )
  
  
  # Plot (if requested)
  if(plot){
    
    results$prob_plot = 
      tibble(x = seq(-0.999,0.999,
                     l = 50)) |> 
      ggplot(aes(x=x)) +
      stat_function(fun = 
                      function(x){
                        dbeta(0.5 * (x + 1.0),
                              prior_shapes[1],
                              prior_shapes[2])
                      },
                    aes(color = "Prior"),
                    linewidth = 2) + 
      
      stat_function(fun = 
                      function(x){
                        dbeta(0.5 * (x + 1.0),
                              dfba_object$a_post,
                              dfba_object$b_post)
                      },
                    aes(color = "Posterior"),
                    linewidth = 2) + 
      
      scale_color_manual(values = c("Prior" = "#440154FF", 
                                    "Posterior" = "#FDE725FF")) +
      theme_classic(base_size = 15) +
      xlab("") + 
      ylab("") + 
      labs(color = "Distribution") + 
      ggtitle("Kendall's tau")
    
    
    print(results$prob_plot)
    
  }
  
  # Add posterior parameters to returned object
  results$posterior_parameters = 
    c(shape_1 = dfba_object$a_post,
      shape_2 = dfba_object$b_post)
  
  
  # Print results
  cat("\n----------\n\nAnalysis of Kendall's tau using Bayesian techniques\n")
  cat("\n----------\n\n")
  cat(paste0("Prior used (phi): Beta(", 
             format(signif(prior_shapes[1], 3), 
                    scientific = FALSE),
             ",",
             format(signif(prior_shapes[2], 3), 
                    scientific = FALSE),
             ")\n\n"))
  cat(paste0("Posterior mean (tau): ", 
             format(signif(results$posterior_mean, 3), 
                    scientific = FALSE),
             "\n\n"))
  cat(paste0(100 * CI_level,
             "% credible interval (tau): (", 
             format(signif(results$CI[1], 3), 
                    scientific = FALSE),
             ", ",
             format(signif(results$CI[2], 3), 
                    scientific = FALSE),
             ")\n\n"))
  cat(paste0("Probability that tau < ",
             format(signif(tau, 3), 
                    scientific = FALSE),
             " = ",
             format(signif(results$Pr_less_than_tau, 3), 
                    scientific = FALSE),
             "\n\n"))
  cat(paste0("Probability that tau is in the ROPE, defined to be (",
             format(signif(ROPE_bounds[1], 3), 
                    scientific = FALSE),
             ",",
             format(signif(ROPE_bounds[2], 3), 
                    scientific = FALSE),
             ") = ",
             format(signif(results$Pr_in_ROPE, 3), 
                    scientific = FALSE),
             "\n\n")) 
  cat("\n----------\n\n")
  
  
  
  
  invisible(results)
  
}
