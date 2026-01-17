#' Paired sign test
#' 
#' Sign test for paired data.
#' 
#' @details
#' The sign test looks at \eqn{z_i:= 1_{[x_i > y_i]}} rather than trying to model the 
#' distribution of \eqn{(x_i,y_i)}.  \code{sign_test_b} then uses the fact that 
#' \deqn{
#'  z_i \overset{iid}{\sim} Bernoulli(p)
#' }
#' and then makes inference on \eqn{p}.  The prior on \eqn{p} is
#' \deqn{
#'  p \sim Beta(a,b),
#' }
#' where \eqn{a} and \eqn{b} are given by the argument \code{prior_shapes}.  If 
#' \code{prior_shapes} is missing and \code{prior = "jeffreys"}, then a 
#' Jeffreys prior will be used (\eqn{Beta(1/2,1/2)}), and if 
#' \code{prior = "uniform"}, then a uniform prior will be used (\eqn{Beta(1,1)}).
#' 
#' @param x Either numeric vector or binary vector.  If the former, 
#' \eqn{z_i = 1_{[x_i > y_i]}} if \code{y} is supplied, else 
#' \eqn{z_i = 1_{[x_i > 0]}}.  If the latter, then \eqn{z_i = x_i}.
#' @param y Optional numeric vector to pair with \code{x}.
#' @param p0 \code{sign_test_b} will return the posterior probability that 
#' \eqn{p < p0}.  Defaults to 0.5, as is most typical in the sign test.
#' @param prior Either "jeffreys" (Beta(1/2,1/2)) or "uniform" (Beta(1,1)).  
#' This is ignored if prior_shapes is provided.
#' @param prior_shapes Vector of length two, giving the shape parameters 
#' for the beta distribution that will act as the prior on the probability 
#' that \eqn{z_i = 1}.
#' @param ROPE positive numeric of length 1 or 2.  If of length 1, then ROPE 
#' is taken to be \eqn{p0\pm }ROPE. Defaults to \eqn{\pm 0.05}.
#' @param CI_level The posterior probability to be contained in the 
#' credible interval for \eqn{p}.
#' @param plot logical.  Should a plot be shown?
#' 
#' @returns (returned invisible) A list with the following:
#' \itemize{
#'  \item posterior mean of the median difference
#'  \item Credible interval for the median difference
#'  \item ROPE bounds
#'  \item Probability of falling in the ROPE
#'  \item Prior and posterior plot
#'  \item posterior parameters
#' }
#' 
#' @examples
#' # Single population
#' sign_test_b(x = rnorm(50))
#' 
#' ## Change ROPE
#' sign_test_b(x = rnorm(50),
#'             ROPE = 0.1)
#' 
#' ## Change the prior
#' sign_test_b(x = rnorm(50),
#'             prior = "uniform")
#' sign_test_b(x = rnorm(50),
#'             prior_shapes = c(2,2))
#' 
#' ## Two populations
#' sign_test_b(x = rnorm(50,1),
#'             y = rnorm(50,0))
#' 
#' ## Change reference probability
#' sign_test_b(x = rnorm(50,1),
#'             y = rnorm(50,0),
#'             p0 = 0.7)
#' 
#' @export



sign_test_b = function(x,
                       y,
                       p0 = 0.5,
                       prior = "jeffreys",
                       prior_shapes,
                       ROPE,
                       CI_level = 0.95,
                       plot = TRUE){
  
  alpha_ci = 1.0 - CI_level
  
  if( (p0 <= 0) | (p0 >= 1) )
    stop("p0 must be between 0 and 1")
  
  # Get ROPE
  if(missing(ROPE)){
    ROPE = 
      min(c(0.05,
            p0 / 2,
            (1.0 - p0) / 2))
    ROPE_bounds = p0 + c(-1,1) * ROPE
  }else{
    if(!(length(ROPE) %in% 1:2))
      stop("ROPE must be of length 1 or 2")
    if(length(ROPE) == 1){
      if(ROPE < 0)
        stop("ROPE must be positive")
      if(ROPE > min(c(p0,1.0 - p0) / 2))
        stop("ROPE extends from p0 beyond 0 or 1")
      ROPE_bounds = p0 + c(-1,1) * ROPE
    }else{
      if(ROPE[1] < 0)
        stop("Lower bound of ROPE cannot be less than 0")
      if(ROPE[2] < 0)
        stop("Upper bound of ROPE cannot be greater than 1")
      ROPE_bounds = ROPE
    }
  }
  
  # Get z
  if(!missing(y)){
    if(length(x) != length(y))
      stop("x and y must be of the same length.")
    
    z = (x > y) * 1
    plot_title = "Pr(x > y)"
    
  }else{
    if( (length(unique(x)) == 2) &&
        all(near(sort(unique(x)),0:1)) ){
      z = x
      plot_title = "Pr(x = 1)"
    }else{
      z = (x > 0) * 1
      plot_title = "Pr(x > 0)"
    }
  }
  
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
  
  
  # Get posterior parameters
  post_shapes = 
    prior_shapes + 
    c(sum(near(z,1.0)),
      length(z) - sum(near(z,1.0)))
  
  # Compute results
  results = 
    list(posterior_mean = 
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
         Pr_less_than_p = 
           pbeta(p0,post_shapes[1],post_shapes[2]),
         ROPE = 
           pbeta(ROPE_bounds[2],post_shapes[1],post_shapes[2]) - 
           pbeta(ROPE_bounds[1],post_shapes[1],post_shapes[2]),
         ROPE_bounds = ROPE_bounds
    )
  
  
  
  # Print results
  cat("\n----------\n\n Non-parametric sign test using Bayesian techniques\n")
  cat("\n----------\n\n")
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
  cat(paste0("Probability that p < ",
             format(signif(p0, 3), 
                    scientific = FALSE),
             ": ",
             format(signif(results$Pr_less_than_p, 3), 
                    scientific = FALSE),
             "\n\n"))
  cat(paste0("Probability that ",
             format(signif(ROPE_bounds[1], 3), 
                    scientific = FALSE),
             " < p < ",
             format(signif(ROPE_bounds[2], 3), 
                    scientific = FALSE),
             ": ",
             format(signif(results$ROPE, 3), 
                    scientific = FALSE),
             "\n\n"))
  cat("\n----------\n\n")
  
  
  # Plot (if requested)
  if(plot){
    
    results$prop_plot = 
      tibble::tibble(x = seq(0.001,0.999,#seq(.Machine$double.eps,1.0 - .Machine$double.eps,
                             l = 50)) |> 
      ggplot(aes(x=x)) +
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
      ggtitle(plot_title)
    
    print(results$prop_plot)
    
  }
  
  # Add posterior parameters to returned object
  results$posterior_parameters = 
    c(shape_1 = post_shapes[1],
      shape_2 = post_shapes[2])
  
  invisible(results)
}