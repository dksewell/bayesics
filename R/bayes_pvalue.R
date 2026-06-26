#' @rdname bayes_pvalue
#' 
#' @title Bayesian P-values for Regression Models
#' 
#' @param object object of class \code{lm_b} or \code{aov_b}
#' @param statistic Statistic used to compute Bayesian p-value.  
#' Either "deviance", or else a function taking in data, expected value, and 
#' if applicable to the family, dispersion (residual variance for 
#' \code{gaussian} and \eqn{\phi} for \code{negbinom}, where 
#' \eqn{Var(y) = \mu + \mu^2/\phi}).
#' @param mc_error The number of posterior draws will ensure that with 
#' 99% probability the estimated Bayesian p-value will be within 
#' \eqn{\pm} \code{mc_error} of the actual Bayesian p-value.
#' @param seed integer.
#' @param ... optional arguments.
#' 
#' @returns named list with:
#' \itemize{
#'  \item bpvalue - numeric between 0 and 1, giving the Bayesian p-value
#'  \item statistic_posterior_draws - tibble with two columns, one for 
#'  \eqn{T(y_{obs},\theta)} and the other for \eqn{T(y_{pred},\theta)}.  
#'  See details.
#' }
#' 
#' @details
#' 
#' \strong{Overview:}
#' 
#' Bayesian p-values are measures of how well the model match the data, as 
#' evaluated through the predictive posterior distribution.  While they can 
#' be extremely flexible - testing very specific aspects of the model being 
#' fitted-, the default setting for \code{bayes_pvalue} is the deviance, 
#' acting to do an overall goodness-of-fit test.  More generally, a
#' Bayesian p-value takes a test statistic of the data and the model 
#' parameters \eqn{T(y,\theta)} and compares the posterior probability that 
#' the test statistic evaluated at the observed data compared to data randomly 
#' generated according to the model.  I.e., the Bayesian p-value is given by
#' \deqn{
#'  \Pr(T(y_{obs},\theta) > T(y_{pred},\theta) | y_{pred}).
#' }
#' 
#' 
#' \strong{MC error:}
#' 
#' The number of posterior samples is determined by the fact that if the true 
#' Bayesian p-value is \eqn{p}, the Monte Carlo estimate of the Bayesian p-value 
#' will have 99% probability of being with \eqn{\approx 2.3 \sqrt{p(1-p)/L}}, 
#' where \eqn{L} is the number of posterior draws.  The worst case scenario, in 
#' terms of MC variance is when the Bayesian p-value is \eqn{p=0.5}.  However, in 
#' such a case, there will be no question that the model fit is adequate and we 
#' can afford a much larger MC error.  It is nearer thresholds (typically  
#' Bayesian p-values less than 0.05 or greater than 0.95 are cause for alarm) 
#' where we need more precise estimates of what the Bayesian p-value actually 
#' is, but at near these boundaries the MC error is much smaller than the worst 
#' case scenario of \eqn{p=0.5}.  Hence we compute the number of posterior 
#' samples to be within the user-specified \code{mc_error} at \eqn{p(1-p)=(0.15)(0.85)}. 
#' 
#' The Monte Carlo error is implemented in such a way so as to obtain as few 
#' posterior samples as necessary to ensure a small probability of incorrectly 
#' determining an inadequate fit when the fit is good and vice versa. 
#' 
#' 
#' @examples
#' \donttest{
#' # Create some data
#' set.seed(2026)
#' N = 500
#' test_data = 
#'   data.frame(x1 = rnorm(N),
#'              x2 = rnorm(N),
#'              x3 = letters[1:5])
#' test_data$outcome = 
#'   rnorm(N,-1 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) )
#' 
#' # Fit a linear regression model
#' fit = 
#'   lm_b(outcome ~ x1 + x2 + x3,
#'        data = test_data)
#' 
#' # Compute the Bayesian p-value based on the deviance
#' bp_deviance = 
#'   bayes_pvalue(fit)
#' bp_deviance$bpvalue # We want a value near 0.5, say in (0.05,0.95).
#' plot(T_y_observed ~ T_y_predicted,
#'      data = bp_deviance$statistic_posterior_draws,
#'      xlab = expression(T(y[pred],theta)),
#'      ylab = expression(T(y[obs],theta)),
#'      pch = 16,
#'      cex = 0.1,
#'      col = gray(0.5,0.25))
#' abline(0,1,
#'        lwd = 2)
#' 
#' # Use a custom test statistic
#' bp_sw = 
#'   bayes_pvalue(fit,
#'                statistic = 
#'                  function(y,mu,dispersion){
#'                    shapiro.test((y - mu)/sqrt(dispersion))$statistic
#'                  }
#'   )
#' bp_sw$bpvalue
#' } 



#' @export
bayes_pvalue = function(object,
                        statistic,
                        mc_error,
                        seed,
                        ...){
  UseMethod("bayes_pvalue")
}


#' @rdname bayes_pvalue
#' @exportS3Method bayes_pvalue lm_b 
bayes_pvalue.lm_b = function(object,
                             statistic,
                             mc_error = 0.005,
                             seed = 1,
                             ...){
  
  # object
  if (!inherits(object, "lm_b"))
    stop("`object` must be an object of class \"lm_b\"",
         call. = FALSE)
  
  # Check to see if a parametric fit is used
  if(object$model_type != "parametric")
    stop("Object should be a parametric fit in order to obtain posterior predicted values.")
  
  # statistic
  if(!missing(statistic) && class(statistic) != "function")
    stop("Is statistic is provided, it must be a function that takes in y, E(y), and if applicable a dispersion parameter, in that order.")
  if(missing(statistic)){
    
    statistic <- function(y, mu, dispersion = NULL) {
      switch(object$family$family,
             gaussian   = shapiro.test((y - mu) / sqrt(dispersion))$statistic,
             binomial   = -2.0 * sum(dbinom(y,object$trials,mu/object$trials,log=T)),
             poisson    = -2.0 * sum(dpois(y,mu,log=T)),
             negbinom   = -2.0 * sum(dnbinom(y,mu = mu,size = dispersion,log=T)),
             stop("Unsupported family")
      )
    }
    
    if(object$family$family == "gaussian"){
      
    }else{
      
    }
  }
  
  
  # mc_error
  if (!is.numeric(mc_error) ||
      length(mc_error) != 1 ||
      mc_error <= 0)
    stop(
      "`mc_error` must be a positive numeric scalar",
      call. = FALSE
    )
  
  # seed
  if (!is.numeric(seed) ||
      length(seed) != 1 ||
      seed %% 1 != 0)
    stop(
      "`seed` must be a single integer value",
      call. = FALSE
    )
  
  
  
  
  
  # Get number of posterior draws required (see details)
  n_draws = 
    ceiling(qnorm(0.99)^2 / mc_error^2 * sqrt(0.15 * 0.85))
  
  # Get posterior draws
  theta_draws = 
    get_posterior_draws(object,
                        n_draws = n_draws,
                        seed = seed)
  
  # Extract 
  mframe = model.frame(object$formula,
                       object$data)
  y = 
    model.response(mframe)
  if(is.character(y))
    y = factor(y)
  if(is.factor(y)){
    y = as.integer(y)
    if(length(unique(y)) == 2) y = y - 1
  }
  X = model.matrix(object$formula,
                   object$data)
  p = ncol(X)
  os = model.offset(mframe)
  N = nrow(X)
  if(is.null(os)) os = numeric(N)
  if(is.null(object$trials)) object$trials = rep(1.0,nrow(object$data))
  if(is.null(object$weights)) # bma_inference does not take in weights currently.
    object$weights = rep(1.0,N)
  
  
  
  # Convert beta to E(y) and get dispersion parameter
  mu_draws = 
    object$trials * 
    object$family$linkinv(os + tcrossprod(X, as.matrix(theta_draws[,1:p])))
  
  if(object$family$family == "gaussian"){
    phi = 
      tcrossprod(1.0 / object$weights,
                 unlist(theta_draws[,p + 1]))
  }else{
    if(object$family$family == "negbinom"){
      phi = 
        matrix(exp(theta_draws[,p + 1]),
               N,n_draws,byrow = TRUE)
    }else{
      phi = NULL
    }
  }
  
  # Draw y_pred
  ## Gaussian
  if(object$family$family == "gaussian"){
    y_pred = 
      sapply(1:n_draws,
             function(draw){
               rnorm(N,
                     mu_draws[,draw],
                     sd = sqrt(phi[,draw]))
             })
  }
  
  ## Binomial
  if(object$family$family == "binomial"){
    y_pred = 
      sapply(1:n_draws,
             function(draw){
               rbinom(N,
                      size = object$trials,
                      prob = mu_draws[,draw] / object$trials)
             })
  }
  
  ## Poisson
  if(object$family$family == "poisson"){
    y_pred = 
      sapply(1:n_draws,
             function(draw){
               rpois(N,
                     lambda = mu_draws[,draw])
             })
  }
  
  ## Negative binomial
  if(object$family$family == "negbinom"){
    y_pred = 
      sapply(1:n_draws,
             function(draw){
               rnbinom(N,
                       mu = mu_draws[,draw],
                       size = phi[,draw])
             })
  }
  
  
  
  # Evaulate T(y,theta) and T(y_pred,theta)
  
  ## Compute posterior draws of test statistic
  if(is.null(phi)){
    T_pred = 
      sapply(1:n_draws,
             function(draw){
               statistic(y_pred[,draw],
                         mu_draws[,draw])
             })
    T_obs = 
      sapply(1:n_draws,
             function(draw){
               statistic(y,
                         mu_draws[,draw])
             })
  }else{
    T_pred = 
      sapply(1:n_draws,
             function(draw){
               statistic(y_pred[,draw],
                         mu_draws[,draw],
                         phi[,draw])
             })
    T_obs = 
      sapply(1:n_draws,
             function(draw){
               statistic(y,
                         mu_draws[,draw],
                         phi[,draw])
             })
  }
  
  # Return bayesian p-value
  return(
    list(bpvalue = 
           mean(T_obs > T_pred),
         statistic_posterior_draws = 
           tibble(T_y_observed = T_obs,
                  T_y_predicted = T_pred))
  )
}







#' @rdname bayes_pvalue
#' @exportS3Method bayes_pvalue aov_b 
bayes_pvalue.aov_b = function(object,
                              statistic,
                              mc_error = 0.005,
                              seed = 1,
                              ...){
  
  # object
  if (!inherits(object, "aov_b"))
    stop("`object` must be an object of class \"aov_b\"",
         call. = FALSE)
  
  # statistic
  if(!missing(statistic) && class(statistic) != "function")
    stop("Is statistic is provided, it must be a function that takes in y, E(y), and if applicable a dispersion parameter, in that order.")
  if(missing(statistic)){
    statistic <- function(y, mu, dispersion = NULL) {
      shapiro.test((y - mu) / sqrt(dispersion))$statistic
    }
  }
  
  # mc_error
  if (!is.numeric(mc_error) ||
      length(mc_error) != 1 ||
      mc_error <= 0)
    stop(
      "`mc_error` must be a positive numeric scalar",
      call. = FALSE
    )
  
  # seed
  if (!is.numeric(seed) ||
      length(seed) != 1 ||
      seed %% 1 != 0)
    stop(
      "`seed` must be a single integer value",
      call. = FALSE
    )
  
  # Get number of posterior draws required (see details)
  n_draws = 
    ceiling(qnorm(0.99)^2 / mc_error^2 * sqrt(0.15 * 0.85))
  
  # Get posterior draws
  theta_draws = 
    get_posterior_draws(object,
                        n_draws = n_draws,
                        seed = seed)
  
  
  # Get new draws of y
  ## Get group assignments
  group_assignment =
    as.integer(object$data$group)
  G = length(object$posterior_parameters$nu_g)
  N = nrow(object$data)
  
  ## Get mu_draws 
  mu_draws = 
    sapply(1:n_draws,
           function(draw) theta_draws[draw,group_assignment])
  
  ## Get variance draws
  heteroscedastic = 
    (length(object$posterior_parameters$a_g) > 1)
  if(heteroscedastic){
    s2_draws = 
      sapply(1:n_draws,
             function(draw) theta_draws[draw,G + group_assignment])
  }else{
    s2_draws = 
      matrix(theta_draws[,G + 1],
             N,n_draws,
             byrow = TRUE)
  }
  
  
  ## Draw y_pred
  y_pred = 
    sapply(1:n_draws,
           function(draw){
             rnorm(N,
                   mu_draws[,draw],
                   sd = sqrt(s2_draws[,draw]))
           })
  
  
  
  # Evaulate T(y,theta) and T(y_pred,theta)
  
  ## Compute posterior draws of test statistic
  T_pred = 
    sapply(1:n_draws,
           function(draw){
             statistic(y_pred[,draw],
                       mu_draws[,draw],
                       s2_draws[,draw])
           })
  T_obs = 
    sapply(1:n_draws,
           function(draw){
             statistic(object$data[[all.vars(object$formula)[1]]],
                       mu_draws[,draw],
                       s2_draws[,draw])
           })
  
  # Return bayesian p-value
  return(
    list(bpvalue = 
           mean(T_obs > T_pred),
         statistic_posterior_draws = 
           tibble(T_y_observed = T_obs,
                  T_y_predicted = T_pred))
  )
}

























