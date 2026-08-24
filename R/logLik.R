#' Extract Log-Likelihood
#'
#' Computes the log-likelihood for fitted model objects.
#'
#' @param object An object of class \code{aov_b}, \code{lm_b},
#'   \code{glm_b}, or \code{lm_b_bma}.
# #' @param parameter_estimates If you want to use evaluate the 
# #' log likelihood at a value other than the posterior mean, then 
# #' use \code{parameter_estimates}.  The regression coefficient 
# #' estimates should go first, followed by the auxiliary parameter 
# #' (\eqn{\sigma^2}, \eqn{\phi}) if applicable.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return
#' An object of class \code{"logLik"} with attributes \code{"df"} and
#' \code{"nobs"}.
#'
#' @seealso \code{\link[stats]{logLik}}
#'

#' @rdname logLik
#' @exportS3Method logLik lm_b
logLik.lm_b <- function(object, ...){
  
  if(object$model_type == "nonparametric")
    stop("Cannot compute likelihood for a non-parametric fit.")
  
  if(object$family$family == "gaussian"){
    
    if(!is.numeric(object$standardized_residuals))
      stop("Must compute residuals to compute log likelihood.")
    
    val = 
      sum(dnorm(object$standardized_residuals,
                log = TRUE))
    
  }else{
    
    # Get llik function
    log_lik_function <- function(y, mu, phi = NULL) {
      switch(object$family$family,
             binomial   = dbinom(y,
                                 object$trials,
                                 mu/object$trials,
                                 log = T),
             poisson    = dpois(y,
                                mu,
                                log = T),
             negbinom   = dnbinom(y,
                                  mu = mu,
                                  size = phi,
                                  log = T),
             stop("Unsupported family")
      )
    }
    
    
    # Get data elements for computing mean
    mframe = 
      model.frame(as.formula(paste0(all.vars(object$formula)[1],
                                    "~ 1")),
                  object$data)
    y = model.response(mframe)
    if(is.character(y)) 
      y = factor(y)
    if(is.factor(y)){
      y = as.integer(y)
      if(length(unique(y)) == 2) y = y - 1
    }
    
    phi = 
      ifelse(
        object$family$family == "gaussian",
        object$sigma_sq["Estimate"],
        ifelse(
          object$family$family == "negbinom",
          exp(object$summary$`Post Mean`[nrow(object$summary)]),
          1.0)
      )
    
    val = 
      sum(
        log_lik_function(y,
                         object$fitted,
                         phi)
      )
  }
  
  attr(val,"nall") = nrow(object$data)
  
  attr(val,"nobs") = nrow(object$data)
  s = 
    capture.output(out <- summary(object,
                                  interpretable_scale = FALSE))
  attr(val,"df") = 
    ifelse("list" %in% class(out)[1],
           nrow(out[[1]]),
           nrow(out))
  
  return(structure(val,
                   class = "logLik"))
}