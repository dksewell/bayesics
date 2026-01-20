#' Predict method for lm_b model fits
#' 
#' 
#' @param object Object of class lm_b
#' @param newdata An optional data.frame in which to look for variables with which 
#' to predict. 
#' @param trials Integer vector giving the number of trials for each 
#' observation if family = binomial().
#' @param CI_level numeric. Credible interval level.
#' @param ... optional arguments.
#' 
#' @returns tibble with estimate, prediction intervals, and credible intervals 
#' for the mean.
#' 
#' @examples
#' \donttest{
#' # Generate some data
#' set.seed(2025)
#' N = 500
#' test_data = 
#'   data.frame(x1 = rnorm(N),
#'              x2 = rnorm(N),
#'              x3 = letters[1:5])
#' test_data$outcome = 
#'   rbinom(N,1,1.0 / (1.0 + exp(-(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) ))))
#' 
#' # Fit the GLM via the (non-parametric) loss-likelihood bootstrap.
#' fit1 <-
#'   np_glm_b(outcome ~ x1 + x2 + x3,
#'            data = test_data,
#'            family = binomial())
#' predict(fit1)
#' }
#' 
#' @exportS3Method predict np_glm_b


predict.np_glm_b = function(object,
                            newdata,
                            trials,
                            CI_level = 0.95,
                            ...){
  
  
  if(missing(newdata)){
    newdata = object$data
  }
  
  if(!is.null(object$xlevels)){
    for(j in names(object$xlevels)){
      if(!("factor" %in% class(newdata[[j]]))){
        newdata[[j]] = 
          factor(newdata[[j]],
                 levels = object$xlevels[[j]])
      }
    }
  }
  
  
  # Extract 
  mframe = model.frame(delete.response(terms(object)),
                       data = newdata)
  X = model.matrix(delete.response(terms(object)),
                   data = newdata)
  os = model.offset(mframe)
  N = nrow(X)
  p = ncol(X)
  alpha = 1 - CI_level
  if(is.null(os)) os = numeric(N)
  
  # Get trials variables sorted
  if(object$family$family == "binomial"){
    if(missing(trials)){
      message("Assuming all observations correspond to Bernoulli, i.e., Binomial with one trial.")
      trials = rep(1.0,N)
    }else{
      if(is.character(trials)) trials = newdata[[trials]]
      trials = as.numeric(trials)
    }
  }else{
    trials = rep(1.0,N)
  }
  
  
  # Assign new values
  ## Get estimates
  yhats = 
    trials * 
    object$family$linkinv(eta = 
                            drop(X %*% 
                                   object$summary$`Post Mean`[1:p]) + 
                            os)
  ## Get CI bounds
  if("posterior_covariance" %in% names(object)){
    
    if(object$family$family == "gaussian"){
      grad_ginv_xbeta = 
        X
    }
    if(object$family$family %in% c("poisson","negbinom")){
      grad_ginv_xbeta = 
        drop(exp(X %*% object$summary$`Post Mean`[1:p] + os)) * X
    }
    if(object$family$family == "binomial"){
      probs = 
        1.0 / (1.0 + drop(exp(-X %*% object$summary$`Post Mean` - os)))
      grad_ginv_xbeta = 
        trials * probs * (1.0 - probs) * X
    }
    
    yhats_covar = 
      tcrossprod(grad_ginv_xbeta %*% object$posterior_covariance[1:p,1:p],
                 grad_ginv_xbeta)
    
    newdata =
      newdata |> 
      dplyr::mutate(`Post Mean` = yhats,
                    CI_lower = 
                      qnorm(alpha / 2.0,
                            yhats,
                            sqrt(diag(yhats_covar))),
                    CI_upper = 
                      qnorm(1.0 - alpha / 2.0,
                            yhats,
                            sqrt(diag(yhats_covar))))
    
  }else{
    yhat_draws = 
      trials * 
      object$family$linkinv(os + tcrossprod(X, object$posterior_draws[,1:p]))
    if(NCOL(yhat_draws) == 1)
      yhat_draws = matrix(yhat_draws,nrow = 1)
    
    newdata =
      newdata |> 
      tibble::as_tibble() |> 
      dplyr::mutate(`Post Mean` = yhats,
                    CI_lower = 
                      yhat_draws |> 
                      apply(1,quantile, probs = alpha / 2.0),
                    CI_upper = 
                      yhat_draws |> 
                      apply(1,quantile, probs = 1.0 - alpha / 2.0))
  }
  
  
  return(newdata)
}
