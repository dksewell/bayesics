#' Predict method for lm_b model fits
#' 
#' 
#' @param object Object of class lm_b
#' @param newdata An optional data.frame in which to look for variables with which 
#' to predict. 
#' @param trials Integer vector giving the number of trials for each 
#' observation if family = binomial().
#' @param CI_level Posterior probability covered by credible interval
#' @param PI_level Posterior probability covered by prediction interval
#' @param seed integer.  Always set your seed!!!
#' @param ... optional arguments.
#' 
#' @return tibble with estimate (posterior mean), prediction intervals, and credible intervals 
#' for the mean.
#' 
#' @exportS3Method predict glm_b


predict.glm_b = function(object,
                         newdata,
                         trials,
                         CI_level = 0.95,
                         PI_level = 0.95,
                         seed = 1,
                         n_draws = 5e3,
                         ...){
  
  alpha_ci = 1.0 - CI_level
  alpha_pi = 1.0 - PI_level
  
  
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
      if(class(trials) == "character") trials = newdata[[trials]]
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
                                   object$summary$`Post Mean`[1:ncol(X)]) + 
                            os)
  
  
  ## Get CI bounds
  if("posterior_covariance" %in% names(object)){
    
    if(object$family$family %in% c("poisson","negbinom")){
      grad_ginv_xbeta = 
        drop(exp(X %*% object$summary$`Post Mean`[1:ncol(X)] + os)) * X
    }
    if(object$family$family == "binomial"){
      probs = 
        1.0 / (1.0 + drop(exp(-X %*% object$summary$`Post Mean` - os)))
      grad_ginv_xbeta = 
        trials * probs * (1.0 - probs) * X
    }
    
    yhats_covar =
      rowSums(
        grad_ginv_xbeta * 
          (grad_ginv_xbeta %*%
             object$posterior_covariance[1:ncol(X),1:ncol(X)]
          )
      )
    
    
    newdata =
      newdata |> 
      mutate(`Post Mean` = yhats,
             CI_lower = 
               qnorm(alpha / 2.0,
                     yhats,
                     sqrt(yhats_covar)),
             CI_upper = 
               qnorm(1.0 - alpha / 2.0,
                     yhats,
                     sqrt(yhats_covar)))
    
    yhats_sds = 
      sqrt(yhats_covar)
    
    
    if(object$family$family == "poisson"){
      y_draws =
        future.apply::future_sapply(1:n_draws, 
                      function(i){
                        rpois(nrow(newdata),
                              pmax(rnorm(nrow(newdata),
                                         newdata$`Post Mean`,
                                         yhats_sds),
                                   .Machine$double.eps)
                        )
                      },
                      future.seed = seed)
      if( ((n_draws > 1) && (NCOL(y_draws) == 1)) || (length(y_draws) == 1) )
        y_draws = matrix(y_draws,nrow = 1)
      
    }
    if(object$family$family == "binomial"){
      y_draws =
        future.apply::future_sapply(1:n_draws, 
                      function(i){
                        rbinom(nrow(newdata),
                               trials,
                               pmin(
                                 pmax(
                                   rnorm(nrow(newdata),
                                         newdata$`Post Mean`,
                                         yhats_sds),
                                   .Machine$double.eps),
                                 1.0 - .Machine$double.eps)
                               )
                      },
                      future.seed = seed)
      if( ((n_draws > 1) && (NCOL(y_draws) == 1)) || (length(y_draws) == 1) )
        y_draws = matrix(y_draws,nrow = 1)
    }
    if(object$family$family == "negbinom"){
      theta_draws =
        mvtnorm::rmvnorm(n_draws,
                         mean = object$summary$`Post Mean`,
                         sigma = object$posterior_covariance)
      y_draws =
        future.apply::future_sapply(1:n_draws, 
                      function(i){
                        rnbinom(nrow(newdata),
                                mu = pmax(drop(exp(X %*% theta_draws[i,1:ncol(X)] + os)),
                                          .Machine$double.eps),
                                size = exp(theta_draws[i,ncol(X) + 1])
                        )
                      },
                      future.seed = seed)
      if( ((n_draws > 1) && (NCOL(y_draws) == 1)) || (length(y_draws) == 1) )
        y_draws = matrix(y_draws,nrow = 1)
    }
    
    if(n_draws > 1){
      PI_bounds = 
        y_draws |> 
        future.apply::future_apply(1,quantile,probs = c(0.5 * alpha_pi,
                                          1.0 - 0.5 * alpha_pi))
      newdata$PI_lower = 
        PI_bounds[1,]
      newdata$PI_upper = 
        PI_bounds[2,]
    }
    
  }else{#End: if asymptotic approx or VB was used.
    
    yhat_draws = 
      trials * 
      object$family$linkinv(os + tcrossprod(X, object$proposal_draws[,1:ncol(X)]))
    
    CI_from_weighted_sample = function(x,w,level = alpha_ci){
      w = cumsum(w[order(x)])
      x = x[order(x)]
      LB = max(which(w <= 0.5 * level))
      UB = min(which(w >= 1.0 - 0.5 * level))
      return(c(lower = x[LB],
               upper = x[UB]))
    }
    CI_bounds = 
      apply(yhat_draws,1,
            CI_from_weighted_sample,
            w = object$importance_sampling_weights)
    
    newdata =
      newdata |> 
      tibble::as_tibble() |> 
      dplyr::mutate(`Post Mean` = yhats,
             CI_lower = 
               CI_bounds["lower",],
             CI_upper = 
               CI_bounds["upper",])
    
    if(object$family$family == "poisson"){
      y_draws = 
        future.apply::future_sapply(1:ncol(yhat_draws),
                      function(i){
                        rpois(nrow(yhat_draws),yhat_draws[,i])
                      },
                      future.seed = seed)
      if(NCOL(y_draws) == 1)
        y_draws = matrix(y_draws,nrow = 1)
    }
    
    if(object$family$family == "binomial"){
      y_draws = 
        future.apply::future_sapply(1:ncol(yhat_draws),
                      function(i){
                        rbinom(nrow(yhat_draws),
                               trials,
                               yhat_draws[,i] / trials)
                      },
                      future.seed = seed)
      if(NCOL(y_draws) == 1)
        y_draws = matrix(y_draws,nrow = 1)
    }
    
    if(object$family$family == "negbinom"){
      y_draws =
        future.apply::future_sapply(1:ncol(yhat_draws),
                      function(i){
                        rnbinom(nrow(newdata),
                                mu = pmax(drop(exp(X %*% object$proposal_draws[i,1:ncol(X)] + os)),
                                          .Machine$double.eps),
                                size = exp(object$proposal_draws[i,ncol(X) + 1])
                        )
                      },
                      future.seed = seed)
      if(NCOL(y_draws) == 1)
        y_draws = matrix(y_draws,nrow = 1)
      
    }
    
    PI_bounds = NULL
    try({
      PI_bounds = 
        y_draws |> 
        future.apply::future_apply(1,CI_from_weighted_sample,
                                   w = object$importance_sampling_weights,
                                   level = alpha_pi)
    },silent=TRUE)
    if(is.null(PI_bounds)){
      PI_bounds = 
        y_draws |> 
        apply(1,CI_from_weighted_sample,
              w = object$importance_sampling_weights,
              level = alpha_pi)
    }
    newdata$PI_lower = 
      PI_bounds[1,]
    newdata$PI_upper = 
      PI_bounds[2,]
    
  }
  
  colnames(y_draws) = paste("y_new",1:ncol(y_draws),sep="")
  newdata =
    newdata |> 
    dplyr::bind_cols(y_draws |> 
                as_tibble())
  
  return(newdata)
}
