#' @name IC
#' @aliases AIC
#' @aliases BIC
#' @aliases DIC
#' @aliases WAIC
#' 
#' @title Compute AIC, BIC, DIC, or WAIC for aov_b or lm_b objects.  (Lower is better.)  
#' 
#' @param object aov_b, lm_b, or glm_b object
#' @param seed integer.  Always set your seed!!!
#' @param mc_error The number of posterior draws will ensure that 
#' with 99% probability the posterior mean of the deviance for DIC will be 
#' within \eqn{\pm}\code{mc_error}E(deviance).
#' @param ... Passed to methods.
#' 
#' @returns Numeric (or in the case of DIC, a numeric vector)
#' 
#' @examples
#' set.seed(2025)
#' N = 500
#' test_data <-
#'   data.frame(x1 = rnorm(N),
#'              x2 = rnorm(N),
#'              x3 = letters[1:5])
#' test_data$outcome <-
#'   rnorm(N,-1 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) )
#' fit1 <-
#'   lm_b(outcome ~ x1 + x2 + x3,
#'        data = test_data)
#' AIC(fit1)
#' BIC(fit1)
#' DIC(fit1)
#' WAIC(fit1)
#' 
#' 
#' 
#' @export


DIC = function(object, ...){
  UseMethod("DIC")
}

#' @export
WAIC = function(object, ...){
  UseMethod("WAIC")
}

#' @rdname IC
#' @exportS3Method BIC lm_b
BIC.lm_b = function(object, ...){
  y = model.frame(object$formula,
                  object$data)[,1]
  
  llik = 
    dnorm(y,
          mean = object$fitted,
          sd = sqrt(0.5 * object$posterior_parameters$b_tilde / 
                      (0.5 * object$posterior_parameters$a_tilde + 1.0) ),
          log = TRUE) |> 
    sum()
  
  -2.0 * llik + log(nrow(object$data)) * (length(object$posterior_parameters$mu_tilde) + 1.0)
}

#' @rdname IC
#' @exportS3Method BIC glm_b
BIC.glm_b = function(object, ...){
  
  y = model.frame(object$formula,
                  object$data)[,1]
  if(object$family$family == "negbinom"){
    object_dev = exp(object$summary$`Post Mean`[nrow(object$summary)])
  }else{
    object_dev = 1.0
  }
  
  object$family$aic(y = y / object$trials,
                    n = object$trials,
                    mu = object$fitted,
                    wt = rep(1.0, NROW(object$fitted)),
                    dev = object_dev) +
    log(nrow(object$data)) * nrow(object$summary)
}

#' @rdname IC
#' @exportS3Method BIC aov_b
BIC.aov_b = function(object, ...){
  G = length(object$posterior_parameters$mu_g)
  nparms = G + length(object$posterior_parameters$a_g)
  
  if(nparms == G+1){
    llik = 
      dnorm(object$data[[all.vars(object$formula)[1]]],
            mean = object$posterior_parameters$mu_g[as.integer(object$data$group)],
            sd = sqrt(0.5 * object$posterior_parameters$b_g / 
                        (0.5 * object$posterior_parameters$a_g + 1.0)),
            log = TRUE) |> 
      sum()
  }else{
    variances = 
      0.5 * object$posterior_parameters$b_g / 
      (0.5 * object$posterior_parameters$a_g + 1.0)
    llik = 
      dnorm(object$data[[all.vars(object$formula)[1]]],
            mean = object$posterior_parameters$mu_g[as.integer(object$data$group)],
            sd = sqrt(variances[as.integer(object$data$group)]),
            log = TRUE) |> 
      sum()
  }
  
  -2.0 * llik + log(nrow(object$data)) * nparms
}


#' @rdname IC
#' @exportS3Method AIC lm_b
AIC.lm_b = function(object, ...){
  y = model.frame(object$formula,
                  object$data)[,1]
  
  llik = 
    dnorm(y,
          mean = object$fitted,
          sd = sqrt(0.5 * object$posterior_parameters$b_tilde / 
                      (0.5 * object$posterior_parameters$a_tilde + 1.0) ),
          log = TRUE) |> 
    sum()
  
  -2.0 * llik + 2.0 * (length(object$posterior_parameters$mu_tilde) + 1.0)
}

#' @rdname IC
#' @exportS3Method AIC glm_b
AIC.glm_b = function(object, ...){
  
  y = model.frame(object$formula,
                  object$data)[,1]
  
  if(object$family$family == "negbinom"){
    object_dev = exp(object$summary$`Post Mean`[nrow(object$summary)])
  }else{
    object_dev = 1.0
  }
  
  object$family$aic(y = y / object$trials,
                    n = object$trials,
                    mu = object$fitted,
                    wt = rep(1.0, NROW(object$fitted)),
                    dev = object_dev) +
    2.0 * nrow(object$summary)
}

#' @rdname IC
#' @exportS3Method AIC aov_b
AIC.aov_b = function(object, ...){
  G = length(object$posterior_parameters$mu_g)
  nparms = G + length(object$posterior_parameters$a_g)
  
  if(nparms == G+1){
    llik = 
      dnorm(object$data[[all.vars(object$formula)[1]]],
            mean = object$posterior_parameters$mu_g[as.integer(object$data$group)],
            sd = sqrt(0.5 * object$posterior_parameters$b_g / 
                        (0.5 * object$posterior_parameters$a_g + 1.0)),
            log = TRUE) |> 
      sum()
  }else{
    variances = 
      0.5 * object$posterior_parameters$b_g / 
      (0.5 * object$posterior_parameters$a_g + 1.0)
    llik = 
      dnorm(object$data[[all.vars(object$formula)[1]]],
            mean = object$posterior_parameters$mu_g[as.integer(object$data$group)],
            sd = sqrt(variances[as.integer(object$data$group)]),
            log = TRUE) |> 
      sum()
  }
  
  -2.0 * llik + 2 * nparms
}


#' @rdname IC
#' @exportS3Method DIC lm_b 
DIC.lm_b = function(object,
                    seed = 1,
                    mc_error = 0.01,
                    ...){
  set.seed(seed)
  
  y = model.frame(object$formula,
                  object$data)[,1]
  X = model.matrix(object$formula,
                   object$data)
  
  p = nrow(object$summary)
  
  
  V_tilde_eig = eigen(object$posterior_parameters$V_tilde)
  if(ncol(X) > 1){
    Vinv_sqrt = tcrossprod(diag(1 / sqrt(V_tilde_eig$values)),
                           V_tilde_eig$vectors)
  }else{
    Vinv_sqrt = drop(V_tilde_eig$vectors) / sqrt(V_tilde_eig$values)
  }
  
  # Get posterior samples 
  ## Get preliminary draws
  post_draws = 
    matrix(0.0,
           500,
           p + 1,
           dimnames = list(NULL,
                           c(object$summary$Variable,"s2")))
  post_draws[,"s2"] = 
    extraDistr::rinvgamma(500,
                          0.5 * object$posterior_parameters$a_tilde,
                          0.5 * object$posterior_parameters$b_tilde)
  post_draws[,1:p] = 
    matrix(1.0,
           500,1) %*% 
    matrix(object$summary$`Post Mean`,nrow=1) +
    matrix(rnorm(500 * p,
                 sd = sqrt(rep(post_draws[,"s2"],p))),
           500,p) %*% 
    Vinv_sqrt
  
  
  llik = 
    future.apply::future_sapply(1:nrow(object$data),
                                function(i){
                                  dnorm(y[i],
                                        mean = drop(tcrossprod(post_draws[,1:p],X[i,,drop=FALSE])),
                                        sd = sqrt(post_draws[,p + 1]),
                                        log = TRUE)
                                })
  E_D_draws = -2.0 * rowSums(llik)
  n_draws = 
    var(E_D_draws) / 
    (mean(E_D_draws) * mc_error)^2 *
    qnorm(0.5 * (1.0 - 0.99))^2
  ## Get remaining draws if needed.
  if(n_draws > 500){
    post_draws = 
      matrix(0.0,
             n_draws,
             p + 1,
             dimnames = list(NULL,
                             c(object$summary$Variable,"s2")))
    post_draws[,"s2"] = 
      extraDistr::rinvgamma(n_draws,
                            0.5 * object$posterior_parameters$a_tilde,
                            0.5 * object$posterior_parameters$b_tilde)
    post_draws[,1:p] = 
      matrix(1.0,
             n_draws,1) %*% 
      matrix(object$summary$`Post Mean`,nrow=1) +
      matrix(rnorm(n_draws * p,
                   sd = sqrt(rep(post_draws[,"s2"],p))),
             n_draws,p) %*% 
      Vinv_sqrt
    
    
    llik = 
      future.apply::future_sapply(1:nrow(object$data),
                                  function(i){
                                    dnorm(y[i],
                                          mean = drop(tcrossprod(post_draws[,1:p],X[i,,drop=FALSE])),
                                          sd = sqrt(post_draws[,p + 1]),
                                          log = TRUE)
                                  })
  }
  
  E_D = -2 * mean(rowSums(llik))
  
  
  # Finish computing DIC
  D_E = 
    -2.0 * 
    dnorm(y,
          mean = object$fitted,
          sd = sqrt(0.5 * object$posterior_parameters$b_tilde / 
                      (0.5 * object$posterior_parameters$a_tilde + 1.0)),
          log = TRUE) |> 
    sum()
  
  p_D = E_D - D_E
  
  c(DIC = D_E + 2 * p_D,
    eff_n_parms = p_D)
}


#' @rdname IC
#' @exportS3Method DIC glm_b  
DIC.glm_b = function(object,
                     seed = 1,
                     ...){
  set.seed(seed)
  
  mframe = model.frame(object$formula, object$data)
  y = model.response(mframe)
  X = model.matrix(object$formula,object$data)
  os = model.offset(mframe)
  if(is.null(os)) os = numeric(nrow(object$data))
  
  n_draws = 1e4
  if(object$family$family == "negbinom"){
    object_dev = exp(object$summary$`Post Mean`[nrow(object$summary)])
  }else{
    object_dev = 1.0
  }
  
  D_E = 
    object$family$aic(y = y / object$trials,
                      n = object$trials,
                      mu = object$fitted,
                      wt = rep(1.0, NROW(object$fitted)),
                      dev = object_dev)
  
  if("posterior_covariance" %in% names(object)){
    post_draws = 
      mvtnorm::rmvnorm(n_draws,
                       mean = object$summary$`Post Mean`,
                       sigma = object$posterior_covariance)
  }else{#End: large sample approx
    # If IS was used, use SIR
    post_draws = 
      object$proposal_draws[sample(1:NROW(object$importance_sampling_weights),
                                   n_draws,
                                   TRUE,
                                   object$importance_sampling_weights),,drop=FALSE]
  }#End: IS approach
  
  if(object$family$family == "negbinom"){
    post_draws_dev = exp(post_draws[,ncol(X) + 1])
  }else{
    post_draws_dev = 1.0
  }
  
  mu = 
    (os + tcrossprod(X,post_draws[,1:ncol(X)])) |> 
    object$family$linkinv()
  deviance_draws = 
    future.apply::future_sapply(1:n_draws,
                                function(i){
                                  object$family$aic(y = y / object$trials,
                                                    n = object$trials,
                                                    mu = mu[,i],
                                                    wt = rep(1.0, NROW(object$fitted)),
                                                    dev = post_draws_dev[i])
                                })
  E_D = mean(deviance_draws)
  
  p_D = E_D - D_E
  
  c(DIC = D_E + 2 * p_D,
    eff_n_parms = p_D)
}

#' @rdname IC
#' @exportS3Method DIC aov_b 
DIC.aov_b = function(object, ...){
  G = length(object$posterior_parameters$mu_g)
  nparms = G + length(object$posterior_parameters$a_g)
  
  if(nparms == G+1){
    D_E = 
      -2.0 * 
      dnorm(object$data[[all.vars(object$formula)[1]]],
            mean = object$posterior_parameters$mu_g[as.integer(object$data$group)],
            sd = sqrt(0.5 * object$posterior_parameters$b_g / 
                        (0.5 * object$posterior_parameters$a_g + 1.0)),
            log = TRUE) |> 
      sum()
  }else{
    variances = 
      0.5 * object$posterior_parameters$b_g / 
      (0.5 * object$posterior_parameters$a_g + 1.0)
    D_E = 
      -2.0 * 
      dnorm(object$data[[all.vars(object$formula)[1]]],
            mean = object$posterior_parameters$mu_g[as.integer(object$data$group)],
            sd = sqrt(variances[as.integer(object$data$group)]),
            log = TRUE) |> 
      sum()
  }
  
  
  if(nparms == G+1){
    llik = 
      sapply(1:nrow(object$data),
             function(i){
               dnorm(object$data[[all.vars(object$formula)[1]]][i],
                     mean = object$posterior_draws[,as.integer(object$data$group)[i]],
                     sd = sqrt(object$posterior_draws[,G + 1]),
                     log = TRUE)
             })
  }else{
    llik = 
      sapply(1:nrow(object$data),
             function(i){
               dnorm(object$data[[all.vars(object$formula)[1]]][i],
                     mean = object$posterior_draws[,as.integer(object$data$group)[i]],
                     sd = sqrt(object$posterior_draws[,G + as.integer(object$data$group)[i]]),
                     log = TRUE)
             })
  }
  E_D = -2 * mean(rowSums(llik))
  
  p_D = E_D - D_E
  
  c(DIC = D_E + 2 * p_D,
    eff_n_parms = p_D)
}


#' @rdname IC
#' @exportS3Method WAIC lm_b 
WAIC.lm_b = function(object,
                     seed = 1,
                     ...){
  set.seed(seed)
  y = model.frame(object$formula,
                  object$data)[,1]
  X = model.matrix(object$formula,
                   object$data)
  
  n_draws = 1e4
  p = nrow(object$summary)
  n = nrow(X)
  
  V_tilde_eig = eigen(object$posterior_parameters$V_tilde)
  if(ncol(X) > 1){
    Vinv_sqrt = tcrossprod(diag(1 / sqrt(V_tilde_eig$values)),
                           V_tilde_eig$vectors)
  }else{
    Vinv_sqrt = drop(V_tilde_eig$vectors) / sqrt(V_tilde_eig$values)
  }
  post_draws = 
    matrix(0.0,n_draws,p + 1,
           dimnames = list(NULL,
                           c(object$summary$Variable,"s2")))
  post_draws[,"s2"] = 
    extraDistr::rinvgamma(n_draws,
                          0.5 * object$posterior_parameters$a_tilde,
                          0.5 * object$posterior_parameters$b_tilde)
  post_draws[,1:p] = 
    matrix(1.0,n_draws,1) %*% matrix(object$summary$`Post Mean`,nrow=1) +
    matrix(rnorm(n_draws*p,
                 sd = sqrt(rep(post_draws[,"s2"],p))),n_draws,p) %*% Vinv_sqrt
  
  betaX = tcrossprod(post_draws[,1:p],X)
  lik_i = 
    matrix(0.0,n_draws,n)
  for(i in 1:n){
    lik_i[,i] = 
      dnorm(y[i],
            mean = betaX[,i],
            sd = sqrt(post_draws[,"s2"]))
  }
  
  lppd = 
    lik_i |>
    colMeans() |>
    log() |>
    sum()
  
  p_waic2 = 
    lik_i |>
    log() |>
    apply(2,var) |>
    sum()
  
  -2.0 * (
    lppd - p_waic2
  )
}

#' @rdname IC
#' @exportS3Method WAIC aov_b 
WAIC.aov_b = function(object,
                      ...){
  G = length(object$posterior_parameters$mu_g)
  nparms = G + length(object$posterior_parameters$a_g)
  n_draws = nrow(object$posterior_draws)
  n = nrow(object$data)
  
  lik_i = 
    matrix(0.0,n_draws,n)
  
  if(nparms == G+1){
    
    for(i in 1:n){
      lik_i[,i] = 
        dnorm(object$data[[all.vars(object$formula)[1]]][i],
              mean = 
                object$posterior_draws[,paste0("mean_",
                                               object$data$group[i])],
              sd = 
                sqrt(object$posterior_draws[,"Var"])
        )
    }
    
  }else{
    
    for(i in 1:n){
      lik_i[,i] = 
        dnorm(object$data[[all.vars(object$formula)[1]]][i],
              mean = 
                object$posterior_draws[,paste0("mean_",
                                               object$data$group[i])],
              sd = 
                sqrt(
                  object$posterior_draws[,paste0("variance_",
                                                 object$data$group[i])])
        )
    }
    
  }
  
  lppd = 
    lik_i |>
    colMeans() |>
    log() |>
    sum()
  
  p_waic2 = 
    lik_i |>
    log() |>
    apply(2,var) |>
    sum()
  
  -2.0 * (
    lppd - p_waic2
  )
  
}


#' @rdname IC
#' @exportS3Method WAIC glm_b 
WAIC.glm_b = function(object,
                      seed = 1,
                      ...){
  set.seed(seed)
  mframe = model.frame(object$formula, object$data)
  y = model.response(mframe)
  X = model.matrix(object$formula,object$data)
  os = model.offset(mframe)
  if(is.null(os)) os = numeric(nrow(object$data))
  
  n_draws = 1e4
  
  if("posterior_covariance" %in% names(object)){
    post_draws = 
      mvtnorm::rmvnorm(n_draws,
                       mean = object$summary$`Post Mean`,
                       sigma = object$posterior_covariance)
  }else{#End: large sample approx
    # If IS was used, use SIR
    post_draws = 
      object$proposal_draws[sample(1:NROW(object$importance_sampling_weights),
                                   n_draws,
                                   TRUE,
                                   object$importance_sampling_weights),,drop=FALSE]
    
  }#End: IS approach
  
  mu = 
    (os + tcrossprod(X,post_draws[,1:ncol(X)])) |> 
    object$family$linkinv()
  if(object$family$family == "poisson"){
    llik_i = 
      future.apply::future_sapply(1:nrow(object$data),
                                  function(i){
                                    dpois(y[i],
                                          mu[i,],
                                          log = TRUE)
                                  })
  }
  if(object$family$family == "binomial"){
    llik_i = 
      future.apply::future_sapply(1:nrow(object$data),
                                  function(i){
                                    dbinom(y[i],
                                           object$trials[i],
                                           mu[i,],
                                           log = TRUE)
                                  })
  }
  if(object$family$family == "negbinom"){
    phi = exp(post_draws[,ncol(X) + 1])
    llik_i = 
      future.apply::future_sapply(1:nrow(object$data),
                                  function(i){
                                    dnbinom(y[i],
                                            mu = mu[i,],
                                            size = phi[i],
                                            log = TRUE)
                                  })
  }
  
  lppd = 
    llik_i |>
    exp() |> 
    colMeans() |>
    log() |>
    sum()
  
  p_waic2 = 
    llik_i |>
    apply(2,var) |>
    sum()
  
  -2.0 * (
    lppd - p_waic2
  )
}
