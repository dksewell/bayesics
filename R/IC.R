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
#' \donttest{
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
#' }
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
  
  ll = logLik(object)
  
  -2.0 * ll + 
    log(nrow(object$data)) * attr(ll,"df")
}


#' @rdname IC
#' @exportS3Method AIC lm_b
AIC.lm_b = function(object, ...){
  
  ll = logLik(object)
  
  -2.0 * ll + 
    2.0 * attr(ll,"df")
}


#' @rdname IC
#' @exportS3Method DIC lm_b 
DIC.lm_b = function(object,
                    seed = 1,
                    mc_error = 0.5,
                    ...){
  
  if(object$model_type == "nonparametric")
    stop("Cannot compute likelihood for a non-parametric fit.")
  
  set.seed(seed)
  
  # Get log likelihood function
  log_lik_function <- function(y, mu, phi = NULL) {
    switch(object$family$family,
           gaussian   = dnorm(y,
                              mu,
                              sqrt(phi),
                              log = TRUE),
           binomial   = dbinom(y,
                               object$trials,
                               mu/object$trials,
                               log = TRUE),
           poisson    = dpois(y,
                              mu,
                              log = TRUE),
           negbinom   = dnbinom(y,
                                mu = mu,
                                size = phi,
                                log = TRUE),
           stop("Unsupported family")
    )
  }
  
  # Extract 
  mframe = model.frame(terms(object),
                       data = object$data)
  
  X = model.matrix(delete.response(terms(object)),
                   data = object$data)
  
  os = model.offset(mframe)
  N = nrow(X)
  p = ncol(X)
  if(is.null(os)) os = numeric(N)
  
  y = model.response(mframe)
  if(is.character(y)) 
    y = factor(y)
  if(is.factor(y)){
    y = as.integer(y)
    if(length(unique(y)) == 2) y = y - 1
  }
  
  if("trials" %in% names(object)){
    trials = object$trials
  }else{
    trials = rep(1.0,N)
  }
  
  
  # Get posterior samples 
  ## Get preliminary draws
  
  ### Get draws
  theta_draws = 
    get_posterior_draws(object,
                        n_draws = 500)
  ### Compute draws of linear predictor
  Xbeta_draws = 
    tcrossprod(X,as.matrix(theta_draws[,1:p]))
  ### Compute draws of phi
  if(ncol(theta_draws) > p){
    phi = as.vector(unlist(theta_draws[,p]))
    if(object$family$family == "negbinom")
      phi = exp(phi)
  }else{
    phi = rep(1.0,500)
  }
  ### Compute deviance
  deviance_draws = 
    -2.0 * 
    future.apply::future_sapply(1:nrow(theta_draws),
                                function(i){
                                  log_lik_function(y,
                                                   trials *
                                                     object$family$linkinv(eta = 
                                                                             drop(Xbeta_draws[,i]) + 
                                                                             os),
                                                   phi[i]) |> 
                                    sum()
                                }) |> 
    na.omit()
  
  E_D = mean(deviance_draws)
  n_draws = 
    var(deviance_draws) / 
    (mc_error)^2 *
    qnorm(0.5 * (1.0 - 0.99))^2
  n_draws = round(n_draws)
  n_prelim_draws = length(deviance_draws)
  
  ## Get remaining draws if needed.
  if(n_draws > n_prelim_draws){
    ### Get draws
    theta_draws = 
      rbind(theta_draws,
            get_posterior_draws(object,
                                n_draws = n_draws - n_prelim_draws)
      )
    ### Compute draws of linear predictor
    Xbeta_draws = 
      cbind(
        Xbeta_draws,
        tcrossprod(X,as.matrix(theta_draws[-c(1:n_prelim_draws),1:p]))
      )
    ### Compute draws of phi
    if(ncol(theta_draws) > p){
      phi = as.vector(unlist(theta_draws[,p]))
      if(object$family$family == "negbinom")
        phi = exp(phi)
    }else{
      phi = rep(1.0,n_draws)
    }
    ### Compute deviance
    deviance_draws = 
      c(deviance_draws,
        -2.0 * 
          future.apply::future_sapply((n_prelim_draws + 1):n_draws,
                                      function(i){
                                        log_lik_function(y,
                                                         trials *
                                                           object$family$linkinv(eta = 
                                                                                   drop(Xbeta_draws[,i]) + 
                                                                                   os),
                                                         phi[i]) |> 
                                          sum()
                                      }) |> 
          na.omit()
      )
    
  }
  
  
  # Finish computing DIC
  E_D = mean(deviance_draws)
  
  D_E = 
    -2.0 * 
    as.numeric(logLik(object))
  
  p_D = E_D - D_E
  
  c(DIC = D_E + 2 * p_D,
    eff_n_parms = p_D)
}

#' @rdname IC
#' @exportS3Method DIC aov_b 
DIC.aov_b = function(object, ...){
  
  if(object$model_type == "nonparametric")
    stop("Cannot compute likelihood for a non-parametric fit.")
  
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
  
  if(object$model_type == "nonparametric")
    stop("Cannot compute likelihood for a non-parametric fit.")
  
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
  
  if(object$model_type == "nonparametric")
    stop("Cannot compute likelihood for a non-parametric fit.")
  
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
  
  if(object$model_type == "nonparametric")
    stop("Cannot compute likelihood for a non-parametric fit.")
  
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
