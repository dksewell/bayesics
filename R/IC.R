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
#' within \eqn{\pm}\code{mc_error}.  For WAIC, this is based on 
#' extrapolating the standard error from the preliminary 
#' posterior samples and may be inaccurate (at least 2000 samples 
#' will be used in final calculation).
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
    phi = as.vector(unlist(theta_draws[,p + 1]))
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
      phi = as.vector(unlist(theta_draws[,p + 1]))
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
  
  # Create helper functions
  ## Compute WAIC from matrix of individual log likelihood draws
  waic_from_matrix = function(llik_sub){
    
    lik_sub <- exp(llik_sub)
    
    lppd_i =
      apply(lik_sub, 1, function(x) log(mean(x)))
    
    # WAIC penalty
    p_waic_i <- apply(llik_sub, 1, var)
    
    # WAIC
    -2.0 * sum(lppd_i - p_waic_i)
  }
  
  ## Estimate MC variance using batches
  waic_mc_variance <- function(llik_i_draws,
                               alpha = 0.2,
                               R = 200){
    
    N <- nrow(llik_i_draws)
    S <- ncol(llik_i_draws)
    
    m <- floor(alpha * S)
    K <- floor(S / m)
    
    var_estimates <- numeric(R)
    
    for(r in 1:R){
      # random permutation of draws
      perm <- sample(S)
      
      waic_batches <- numeric(K)
      
      for(k in 1:K){
        idx <- perm[((k - 1) * m + 1):(k * m)]
        llik_sub <- llik_i_draws[, idx, drop = FALSE]
        waic_batches[k] <- waic_from_matrix(llik_sub)
      }
      
      var_estimates[r] <- var(waic_batches)
    }
    
    list(
      m = m,
      K = K,
      R = R,
      mean_var = mean(var_estimates),
      sd_var = sd(var_estimates),
      all_var_estimates = var_estimates
    )
  }
  
  
  # Get posterior samples 
  ## Get preliminary draws
  
  ### Get draws
  theta_draws = 
    get_posterior_draws(object,
                        n_draws = 2000)
  ### Compute draws of linear predictor
  Xbeta_draws = 
    tcrossprod(X,as.matrix(theta_draws[,1:p]))
  ### Compute draws of phi
  if(ncol(theta_draws) > p){
    phi = as.vector(unlist(theta_draws[,p + 1]))
    if(object$family$family == "negbinom")
      phi = exp(phi)
  }else{
    phi = rep(1.0,2000)
  }
  ### Compute individual likelihood values
  llik_i_draws = 
    future.apply::future_sapply(1:nrow(theta_draws),
                                function(i){
                                  log_lik_function(y,
                                                   trials *
                                                     object$family$linkinv(eta = 
                                                                             drop(Xbeta_draws[,i]) + 
                                                                             os),
                                                   phi[i])
                                })
  
  
  ### Estimate the variance of smaller batches
  preliminary_results = 
    waic_mc_variance(llik_i_draws)
  
  ### Pick the 90th percentile of variance estimates from all splits
  conservative_est_of_mc_var = 
    quantile(preliminary_results$all_var_estimates,0.9)
  
  ### Estimate the number of samples needed for accurate WAIC estimation
  #     Note: This is unjustified extrapolation.  
  n_prelim_draws = ncol(llik_i_draws)
  n_draws = 
    (0.2 * n_prelim_draws) *
    conservative_est_of_mc_var *
    qnorm(0.5 * (1.0 - 0.99))^2 /
    mc_error^2
  n_draws = round(n_draws)
  
  
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
      phi = as.vector(unlist(theta_draws[,p + 1]))
      if(object$family$family == "negbinom")
        phi = exp(phi)
    }else{
      phi = rep(1.0,ncol(Xbeta_draws))
    }
    ### Compute individual likelihood values
    llik_i_draws = 
      cbind(
        llik_i_draws,
        future.apply::future_sapply((n_prelim_draws + 1):n_draws,
                                    function(i){
                                      log_lik_function(y,
                                                       trials *
                                                         object$family$linkinv(eta = 
                                                                                 drop(Xbeta_draws[,i]) + 
                                                                                 os),
                                                       phi[i])
                                    })
      )
  }
  
  # Compute WAIC
  waic_from_matrix(llik_i_draws)
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

