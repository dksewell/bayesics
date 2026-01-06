#' Bayesian Linear Models
#' 
#' lm_b is used to fit linear models.  It can be used to carry out 
#' regression, single stratum analysis of variance and analysis of 
#' covariance (although aov_b may provide a more convenient 
#' interface for ANOVA.)
#' 
#' @param formula A formula specifying the model.
#' @param data A data frame in which the variables specified in the formula 
#' will be found. If missing, the variables are searched for in the standard way.  
#' However, it is strongly recommended that you use this argument so that other 
#' generics for bayesics objects work correctly.
#' @param weights an optional vector of weights to be used in the fitting process. 
#' Should be NULL or a numeric vector. If non-NULL, it is assumed that the 
#' variance of \eqn{y_i} can be written as \eqn{Var(y_i) = \sigma^2/w_i}. While the 
#' estimands remain the same, the estimation is done by performing unweighted 
#' lm_b on \eqn{W^{\frac{1}{2}}y} and \eqn{W^{\frac{1}{2}}X}, where \eqn{W} is the 
#' diagonal matrix of weights.  Note that this then affects the zellner prior.
#' @param prior character.  One of "zellner", "conjugate", or "improper", giving 
#' the type of prior used on the regression coefficients.
#' @param zellner_g numeric.  Positive number giving the value of "g" in Zellner's
#' g prior.  Ignored unless prior = "zellner".
#' @param prior_beta_mean numeric vector of same length as regression coefficients 
#' (denoted p). Unless otherwise specified, automatically set to rep(0,p).  Ignored 
#' unless prior = "conjugate".
#' @param prior_beta_precision pxp matrix giving a priori precision matrix to be 
#' scaled by the residual precision.
#' @param prior_var_shape numeric. Twice the shape parameter for the inverse gamma prior on
#' the residual variance(s).  I.e., \eqn{\sigma^2\sim IG}(prior_var_shape/2,prior_var_rate/2).
#' @param prior_var_rate  numeric. Twice the rate parameter for the inverse gamma prior on
#' the residual variance(s).  I.e., \eqn{\sigma^2\sim IG}(prior_var_shape/2,prior_var_rate/2).
#' @param ROPE vector of positive values giving ROPE boundaries for each regression 
#' coefficient.  Optionally, you can not include a ROPE boundary for the intercept. 
#' If missing, defaults go to those suggested by Kruchke (2018).
#' @param CI_level numeric. Credible interval level.
#' 
#' @return lm_b() returns an object of class "lm_b", which behaves as a list with 
#' the following elements:
#' \itemize{
#'  \item summary - tibble giving results for regression coefficients
#'  \item posterior_parameters - list giving the posterior parameters
#'  \item hyperparameters - list giving the user input (or default) hyperparameters used
#'  \item fitted - posterior mean of the individuals' means
#'  \item residuals - posterior mean of the residuals
#'  \item formula, data - input by user
#' }
#' 
#' @details
#' 
#' \strong{MODEL:}
#' 
#' The likelihood is given by 
#' \deqn{
#'  y_i \overset{ind}{\sim} N(x_i'\beta,\sigma^2).
#' }
#' The prior is given by 
#' \deqn{
#'  \beta|\sigma^2 \sim N\left( \mu, \sigma^2 V^{-1} \right) \\
#'  \sigma^2 \sim \Gamma^{-1}(a/2,b/2).
#' }
#' \itemize{
#'  \item For Zellner's g prior, \eqn{\frac{1}{g}X'X}.
#'  \item The default for the conjugate prior is based on arguments from 
#'  standardized regression.  The default \eqn{V} is dicated by saying, "a priori, 
#'  we are 95% certain that a standard deviation increase in \eqn{X} will not lead 
#'  to more than a 5 standard deviation in the mean of \eqn{y}."  If we then 
#'  set the prior on the intercept to be flat, this leads to
#'  \deqn{
#'    V^{1/2} = diag(1,s_{x_1},\ldots,s_{x_p}) / 2.5,
#'  }
#'  where \eqn{s_y} is the standard deviation of \eqn{y}, and \eqn{s_{x_j}} is 
#'  the standard deviation of the \eqn{j^{th}} covariate.
#'  \item If \code{prior = "improper"}, then the prior is
#'  \deqn{
#'    \pi(\beta,\sigma^2) \propto \frac{1}{\sigma^2}.
#'  }
#' }
#' 
#' 
#' 
#' 
#' \strong{ROPE:}
#' 
#' If missing, the ROPE bounds will be given under the principle of "half of a 
#' small effect size."  Using Cohen's D of 0.2 as a small effect size, the ROPE is 
#'  built under the principle that moving the full range of X (i.e., \eqn{\pm 2 s_x}) 
#'  will not move the mean of y by more than the overall mean of \eqn{y} 
#'  minus \eqn{0.1s_y} to the overall mean of \eqn{y} plus \eqn{0.1s_y}.  
#'  The result is a ROPE equal to \eqn{|\beta_j| < 0.05s_y/s_j}.  If the covariate is 
#'  binary, then this is simply \eqn{|\beta_j| < 0.2s_y}.
#' 
#' @export


lm_b = function(formula,
                data,
                weights,
                prior = c("zellner","conjugate","improper")[1],
                zellner_g,
                prior_beta_mean,
                prior_beta_precision,
                prior_var_shape = 0.001,
                prior_var_rate = 0.001,
                ROPE,
                CI_level = 0.95){
  
  
  # Extract 
  if(missing(data)){
    m = 
      model.frame(formula)
  }else{
    m = 
      model.frame(formula,data)
  }
  if(missing(weights)) 
    weights = rep(1.0,nrow(m))
  if(is.character(weights))
    weights = data[[weights]]
  
  if(any(weights <=0 ))
    stop("weights must be strictly positive.")
  w_sqrt = sqrt(weights)
  y = 
    model.response(m)
  X = 
    model.matrix(formula,m)
  p = ncol(X)
  alpha = 1 - CI_level
  N = nrow(X)
  s_y = sd(y)
  if(ncol(X) > 1)
    s_j = apply(X[,-1,drop = FALSE],2,sd)
  
  # Get ROPE 
  if(missing(ROPE)){
    if(ncol(X) == 1){
      ROPE = NA
    }else{
      ROPE = 
        c(NA,
          0.2 * s_y / ifelse(apply(X[,-1,drop=FALSE],2,
                                   function(z) isTRUE(all.equal(0:1,
                                                                sort(unique(z))))),
                             1.0,
                             4.0 * s_j))
    }
    # From Kruchke (2018) on standardized regression
    # This considers a small change in y to be \pm 0.1s_y (half of Cohen's D small effect).
    # So this is the change due to moving through the range of x (\pm 2s_X).
    # For binary (or one-hot) use 1.
    
    
  }else{
    if( !(length(ROPE) %in% (ncol(X) - 1:0))) stop("Length of ROPE, if supplied, must match the number of regression coefficients (or one less, if no ROPE is for the intercept).")
    if( any(ROPE) <= 0 ) stop("User supplied ROPE values must be positive.  The ROPE is assumed to be +/- the user supplied values.")
    if(length(ROPE) == ncol(X) - 1) ROPE = c(NA,ROPE)
  }
  ROPE_bounds = 
    paste("(",-round(ROPE,3),",",round(ROPE,3),")",sep="")
  
  
  
  # Adjust for weights and precompute useful quantities
  y = w_sqrt * y
  X = diag(w_sqrt) %*% X
  XtX = crossprod(X)
  XtX_inv = qr.solve(XtX)
  s_y = sd(y)
  if(ncol(X) > 1)
    s_j = apply(X[,-1,drop = FALSE],2,sd)
  
  
  prior = 
    c("zellner","conjugate","improper")[pmatch(tolower(prior),
                                               c("zellner","conjugate","improper"),duplicates.ok = FALSE)]
  
  if(prior == "zellner"){
    if(missing(zellner_g)){
      message("\nThe g hyperparameter in Zellner's g prior is not specified.  It will be set automatically to n.\n")
      zellner_g = N
    }
    if(missing(prior_var_shape)){
      message("\nThe a hyperparameter in the inverse gamma prior is not specified.  It will be set automatically to 0.001.\n")
      a = 0.001
    }else{
      a = prior_var_shape
    }
    if(missing(prior_var_rate)){
      message("\nThe b hyperparameter in the inverse gamma prior is not specified.  It will be set automatically to 0.001.\n")
      b = 0.001
    }else{
      b = prior_var_rate
    }
    
    V = 1/zellner_g * XtX
    
    V_tilde = (zellner_g + 1.0) / zellner_g * XtX
    V_tilde_inv = qr.solve(V_tilde)
    mu_tilde = V_tilde_inv %*% crossprod(X,y)
    a_tilde = a + N
    b_tilde = b + crossprod(y) - crossprod(mu_tilde, V_tilde %*% mu_tilde)
    b_tilde = drop(b_tilde)
    
    results = 
      tibble::tibble(Variable = rownames(mu_tilde),
                     `Post Mean` = mu_tilde[,1],
                     Lower = extraDistr::qlst(alpha/2,
                                              a_tilde,
                                              mu_tilde,
                                              sqrt(b_tilde/a_tilde * diag(V_tilde_inv))),
                     Upper = extraDistr::qlst(1 - alpha/2,
                                              a_tilde,
                                              mu_tilde,
                                              sqrt(b_tilde/a_tilde * diag(V_tilde_inv))),
                     `Prob Dir` = 
                       sapply(extraDistr::plst(0,
                                               a_tilde,
                                               mu_tilde,
                                               sqrt(b_tilde/a_tilde * diag(V_tilde_inv))),
                              function(x) max(x,1-x)),
                     ROPE = 
                       extraDistr::plst(ROPE,
                                        a_tilde,
                                        mu_tilde,
                                        sqrt(b_tilde/a_tilde * diag(V_tilde_inv))) -
                       extraDistr::plst(-ROPE,
                                        a_tilde,
                                        mu_tilde,
                                        sqrt(b_tilde/a_tilde * diag(V_tilde_inv))),
                     `ROPE bounds` = ROPE_bounds
      )
    results$ROPE[1] = NA
    
    return_object = 
      list(summary = results,
           posterior_parameters = list(mu_tilde = mu_tilde,
                                       V_tilde = V_tilde,
                                       a_tilde = a_tilde,
                                       b_tilde = b_tilde),
           hyperparameters = list(mu = rep(0,p),
                                  V = V,
                                  a = a,
                                  b = b))
  }
  
  
  if(prior == "conjugate"){
    
    if(missing(prior_beta_mean)){
      message("\nThe mu hyperparameter in the normal prior is not specified.  It will be set automatically to 0.\n")
      mu = rep(0,p)
      if(all(dplyr::near(X[,1],1.0)))
        mu[1] = mean(y)
    }else{
      mu = prior_beta_mean
    }
    if(missing(prior_var_shape)){
      message("\nThe a hyperparameter in the inverse gamma prior is not specified.  It will be set automatically to 0.001.\n")
      a = 0.001
    }else{
      a = prior_var_shape
    }
    if(missing(prior_var_rate)){
      message("\nThe b hyperparameter in the inverse gamma prior is not specified.  It will be set automatically to 0.001.\n")
      b = 0.001
    }else{
      b = prior_var_rate
    }
    if(missing(prior_beta_precision)){
      message("\nThe V hyperparameter in the normal prior is not specified.  It will be set automatically to 4/25Diag(s^2_{X_j})")
      if(ncol(X) > 1){
        V = 
          diag(1.0 / c(2.5,2.5 / s_j)^2)
      }else{
        V = 
          matrix(1.0 / (2.5)^2,1,1)
      }
    }else{
      V = prior_beta_precision
      if(ncol(X) == 1){
        V = matrix(V,1,1)
      }
    }
    
    V_tilde = V + XtX
    V_tilde_inv = qr.solve(V_tilde)
    mu_tilde = V_tilde_inv %*% (V%*%mu + crossprod(X,y))
    a_tilde = a + N
    b_tilde = b + crossprod(y) - crossprod(mu_tilde, V_tilde %*% mu_tilde)
    b_tilde = drop(b_tilde)
    
    results = 
      tibble::tibble(Variable = rownames(mu_tilde),
                     `Post Mean` = mu_tilde[,1],
                     Lower = extraDistr::qlst(alpha/2,
                                              a_tilde,
                                              mu_tilde,
                                              sqrt(b_tilde/a_tilde * diag(V_tilde_inv))),
                     Upper = extraDistr::qlst(1 - alpha/2,
                                              a_tilde,
                                              mu_tilde,
                                              sqrt(b_tilde/a_tilde * diag(V_tilde_inv))),
                     `Prob Dir` = 
                       sapply(extraDistr::plst(0,
                                               a_tilde,
                                               mu_tilde,
                                               sqrt(b_tilde/a_tilde * diag(V_tilde_inv))),
                              function(x) max(x,1-x)),
                     ROPE = 
                       extraDistr::plst(ROPE,
                                        a_tilde,
                                        mu_tilde,
                                        sqrt(b_tilde/a_tilde * diag(V_tilde_inv))) -
                       extraDistr::plst(-ROPE,
                                        a_tilde,
                                        mu_tilde,
                                        sqrt(b_tilde/a_tilde * diag(V_tilde_inv))),
                     `ROPE bounds` = ROPE_bounds
      )
    results$ROPE[1] = NA
    
    return_object = 
      list(summary = results,
           posterior_parameters = list(mu_tilde = mu_tilde,
                                       V_tilde = V_tilde,
                                       a_tilde = a_tilde,
                                       b_tilde = b_tilde),
           hyperparameters = list(mu = mu,
                                  V = V,
                                  a = a,
                                  b = b))
  }
  
  
  if(prior == "improper"){
    
    mod = lm(formula,m)
    mu_tilde = coef(mod)
    Sigma = sigma(mod)^2 * XtX_inv
    results = 
      tibble::tibble(Variable = names(mu_tilde),
                     `Post Mean` = mu_tilde,
                     Lower = extraDistr::qlst(alpha/2,
                                              N - p,
                                              mu_tilde,
                                              sqrt(diag(Sigma))),
                     Upper = extraDistr::qlst(1 - alpha/2,
                                              N - p,
                                              mu_tilde,
                                              sqrt(diag(Sigma))),
                     `Prob Dir` = 
                       sapply(extraDistr::plst(0,
                                               N - p,
                                               mu_tilde,
                                               sqrt(diag(Sigma))),
                              function(x) max(x,1-x)),
                     ROPE = 
                       extraDistr::plst(ROPE,
                                        N - p,
                                        mu_tilde,
                                        sqrt(diag(Sigma))) -
                       extraDistr::plst(-ROPE,
                                        N - p,
                                        mu_tilde,
                                        sqrt(diag(Sigma))),
                     `ROPE bounds` = ROPE_bounds
      )
    results$ROPE[1] = NA
    
    return_object = 
      list(summary = results,
           posterior_parameters = list(mu_tilde = mu_tilde,
                                       V_tilde = XtX,
                                       Sigma = Sigma,
                                       a_tilde = N - p,
                                       b_tilde = sum(resid(mod)^2)),
           hyperparameters = NA)
  }
  
  
  return_object$fitted = 
    drop(X %*% return_object$summary$`Post Mean`)
  return_object$residuals = 
    drop(y - return_object$fitted)
  
  return_object$formula = formula
  if(missing(data)){
    return_object$data = m
  }else{
    return_object$data = data
  }
  return_object$prior = prior
  return_object$ROPE = ROPE
  return_object$CI_level = CI_level
  
  rownames(return_object$summary) = NULL
  
  # attach helpers for generics
  return_object$terms = terms(m)
  if(any(attr(return_object$terms,"dataClasses") %in% c("factor","character"))){
    return_object$xlevels = list()
    factor_vars = 
      names(attr(return_object$terms,"dataClasses"))[attr(return_object$terms,"dataClasses") %in% c("factor","character")]
    for(j in factor_vars){
      return_object$xlevels[[j]] = 
        unique(return_object$data[[j]])
    }
  }
  
  
  return(structure(return_object,
                   class = "lm_b"))
}