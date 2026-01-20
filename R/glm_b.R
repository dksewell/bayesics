#' Bayesian Generalized Linear Models
#' 
#' glm_b is used to fit linear models.  It can be used to carry out 
#' regression for gaussian, binomial, and poisson data.  Note that if 
#' the family is gaussian, this is just a wrapper for \code{lm_b}.
#' 
#' 
#' @param formula A formula specifying the model.
#' @param data A data frame in which the variables specified in the formula 
#' will be found. If missing, the variables are searched for in the standard way. 
#' However, it is strongly recommended that you use this argument so that other 
#' generics for bayesics objects work correctly.
#' @param family A description of the error distribution and link function 
#' to be used in the model. See \code{?}\link[stats]{glm} for more information. 
#' Currently implemented families are \code{binomial()}, \code{poisson()}, 
#' \code{negbinom()}, and \code{gaussian()} (this last acts as a wrapper for 
#' \code{lm_b}. If missing \code{family}, \code{glm_b} will try to infer 
#' the data type; negative binomial will be used for count data.
#' @param trials Either character naming the variable in \code{data} that 
#' corresponds to the number of trials in the binomial observations, or else 
#' an integer vector giving the number of trials for each observation.
#' @param prior character.  One of "zellner", "normal", or "improper", giving 
#' the type of prior used on the regression coefficients.
#' @param zellner_g numeric.  Positive number giving the value of "g" in Zellner's
#' g prior.  Ignored unless prior = "zellner". Default is the number of observations.
#' @param prior_beta_mean numeric vector of same length as regression coefficients 
#' (denoted p). Unless otherwise specified, automatically set to rep(0,p).  Ignored 
#' unless prior = "normal".
#' @param prior_beta_precision pxp matrix giving a priori precision matrix to be 
#' scaled by the residual precision.  Ignored 
#' unless prior = "normal".
#' @param prior_phi_mean For negative binomial distributed outcomes, an 
#' exponential distribution is used for the prior of the dispersion parameter 
#' \eqn{phi}, parameterized such that \eqn{\text{Var}(y) = \mu + \frac{\mu^2}{\phi}}, 
#' so that the prior on \eqn{\phi} is \eqn{\lambda e^{-\lambda \phi}}, where 
#' \eqn{\lambda} equals \eqn{1/}\code{prior_phi_mean}.
#' @param ROPE vector of positive values giving ROPE boundaries for each regression 
#' coefficient.  Optionally, you can not include a ROPE boundary for the intercept. 
#' If missing, defaults go to those suggested by Kruchke (2018).
#' @param CI_level numeric. Credible interval level.
#' @param algorithm Either "VB" (default) for fixed-form variational Bayes, 
#' "IS" for importance sampling, or "LSA" for large sample approximation.
#' @param vb_maximum_iterations if \code{algorithm = "VB"}, the number of 
#' iterations used in the fixed-form VB algorithm.
#' @param proposal_df degrees of freedom used in the multivariate t proposal 
#' distribution if \code{algorithm = "IS"}.
#' @param seed integer.  Always set your seed!!!  Not used for 
#'  \code{algorithm = LSA}. 
#' @param mc_error If importance sampling is used, the number of posterior 
#' draws will ensure that with 99% probability the bounds of the credible 
#' intervals will be within \eqn{\pm} \code{mc_error}.
#' @param save_memory logical.  If TRUE, a more memory efficient approach 
#' will be taken at the expense of computataional time (for important 
#' sampling only.  But if memory is an issue, it's probably because you have a 
#' large sample size, in which case the normal approximation sans IS should 
#' probably work.)
#' 
#' 
#' @returns glm_b() returns an object of class "glm_b", which behaves as a list with 
#' the following elements:
#' \itemize{
#'  \item \code{summary} - tibble giving results for regression coefficients
#'  \item \code{posterior_draws}
#'  \item \code{ROPE}
#'  \item \code{hyperparameters} - list giving the user input or default hyperparameters used
#'  \item \code{fitted} - posterior mean of the individuals' means
#'  \item \code{residuals} - posterior mean of the residuals
#'  \item If \code{algorithm = "IS"}, the following:
#'  \itemize{
#'    \item \code{proposal_draws} - draws from 
#'  the importance sampling proposal distribution (i.e., multivariate 
#'  t centered at the posterior mode with precision equal to the 
#'  negative hessian, and degrees of freedom set to the user input \code{proposal_df}.
#'    \item \code{importance_sampling_weights} - importance sampling 
#'  weights that match the rows of the returned \code{proposal_draws}.
#'    \item \code{effective_sample_size} 
#'    \item \code{mc_error} 
#'  }
#'  \item other inputs into \code{glm_b}
#' }
#' 
#' \strong{Importance sampling:}
#' 
#' \code{glm_b} will, unless \code{use_importance_sampling = FALSE}, perform importance sampling. 
#' The proposal will use a multivariate t distribution, centered at the 
#' posterior mode, with the negative hessian as its precision matrix.  Do NOT
#' treat the proposal_draws as posterior draws.
#' 
#' \strong{Priors:}
#' 
#' If the prior is set to be either "zellner" or "normal", a normal distribution 
#' will be used as the prior of \eqn{\beta}, specifically 
#' \deqn{\beta \sim N(\mu, V)}
#' where \eqn{\mu} is the prior_beta_mean and V is the prior_beta_precision (not covariance) matrix. 
#' \itemize{
#'  \item{\code{zellner}: \code{glm_b} sets \eqn{\mu=0} and \eqn{V = \frac{1}{g} X^{\top} X}.}
#'  \item{\code{normal}: If missing \code{prior_beta_mean}, \code{glm_b} sets \eqn{\mu=0}, 
#'  and if missing \code{prior_beta_precision} V will be a diagonal matrix.  The first element, 
#'  corresponding to the intercept, will be \eqn{(2.5\times \max{\tilde{s}_y,1})^{-2}}, where 
#'  \eqn{\tilde{s}_y} is max of 1 and the standard deviation of \eqn{y}.  Remaining diagonal elements 
#'  will equal \eqn{(2.5 s_y/s_x)^{-2}}, where \eqn{s_x} is the standard deviations 
#'  of the covariates.  This equates to being 95% certain a priori that a change in 
#'  x by one standard deviation (\eqn{s_x}) would not lead to a change in the linear predictor of 
#'  more than 5 standard deviations (\eqn{5s_y}).  This imposes weak regularization that adapts to the scale 
#'  of the data elements.}
#' }
#' 
#' \strong{ROPE:}
#' 
#' If missing, the ROPE bounds will be given under the principle of "half of a 
#' small effect size."  
#' \itemize{
#'  \item{Gaussian.  Using Cohen's D of 0.2 as a small effect size, the ROPE is 
#'  built under the principle that moving the full range of X (i.e., \eqn{\pm 2 s_x}) 
#'  will not move the mean of y by more than the overall mean of \eqn{y} 
#'  minus \eqn{0.1s_y} to the overall mean of \eqn{y} plus \eqn{0.1s_y}.  
#'  The result is a ROPE equal to \eqn{|\beta_j| < 0.05s_y/s_j}.  If the covariate is 
#'  binary, then this is simply \eqn{|\beta_j| < 0.2s_y}.}
#'  \item{Poisson. FDA guidance suggests a small effect is a rate ratio less 
#'  than 1.25.  We use half this effect: 1.125, and consider ROPE to indicate 
#'  that a moving the full range of X (\eqn{\pm 2s_x} will not change the rate 
#'  ratio by more than this amount.  Thus the ROPE for the regression 
#'  coefficient equals \eqn{|\beta| < \frac{\log(1.125)}{4s_x}}. For binary 
#'  covariates, this is simply \eqn{|\beta| < \log(1.125)}.}
#' }
#' 
#' 
#' @examples
#' \donttest{
#' # Generate some negative-binomial data
#' set.seed(2025)
#' N = 500
#' test_data =
#'   data.frame(x1 = rnorm(N),
#'              x2 = rnorm(N),
#'              x3 = letters[1:5],
#'              time = rexp(N))
#' test_data$outcome =
#'   rnbinom(N,
#'           mu = exp(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e"))) * test_data$time,
#'           size = 0.7)
#' 
#' # Fit using variational Bayes (default)
#' fit_vb1 <-
#'   glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
#'         data = test_data,
#'         family = negbinom(),
#'         seed = 2025)
#' fit_vb1
#' summary(fit_vb1,
#'         CI_level = 0.9)
#' plot(fit_vb1)
#' coef(fit_vb1)
#' credint(fit_vb1,
#'         CI_level = 0.99)
#' bayes_factors(fit_vb1,
#'               by = "v")
#' preds =
#'   predict(fit_vb1)
#' 
#' # Try different priors
#' fit_vb2 <-
#'   glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
#'         data = test_data,
#'         family = negbinom(),
#'         seed = 2025,
#'         prior = "normal")
#' fit_vb2
#' fit_vb3 <-
#'   glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
#'         data = test_data,
#'         family = negbinom(),
#'         seed = 2025,
#'         prior = "improper")
#' fit_vb3
#' 
#' # Use Importance sampling instead of VB
#' fit_is <-
#'   glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
#'         data = test_data,
#'         family = negbinom(),
#'         algorithm = "IS",
#'         seed = 2025)
#' summary(fit_is)
#' 
#' # Use large sample approximation instead of VB
#' fit_lsa <-
#'   glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
#'         data = test_data,
#'         family = negbinom(),
#'         algorithm = "LSA",
#'         seed = 2025)
#' summary(fit_lsa)
#' }
#' 
#' @references 
#' 
#' Kruschke JK. Rejecting or Accepting Parameter Values in Bayesian Estimation. Advances in Methods and Practices in Psychological Science. 2018;1(2):270-280. doi:10.1177/2515245918771304
#' 
#' Tim Salimans. David A. Knowles. "Fixed-Form Variational Posterior Approximation through Stochastic Linear Regression." Bayesian Anal. 8 (4) 837 - 882, December 2013. https://doi.org/10.1214/13-BA858
#' @export


glm_b = function(formula,
                 data,
                 family,
                 trials,
                 prior = c("zellner","normal","improper")[1],
                 zellner_g,
                 prior_beta_mean,
                 prior_beta_precision,
                 prior_phi_mean = 1.0,
                 ROPE,
                 CI_level = 0.95,
                 vb_maximum_iterations = 1000,
                 algorithm = "VB",
                 proposal_df = 5,
                 seed = 1,
                 mc_error = 0.01,
                 save_memory = FALSE){
  
  set.seed(seed)
  
  # Extract
  if(missing(data)){
    mframe = 
      model.frame(formula)
  }else{
    mframe = 
      model.frame(formula,data)
  }
  y = model.response(mframe)
  if(is.factor(y)){
    y = as.integer(y)
    if(length(unique(y)) == 2) y = y - 1
  }
  
  X = model.matrix(formula,mframe)
  os = model.offset(mframe)
  N = nrow(X)
  p = ncol(X)
  if(ncol(X) > 1)
    s_j = apply(X[,-1,drop = FALSE],2,sd)
  s_y = sd(y)
  alpha = 1 - CI_level
  if(is.null(os)) os = numeric(N)
  
  
  # Get correct family
  ## Check data type for missing family
  if(missing(family)){
    if(isTRUE(all.equal(sort(unique(y)),
                        0:1))){
      family = binomial()
      warning("Family was not supplied.  Using binomial().")
    }else{
      if(isTRUE(all.equal(y,round(y))) & !any(y < 0)){
        family = negbinom()
        warning("Family was not supplied.  Using negbinom().")
      }else{
        family = gaussian()
        warning("Family was not supplied.  Using gaussian().")
      }
    }
  }
  
  ## If family is supplied, finesse it.
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  # Get prior rate for phi
  if(family$family == "negbinom")
    prior_phi_rate = 1.0 / prior_phi_mean
  
  # Get prior on $\beta$
  prior =
    c("zellner","normal","improper")[pmatch(tolower(prior),
                                            c("zellner","normal","improper"),duplicates.ok = FALSE)]
  
  # Get algorithm
  algorithm =
    c("VB","IS","LSA")[pmatch(toupper(algorithm),
                              c("VB","IS","LSA"),
                              duplicates.ok = FALSE)]
  
  if( ((family$family == "poisson") & (family$link != "log")) | 
      ((family$family == "binomial") & (family$link != "logit")) |
      !(family$family %in% c("poisson","binomial","negbinom")) ){
    algorithm = "IS"
  }
  
  # Send to lm_b if gaussian.  Else proceed.
  if(family$family == "gaussian"){
    # if(prior == "zellner"){
    #   if(missing(zellner_g)){
    #     message("\nThe g hyperparameter in Zellner's g prior is not specified.  It will be set automatically to n.\n")
    #     zellner_g = N
    #   }
    #   return(
    #     lm_b(formula = formula,
    #          data = data,
    #          prior = "zellner",
    #          zellner_g = zellner_g,
    #          )
    #   )
    # }
    # 
    prior = ifelse(prior == "normal","conjugate",prior)
    # return(
    #   lm_b(formula = formula,
    #        data = data,
    #        prior = prior,)
    # )
    lm_b_args = 
      list(formula = formula,
           prior = prior,
           CI_level = CI_level)
    if(!missing(data)) lm_b_args$data = data
    if(!missing(zellner_g)) lm_b_args$zellner_g = zellner_g
    if(!missing(prior_beta_mean)) lm_b_args$prior_beta_mean = prior_beta_mean
    if(!missing(prior_beta_precision)) lm_b_args$prior_beta_precision = prior_beta_precision
    if(!missing(ROPE)) lm_b_args$ROPE = ROPE
    return(
      do.call(lm_b,
              lm_b_args)
    )
    
  }else{
    
    # Get ROPE 
    if(missing(ROPE)){
      
      if( ((family$family == "poisson") & (family$link == "log")) | 
          ((family$family == "binomial") & (family$link == "logit")) | 
          ((family$family == "negbinom") & (family$link == "log")) ){
        ROPE = NA
        if(ncol(X) > 1){
          ROPE = 
            c(ROPE,
              log(1.0 + 0.25/2) / ifelse(apply(X[,-1,drop = FALSE],2,
                                               function(z) isTRUE(all.equal(0:1,
                                                                            sort(unique(z))))),
                                         1.0,
                                         4 * s_j))
        } 
        # From Kruchke (2018) on rate ratios from FDA <1.25. (Use half of small effect size for ROPE, hence 0.25/2) 
        # Use the same thing for odds ratios.
        # So this is the change due to moving through the range of x (\pm 2s_X).
        # For binary (or one-hot) use 1.
      }else{
        ROPE = NA
      }
      
    }else{
      if( !(length(ROPE) %in% (ncol(X) - 1:0))) stop("Length of ROPE, if supplied, must match the number of regression coefficients (or one less, if no ROPE is for the intercept).")
      if( any(ROPE) <= 0 ) stop("User supplied ROPE values must be positive.  The ROPE is assumed to be +/- the user supplied values.")
      if(length(ROPE) == ncol(X) - 1) ROPE = c(NA,ROPE)
    }
    if(family$family == "negbinom"){
      ROPE = c(ROPE, NA)
    }
    ROPE_bounds = 
      paste("(",-round(ROPE,3),",",round(ROPE,3),")",sep="")
    
    
    if(prior == "zellner"){
      if(missing(zellner_g)){
        message("\nThe g hyperparameter in Zellner's g prior is not specified.  It will be set automatically to n.\n")
        zellner_g = N
      }
      
      prior_beta_mean = rep(0,p)
      prior_beta_precision = 1/zellner_g * crossprod(X)
      
    }
    
    if(prior == "normal"){
      
      if(missing(prior_beta_mean)){
        message("\nThe mu hyperparameter in the normal prior is not specified.  It will be set automatically to 0.\n")
        prior_beta_mean = rep(0,p)
      }
      
      if(missing(prior_beta_precision)){
        message("\nThe V hyperparameter in the normal prior is not specified.  It will be set automatically to 4/25Diag(s^2_{X_j})")
        
        if(ncol(X) > 1){
          prior_beta_precision =
            diag(1.0 / c(2.5 * max(s_y,1),2.5 * s_y / s_j)^2)
        }else{
          prior_beta_precision = 
            matrix(1.0 / (2.5 * s_y)^2,1,1)
        }
        
      }else{
        if(ncol(X) == 1){
          prior_beta_precision = matrix(prior_beta_precision,1,1)
        }
      }
      
    }
    
    if(prior == "improper"){
      prior_beta_mean = prior_beta_precision = NA
    }
    
    
    # Check for errors in family type and outcome
    ## Binomial
    if(family$family == "binomial"){
      if(length(unique(y)) == 2){
        if( (class(y) %in% c("numeric","integer") ) &
            ( !isTRUE(all.equal(0:1,sort(unique(y)))) ) ){
          paste0("Treating ",
                 min(y),
                 " as '0' and ",
                 max(y),
                 " as '1'") |>
            message()
          y = ifelse(y == min(y),0,1)
        }
        if(is.character(y)) y = factor(y)
        if(is.factor(y)){
          paste0("Treating ",
                 levels(y)[2],
                 " as '1' and ",
                 levels(y)[1],
                 " as '0'") |>
            message()
          y = ifelse(y == levels(y)[1],0,1)
        }
      }
      if(missing(trials)){
        message("Assuming all observations correspond to Bernoulli, i.e., Binomial with one trial.")
        trials = rep(1.0,N)
      }else{
        if(is.character(trials)) trials = data[[trials]]
        trials = as.numeric(trials)
      }
    }
    
    ## Poisson/NB
    if(family$family %in% c("poisson","negbinom")){
      if( !(class(y) %in% c("numeric","integer")) ) stop("Outcome must be numeric or integer.")
      if(min(y) < 0) stop("Minumum of outcome must not be negative")
      trials = rep(1.0,N)
    }
    
    
    
    # Get initial start from glm
    ## Binomial
    if(family$family == "binomial"){
      if(ncol(X) > 1){
        init = glm(cbind(y,trials) ~ X[,-1] + offset(os),
                   family)
      }else{
        init = glm(cbind(y,trials) ~ 1 + offset(os),
                   family)
      }
    }
    ## Poisson
    if(family$family == "poisson"){
      if(ncol(X) > 1){
        init = glm(y ~ X[,-1] + offset(os),
                   family)
      }else{
        init = glm(y ~ 1 + offset(os),
                   family)
      }
    }
    ## Negative Binomial
    if(family$family == "negbinom"){
      if(ncol(X) > 1){
        init = glm(y ~ X[,-1] + offset(os),
                   family = poisson())
      }else{
        init = glm(y ~ 1 + offset(os),
                   poisson())
      }
    }
    
    # Get actual posterior mode from this
    ## Create posterior function
    log_posterior = function(x){
      beta_coefs = x[1:c(length(x) - 
                           (family$family == "negbinom"))]
      eta = drop(X %*% beta_coefs) + os
      mu = family$linkinv(eta)
      
      if(family$family == "binomial"){
        lpost =
          sum(
            dbinom(y,trials,mu,log = TRUE)
          )
      }
      if(family$family == "poisson"){
        lpost =
          sum(
            dpois(y,mu,log = TRUE)
          )
      }
      if(family$family == "negbinom"){
        phi = exp(x[length(x)]) # Use log transformation of dispersion for optim()
        lpost =
          sum(
            dnbinom(y,mu = mu, size = phi,log = TRUE)
          )
      }
      
      if(prior != "improper"){
        lpost =
          lpost -
          0.5 *
          drop(crossprod(
            beta_coefs - prior_beta_mean,
            prior_beta_precision %*% (beta_coefs - prior_beta_mean)
          ))
      }
      if(family$family == "negbinom"){
        lpost =
          lpost -
          prior_phi_rate * phi
      }
      
      return(lpost)
    }
    
    ## Create gradient
    nabla_log_posterior = function(x){
      beta_coefs = x[1:c(length(x) - 
                           (family$family == "negbinom"))]
      eta = drop(X %*% beta_coefs) + os
      mu = family$linkinv(eta)
      
      if(family$family == "binomial"){
        grad_lpost =
          drop( (y - trials * mu) %*% X )
      }
      if(family$family == "poisson"){
        grad_lpost =
          drop( (y - mu) %*% X )
      }
      if(family$family == "negbinom"){
        phi = exp(x[length(x)]) # Use log transformation of dispersion for optim()
        grad_lpost =
          c(
            drop(( phi * (y - mu) / (mu + phi) ) %*% X),
            sum(
              digamma(y + phi) - 
                (y + phi) / (mu + phi)  - 
                log(mu + phi)
            ) + 
              N * (1.0 - 
                     digamma(phi) + 
                     log(phi)
              ) -
              prior_phi_rate
          )
      }
      
      if(prior != "improper"){
        grad_lpost[1:c(length(x) - 
                         (family$family == "negbinom"))] =
          grad_lpost[1:c(length(x) - 
                           (family$family == "negbinom"))] -
          drop( prior_beta_precision %*% (beta_coefs - prior_beta_mean) )
      }
      
      return(drop(grad_lpost))
    }
    
    # Get posterior mode and hessian
    ## Run optim
    if( ( (family$family == "binomial") & (family$link == "logit") ) | 
        ( (family$family == "poisson") & (family$link == "log") ) | 
        ( (family$family == "negbinom") & (family$link == "log") ) ){
      opt = 
        optim(c(coef(init),
                log(prior_phi_mean))[1:c(ncol(X) + 
                                           (family$family == "negbinom"))],
              log_posterior,
              nabla_log_posterior,
              method = "BFGS",
              hessian = TRUE,
              control = list(fnscale = -1))
    }else{
      opt = 
        optim(c(coef(init),
                log(prior_phi_mean))[1:c(ncol(X) + 
                                           (family$family == "negbinom"))],
              log_posterior,
              method = "Nelder-Mead",
              hessian = TRUE,
              control = list(fnscale = -1))
    }
    ## Compute covariance matrix
    covmat = NULL
    try({
      covmat = 
        chol2inv(chol(-opt$hessian))
    }, silent=TRUE)
    if(is.null(covmat)){
      try({
        covmat = 
          qr.solve(-opt$hessian)
      }, silent=TRUE)
    }
    if(is.null(covmat)){
      try({
        covmat = 
          solve(-opt$hessian)
      }, silent=TRUE)
    }
    if(is.null(covmat)) stop("Hessian is not invertible.")
    
    
    if(algorithm == "LSA"){
      # Collate results from large sample approx to return
      return_object = list()
      ## Summary
      return_object$summary = 
        tibble::tibble(Variable = 
                         c(colnames(X),
                           "log(phi)")[1:c(ncol(X) + 
                                             (family$family == "negbinom"))],
                       `Post Mean` = opt$par,
                       Lower = 
                         qnorm(alpha / 2,
                               opt$par,
                               sd = sqrt(diag(covmat))),
                       Upper = 
                         qnorm(1 - alpha / 2,
                               opt$par,
                               sd = sqrt(diag(covmat))),
                       `Prob Dir` = 
                         pnorm(0,
                               opt$par,
                               sd = sqrt(diag(covmat))),
                       ROPE = 
                         pnorm(ROPE,
                               opt$par,
                               sqrt(diag(covmat))) - 
                         pnorm(-ROPE,
                               opt$par,
                               sqrt(diag(covmat))),
                       `ROPE bounds` = ROPE_bounds)
      
      return_object$summary$`Prob Dir` = 
        ifelse(return_object$summary$`Prob Dir` > 0.5,
               return_object$summary$`Prob Dir`,
               1.0 - return_object$summary$`Prob Dir`)
      
      ## Posterior covariance matrix
      return_object$posterior_covariance = covmat
      
      
    }#End: large sample summary
    if(algorithm == "IS"){
      # Perform importance sampling
      
      ## Get faster posterior function
      log_posterior_multiple_samples = function(draws){
        eta = tcrossprod(X, draws[,1:ncol(X)]) + os
        mu = family$linkinv(eta)
        
        if(family$family == "binomial"){
          lpost = 
            drop(
              y %*% log(mu) + 
                (trials - y) %*% log(1.0 - mu)
            )
          
        }
        if(family$family == "poisson"){
          lpost =
            drop(y %*% log(mu)) - colSums(mu)
        }
        
        if(prior != "improper"){
          if(ncol(X) > 1){
            lpost =
              lpost + 
              mvtnorm::dmvnorm(draws[,1:ncol(X)],
                               prior_beta_mean,
                               qr.solve(prior_beta_precision),
                               log = TRUE)
          }else{
            lpost =
              lpost + 
              dnorm(draws[,1:ncol(X)],
                    prior_beta_mean,
                    sqrt(1.0/drop(prior_beta_precision)),
                    log = TRUE)
          }
        }
        
        return(lpost)
      }
      
      
      ## Perform preliminary draws
      ### Get draws from the proposal
      proposal_draws = 
        mvtnorm::rmvt(500,
                      delta = opt$par,
                      sigma = covmat,
                      df = proposal_df)
      ### Get IS weights
      if( save_memory | 
          (family$family == "negbinom") ){
        is_weights = 
          sapply(1:500,
                 function(i){
                   log_posterior(proposal_draws[i,])
                 }
          )
      }else{
        is_weights = 
          log_posterior_multiple_samples(proposal_draws)
      }
      is_weights = 
        is_weights -
        sum(
          mvtnorm::dmvt(proposal_draws,
                        delta = opt$par,
                        sigma = covmat,
                        df = proposal_df,
                        log = TRUE)
        )
      
      ### Stabilize and Normalize
      is_weights = exp(is_weights - max(is_weights))
      is_weights = is_weights / sum(is_weights)
      ### Perform SIR (providing upper bound on number of post draws needed)
      new_draws = 
        proposal_draws[sample(500,500,T,is_weights),]
      if(NCOL(new_draws) == 1) new_draws = matrix(new_draws,ncol=1)
      fhats = 
        future.apply::future_lapply(1:NCOL(new_draws),
                                    function(i){
                                      density(new_draws[,i],adjust = 2)
                                    })
      
      n_draws = 
        future.apply::future_sapply(1:NCOL(new_draws),
                                    function(i){
                                      0.5 * alpha * (1.0 - 0.5 * alpha) *
                                        (
                                          qnorm(0.5 * (1.0 - 0.99)) / 
                                            mc_error /
                                            fhats[[i]]$y[which.min(abs(fhats[[i]]$x - 
                                                                         quantile(new_draws[,i], 0.5 * alpha)))]
                                        )^2
                                    }) |> 
        max() |> 
        round()
      
      
      ## Perform final draws
      ## Get draws from the proposal
      proposal_draws = 
        mvtnorm::rmvt(n_draws,
                      delta = opt$par,
                      sigma = covmat,
                      df = proposal_df)
      # Get IS weights
      if( save_memory | 
          (family$family == "negbinom") ){
        is_weights = 
          future.apply::future_sapply(1:n_draws,
                                      function(i){
                                        log_posterior(proposal_draws[i,])
                                      }
          )
      }else{
        is_weights = 
          log_posterior_multiple_samples(proposal_draws)
        gc()
      }
      is_weights = 
        is_weights -
        sum(
          mvtnorm::dmvt(proposal_draws,
                        delta = opt$par,
                        sigma = covmat,
                        df = proposal_df,
                        log = TRUE)
        )
      
      ## Stabilize and Normalize
      is_weights = exp(is_weights - max(is_weights))
      is_weights = is_weights / sum(is_weights)
      
      ## Get effective sample size
      ESS = 1.0 / sum( is_weights^2 )
      
      # Create helper function to get CI bounds
      CI_from_weighted_sample = function(x,w){
        w = cumsum(w[order(x)])
        x = x[order(x)]
        LB = max(which(w <= 0.5 * alpha))
        UB = min(which(w >= 1.0 - 0.5 * alpha))
        return(c(lower = x[LB],
                 upper = x[UB]))
      }
      ## Go ahead and get bounds
      CI_bounds = 
        apply(proposal_draws,2,
              CI_from_weighted_sample,
              w = is_weights)
      
      # Create helper function to get ROPE
      ROPE_from_weighted_sample = function(x,w,R){
        sum(w[(x > -R) & (x < R)])
      }
      
      # Create helper function to get PDir
      PDir_from_weighted_sample = function(x,w){
        max(sum(w[x > 0]),
            sum(w[x < 0]))
      }
      
      results =
        tibble::tibble(Variable = 
                         c(colnames(X),
                           "log(phi)")[1:c(ncol(X) + 
                                             (family$family == "negbinom"))], 
                       `Post Mean` = apply(proposal_draws,2,
                                           weighted.mean,
                                           w = is_weights),
                       Lower = CI_bounds["lower",],
                       Upper = CI_bounds["upper",],
                       `Prob Dir` = 
                         apply(proposal_draws,2,
                               PDir_from_weighted_sample,
                               w = is_weights),
                       ROPE =
                         sapply(1:ncol(proposal_draws),
                                function(j){
                                  ROPE_from_weighted_sample(proposal_draws[,j],
                                                            is_weights,
                                                            ROPE[j])
                                }),
                       `ROPE bounds` = ROPE_bounds
        )
      
      return_object =
        list(summary = results,
             proposal_draws = proposal_draws,
             importance_sampling_weights = is_weights,
             effective_sample_size = ESS,
             mc_error = mc_error)
      
      
    }#End: importance sampling
    if(algorithm == "VB"){
      
      ## Get Hessian of log posterior
      hessian_log_posterior = function(x){ 
        beta_coefs = x[1:c(length(x) - 
                             (family$family == "negbinom"))]
        eta = drop(X %*% beta_coefs) + os
        mu = family$linkinv(eta)
        
        if(family$family == "binomial"){
          hessian_lpost =
            crossprod(X,
                      Diagonal(x = -trials * mu * (1.0 - mu)) %*% X)
        }
        if(family$family == "poisson"){
          hessian_lpost =
            crossprod(X,
                      Diagonal(x = -mu) %*% X)
        }
        if(family$family == "negbinom"){
          phi = exp(x[length(x)])
          d2_bb = 
            crossprod(X,
                      Diagonal(x = -phi * mu *(phi + y) / (phi + mu)^2) %*% X)
          d2_pp = 
            sum(trigamma(y + phi)) -
            N * trigamma(phi) + 
            N / phi -
            sum( 2.0 / (mu + phi) ) + 
            sum( (y + phi) / (mu + phi)^2 )
          d2_bp = 
            ( mu * (y - mu) / (mu + phi)^2 ) %*% X
          hessian_lpost = matrix(0.0,ncol(X)+1,ncol(X)+1)
          hessian_lpost[1:ncol(X),1:ncol(X)] = as.matrix(d2_bb)
          hessian_lpost[ncol(X) + 1, ncol(X) + 1] = d2_pp
          hessian_lpost[1:ncol(X),ncol(X) + 1] = drop(d2_bp)
          hessian_lpost[ncol(X) + 1,1:ncol(X)] = drop(d2_bp)
        }
        
        if(prior != "improper"){
          hessian_lpost[1:c(length(x) - 
                              (family$family == "negbinom")),
                        1:c(length(x) - 
                              (family$family == "negbinom"))]  =
            hessian_lpost[1:c(length(x) - 
                                (family$family == "negbinom")),
                          1:c(length(x) - 
                                (family$family == "negbinom"))]  - 
            prior_beta_precision
        }
        
        return(hessian_lpost)
      }
      
      ## Initialize
      m = opt$par
      V = covmat
      z = m
      P = qr.solve(V)
      a = 0.0
      zbar = 0.0
      Pbar = matrix(0.0,nrow(P),ncol(P))
      abar = 0.0
      
      ## Run algo 2
      step_size = 1.0 / sqrt(vb_maximum_iterations)
      for(i in 1:vb_maximum_iterations){
        beta_draw = 
          mvtnorm::rmvnorm(1,
                           m,
                           V) |> 
          drop()
        
        g = nabla_log_posterior(beta_draw)
        H = hessian_log_posterior(beta_draw)
        
        a = (1.0 - step_size) * a + step_size * g
        P = (1.0 - step_size) * P - step_size * H
        z = (1.0 - step_size) * z + step_size * beta_draw
        V = qr.solve(P)
        m = V %*% a + z
        
        if(i > 0.5 * vb_maximum_iterations){
          abar = abar + 2.0 / vb_maximum_iterations * g
          Pbar = Pbar - 2.0 / vb_maximum_iterations * H
          zbar = zbar + 2.0 / vb_maximum_iterations * beta_draw
        }
      }
      V = qr.solve(Pbar)
      m = drop(V %*% abar + zbar)
      
      
      return_object = list()
      ## Summary
      return_object$summary = 
        tibble(Variable = 
                 c(colnames(X),
                   "log(phi)")[1:c(ncol(X) + 
                                     (family$family == "negbinom"))],
               `Post Mean` = m,
               Lower = 
                 qnorm(alpha / 2,
                       m,
                       sd = sqrt(diag(V))),
               Upper = 
                 qnorm(1 - alpha / 2,
                       m,
                       sd = sqrt(diag(V))),
               `Prob Dir` = 
                 pnorm(0,
                       m,
                       sd = sqrt(diag(V))),
               ROPE = 
                 pnorm(ROPE,
                       m,
                       sqrt(diag(V))) - 
                 pnorm(-ROPE,
                       m,
                       sqrt(diag(V))),
               `ROPE bounds` = ROPE_bounds)
      
      return_object$summary$`Prob Dir` = 
        ifelse(return_object$summary$`Prob Dir` > 0.5,
               return_object$summary$`Prob Dir`,
               1.0 - return_object$summary$`Prob Dir`)
      
      ## Posterior covariance matrix
      return_object$posterior_covariance = V
      
    }
    
    
    
    
    # Return hyperparameters and other user-inputs
    if(prior != "improper"){
      return_object$hyperparameters = 
        list(prior_beta_mean = prior_beta_mean,
             prior_beta_precision = prior_beta_precision)
    }else{
      return_object$hyperparameters = NA
    }
    return_object$formula = formula
    if(missing(data)){
      return_object$data = mframe
    }else{
      return_object$data = data
    }
    return_object$family = family
    return_object$algorithm = algorithm
    return_object$prior = prior
    return_object$ROPE = ROPE
    return_object$trials = trials
    return_object$CI_level = CI_level
    
    
    # Get fitted objects
    eta = drop(X %*% return_object$summary$`Post Mean`[1:ncol(X)]) + os
    mu = family$linkinv(eta)
    return_object$fitted =
      trials * mu
    
    # Get Pearson residuals
    if(family$family == 'negbinom'){
      SDs = 
        sqrt(
          trials * family$variance(mu,
                                   exp(return_object$summary$`Post Mean`[ncol(X)]))
        )
    }else{
      SDs = 
        sqrt(
          trials * family$variance(mu)
        )
    }
    return_object$residuals = 
      (y - trials * mu) / SDs
    
    rownames(return_object$summary) = NULL
    
    
    # attach helpers for generics
    return_object$terms = terms(mframe)
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
                     class = "glm_b"))
  }
}













