#' Bayesian model averaging 
#' 
#' Estimates and CIs from BMA
#' 
#' \code{bma_inference} leverages the \code{bms} function from its 
#' eponymous R package, and then uses \code{lm_b} to obtain inference 
#' on the regression coefficients for Bayesian model averaging.
#' 
#' @param formula A formula specifying the model.
#' @param data Data used in linear regression model
#' @param zellner_g numeric.  Positive number giving the value of "g" in Zellner's
#' g prior.  
#' @param CI_level Level for credible interval
#' @param ROPE vector of positive values giving ROPE boundaries for each regression 
#' coefficient.  Optionally, you can not include a ROPE boundary for the intercept. 
#' If missing, defaults go to those suggested by Kruchke (2018).
#' @param mcmc_draws Integer. Number of draws passed into \code{\link[BMS]{bms}}
#' @param mc_error The number of posterior draws will ensure that with 99% 
#' probability the bounds of the credible intervals will be within \eqn{\pm} 
#' \code{mc_error}\eqn{\times 4s_y}, that is, within 100\code{mc_error}% of the 
#' trimmed range of y.
#' @param seed Integer. Always set your seed!!!
#' @param ... Other arguments for \code{\link[BMS]{bms}}.
#' 
#' @return A list with the following elements:
#' \itemize{
#'  \item summary Tibble with point and interval estimates
#'  \item lm_b_fits A list of lm_b fits using zellner's g prior for
#'  all the top models from \code{\link[BMS]{bms}}
#'  \item hyperparameters A named list with the user-specified zellner's g value.
#'  \item posterior_draws matrix of posterior draws of the regression parameters, 
#'  marginalizing out the model
#' }
#' 
#' @examples
#' # Create data
#' set.seed(2025)
#' N = 500
#' test_data = 
#'   data.frame(x1 = rnorm(N),
#'              x2 = rnorm(N),
#'              x3 = letters[1:5],
#'              x4 = rnorm(N),
#'              x5 = rnorm(N),
#'              x6 = rnorm(N),
#'              x7 = rnorm(N),
#'              x8 = rnorm(N),
#'              x9 = rnorm(N),
#'              x10 = rnorm(N))
#' test_data$outcome = 
#'   rnorm(N,-1 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) )
#' 
#' # Fit linear model using Bayesian model averaging
#' fit <-
#'   bma_inference(outcome ~ .,
#'                 test_data,
#'                 user.int = FALSE)
#' summary(fit)
#' coef(fit)
#' credint(fit)
#' plot(fit)
#' 
#' 
#' @export


bma_inference = function(formula,
                         data,
                         zellner_g = nrow(data),
                         CI_level = 0.95,
                         ROPE,
                         mcmc_draws = 1e4,
                         mc_error = 0.001,
                         seed = 1,
                         ...){
  
  alpha = 1.0 - CI_level
  
  
  # Use BMS package to get top (a posteriori) models
  if(missing(data)){
    m = 
      model.frame(formula)
  }else{
    m = 
      model.frame(formula,data)
  }
  
  X = model.matrix(formula,m)
  X.data = 
    cbind(
      model.response(m),
      X[,-1]
    )
  colnames(X.data)[1] = all.vars(formula)[1]
  
  bms_fit = 
    BMS::bms(X.data,
             g = zellner_g,
             iter = mcmc_draws,
             ...)
  
  
  # Get the posterior probabilities of the models
  model_post_probs = 
    bms_fit$topmod$ncount()
  model_post_probs =
    model_post_probs / sum(model_post_probs)
  
  
  # Get initial set of posterior draws
  ## Get number of posterior samples per model
  mc_draws_by_model = 
    round(500 * model_post_probs)
  
  ## Extract a matrix showing us which variables are used 
  #   in each of our top-most models
  var_inclusion = 
    bms_fit$topmod$bool_binary()
  if(any(mc_draws_by_model == 0)) var_inclusion = var_inclusion[,-which(mc_draws_by_model == 0)]
  
  ## Fit top models
  X.data = 
    as.data.frame(X.data) |> 
    janitor::clean_names()
  full_fits = 
    future.apply::future_lapply(1:ncol(var_inclusion),
                                function(i){
                                  suppressWarnings(suppressPackageStartupMessages(library(bayesics)))
                                  suppressMessages(
                                    lm_b(paste0(colnames(X.data)[1], " ~ ", 
                                                paste(colnames(X.data)[-1][as.logical(var_inclusion[,i])],
                                                      collapse = " + ")) |> 
                                           as.formula(),
                                         data = X.data,
                                         prior = "zellner",
                                         zellner_g = zellner_g)
                                  )
                                },
                                future.seed = seed)
  
  ## Get posterior samples
  post_samples = 
    future.apply::future_lapply(1:length(full_fits),
                                function(i){
                                  suppressWarnings(suppressPackageStartupMessages(library(bayesics)))
                                  samples = 
                                    get_posterior_draws(full_fits[[i]],
                                                        n_draws = mc_draws_by_model[i]) |> 
                                    tibble::as_tibble()
                                  if(ncol(samples) < ncol(X.data)){ # Re dimension: Yes, X.data includes y, but the samples also ought to include s2
                                    for(j in setdiff(colnames(X.data)[-1],
                                                     colnames(samples))) samples[[j]] = 0.0
                                  }
                                  return(samples)
                                },
                                future.seed = seed)
  post_samples = 
    do.call(bind_rows,post_samples)
  
  post_samples =
    post_samples |> 
    na.omit()
  
  post_samples = 
    post_samples |> 
    dplyr::relocate(all_of(c(colnames(X.data)[-1],"s2")),
                    .after = "(Intercept)")
  
  ## Find the number of final draws necessary
  epsilon = mc_error * 4 * sd(model.response(m))
  post_samples = as.matrix(post_samples)
  vars_to_consider = 
    which(
      apply(post_samples[,-ncol(post_samples)],2,
            function(x) mean(dplyr::near(x,0.0))) <= 0.9
    )
  fhats = 
    future.apply::future_lapply(vars_to_consider,
                                function(i){
                                  stats::density(unlist(post_samples[,i]))
                                })
  mc_draws = 
    future.apply::future_sapply(vars_to_consider,
                                function(i){
                                  0.5 * alpha * (1.0 - 0.5 * alpha) *
                                    (
                                      qnorm(0.5 * (1.0 - 0.99)) / 
                                        epsilon /
                                        fhats[[i]]$y[which.min(abs(fhats[[i]]$x - 
                                                                     quantile(post_samples[,i], 0.5 * alpha)))]
                                    )^2
                                }) |> 
    max() |> 
    round()
  
  
  
  # Get final posterior draws
  ## Get number of posterior samples per model
  mc_draws_by_model = 
    round(mc_draws * model_post_probs)
  
  ## Extract a matrix showing us which variables are used 
  #   in each of our top-most models
  var_inclusion = 
    bms_fit$topmod$bool_binary()
  if(any(mc_draws_by_model == 0)) var_inclusion = var_inclusion[,-which(mc_draws_by_model == 0)]
  
  ## Fit top models
  X.data = 
    as.data.frame(X.data) |> 
    janitor::clean_names()
  full_fits = 
    future.apply::future_lapply(1:ncol(var_inclusion),
                  function(i){
                    suppressWarnings(suppressPackageStartupMessages(library(bayesics)))
                    suppressMessages(
                      lm_b(paste0(colnames(X.data)[1], " ~ ", 
                                  paste(colnames(X.data)[-1][as.logical(var_inclusion[,i])],
                                        collapse = " + ")) |> 
                             as.formula(),
                           data = X.data,
                           prior = "zellner",
                           zellner_g = zellner_g)
                    )
                  },
                  future.seed = seed)
  
  ## Get posterior samples
  post_samples = 
    future.apply::future_lapply(1:length(full_fits),
                  function(i){
                    suppressWarnings(suppressPackageStartupMessages(library(bayesics)))
                    samples = 
                      get_posterior_draws(full_fits[[i]],
                                          n_draws = mc_draws_by_model[i]) |> 
                      tibble::as_tibble()
                    if(ncol(samples) < ncol(X.data)){ # Re dimension: Yes, X.data includes y, but the samples also ought to include s2
                      for(j in setdiff(colnames(X.data)[-1],
                                       colnames(samples))) samples[[j]] = 0.0
                    }
                    return(samples)
                  },
                  future.seed = seed)
  post_samples = 
    do.call(bind_rows,post_samples)
  
  post_samples =
    post_samples |> 
    na.omit()
  
  post_samples = 
    post_samples |> 
    dplyr::relocate(all_of(c(colnames(X.data)[-1],"s2")),
             .after = "(Intercept)")
  
  # Summarize results
  ## Get values for ROPE
  s_y = sd(X.data[,1])
  s_j = apply(X.data[,-1,drop = FALSE],2,sd)
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
    c(
      paste("(",-round(ROPE,3),",",round(ROPE,3),")",sep=""),
      "(NA,NA)")
  
  
  boundaries = 
    matrix(c(ROPE,NA),
           nrow = nrow(post_samples),
           ncol = ncol(post_samples),
           byrow = TRUE)
  
  
  ## Compile results
  results = 
    tibble::tibble(Variable = c(colnames(X),"s2"),
           `Post Mean` = colMeans(post_samples),
           Lower = 
             apply(post_samples,2,quantile,probs = alpha/2),
           Upper =
             apply(post_samples,2,quantile,probs = 1.0 - alpha/2),
           `Prob Dir` = 
             c(apply(post_samples[,-ncol(post_samples)],
                     2,
                     function(x) max(mean(x < 0),
                                     mean(x > 0))),
               NA),
           ROPE = 
             colMeans(
               (-boundaries < post_samples) & 
                 (boundaries > post_samples)
             ),
           `ROPE bounds` = ROPE_bounds
    )
  
  # Get fitted values and NOT faux residuals
  fitted = 
    drop(
      X %*% results$`Post Mean`[-nrow(results)]
    )
  # residuals = 
  #   drop(
  #     X.data[,1] - fitted
  #   ) # This might mislead folks.  The data are NOT normally distributed.
  
  return_object = 
    list(summary = results,
         lm_b_fits = full_fits,
         bms_fit = bms_fit,
         hyperparms = list(zellner_g = zellner_g),
         posterior_draws = post_samples,
         fitted = fitted,
         formula = formula,
         data = data,
         CI_level = CI_level,
         terms = terms(m))
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
                   class = "lm_b_bma"))
}
