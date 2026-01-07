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
#' @param mc_draws Integer. Number of draws in Monte Carlo integration. Passed into \code{\link[BMS]{bms}}
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
#' 
#' @export


bma_inference = function(formula,
                         data,
                         zellner_g = nrow(data),
                         CI_level = 0.95,
                         mc_draws = 1e4,
                         seed = 1,
                         ...){
  
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
             iter = mc_draws,
             ...)
  
  
  # Get the posterior probabilities of the models
  model_post_probs = 
    bms_fit$topmod$ncount()
  model_post_probs =
    model_post_probs / sum(model_post_probs)
  
  # Get number of posterior samples per model
  mc_draws_by_model = 
    round(mc_draws * model_post_probs)
  
  # Extract a matrix showing us which variables are used 
  # in each of our top-most models
  var_inclusion = 
    bms_fit$topmod$bool_binary()
  if(any(mc_draws_by_model == 0)) var_inclusion = var_inclusion[,-which(mc_draws_by_model == 0)]
  
  # Fit top models
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
  
  # Get posterior samples
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
  alpha = 1.0 - CI_level
  ## Get values for ROPE
  s_y = sd(X.data[,1])
  s_j = apply(X.data[,-1,drop = FALSE],2,sd)
  boundaries = 
    matrix(0.05 * s_y / c(NA,s_j,NA),
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
           ROPE = 
             colMeans(
               (-boundaries < post_samples) & 
                 (boundaries > post_samples)
             ),
           `Prob Dir` = 
             c(apply(post_samples[,-ncol(post_samples)],
                     2,
                     function(x) max(mean(x < 0),
                                     mean(x > 0))),
               NA)
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
