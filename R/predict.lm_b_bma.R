#' Predict method for bma model fits
#' 
#' 
#' @param object Object of class bma
#' @param newdata An optional data.frame in which to look for variables with which 
#' to predict. 
#' @param CI_level Posterior probability covered by credible interval
#' @param PI_level Posterior probability covered by prediction interval
#' @param seed integer. Always set your seed!!!
#' @param ... optional arguments.
#' 
#' @returns list.
#' \itemize{
#'  \item newdata tibble with estimate, prediction intervals, and credible intervals 
#' for the mean.
#'  \item posterior_draws
#'    \itemize{
#'      \item mean_of_ynew draws of \eqn{E(y)}, marginalizing out the model
#'      \item posterior draws of ynew
#'    }
#'  }
#'  
#'  @examples
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
#' predict(fit)
#' 
#' 
#' @exportS3Method predict lm_b_bma

predict.lm_b_bma = function(object,
                            newdata,
                            CI_level = 0.95,
                            PI_level = 0.95,
                            seed = 1,
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
  
  m = model.frame(delete.response(terms(object)),
                  data = newdata)
  
  X = model.matrix(delete.response(terms(object)),
                   data = newdata)
  N = nrow(X)
  p = ncol(X)
  
  # Get means of E(y|X)
  mu_draws = 
    tcrossprod(as.matrix(object$posterior_draws[,1:ncol(X)]),
               X)
  
  # Get draws of y_new
  y_draws = 
    mu_draws + 
    sqrt(object$posterior_draws$s2) * 
    matrix(rnorm(prod(dim(mu_draws))),
           nrow(mu_draws),
           ncol(mu_draws))
  
  # Compile results
  newdata =
    newdata |> 
    dplyr::mutate(`Post Mean` = 
                    colMeans(mu_draws),
                  CI_lower = 
                    apply(mu_draws,2,quantile,probs = 0.5 * alpha_ci),
                  CI_upper = 
                    apply(mu_draws,2,quantile,probs = 1.0 - 0.5 * alpha_ci),
                  PI_lower = 
                    apply(y_draws,2,quantile,probs = 0.5 * alpha_pi),
                  PI_upper = 
                    apply(y_draws,2,quantile,probs = 1.0 - 0.5 * alpha_pi))
  
  return(list(newdata = newdata,
              posterior_draws = 
                list(mean_of_ynew = mu_draws,
                     ynew = y_draws)))
}
