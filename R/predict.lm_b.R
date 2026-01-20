#' Predict method for lm_b model fits
#' 
#' 
#' @param object Object of class lm_b
#' @param newdata An optional data.frame in which to look for variables with which 
#' to predict. 
#' @param CI_level Posterior probability covered by credible interval
#' @param PI_level Posterior probability covered by prediction interval
#' @param n_draws If desired, the number of posterior samples drawn.
#' @param ... optional arguments.
#' 
#' @returns tibble with estimate (posterior mean), prediction intervals, and credible intervals 
#' for the mean.
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
#' predict(fit1)
#' }
#' 
#' @exportS3Method predict lm_b

predict.lm_b = function(object,
                        newdata,
                        CI_level = 0.95,
                        PI_level = 0.95,
                        n_draws = 0,
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
  
  V_eig = eigen(object$posterior_parameters$V_tilde)
  XULneghalf = X %*% V_eig$vectors %*% diag(1 / sqrt(V_eig$values))
  X_VInverse_X = 
    apply(XULneghalf,1,crossprod)
  
  temporary_mean = drop(X %*% object$posterior_parameters$mu_tilde)
  
  newdata =
    newdata |> 
    dplyr::mutate(`Post Mean` = temporary_mean) |> 
    dplyr::mutate(PI_lower = 
                    extraDistr::qlst(alpha_pi/2.0,
                                     df = object$posterior_parameters$a_tilde,
                                     mu = .data$`Post Mean`,
                                     sigma = sqrt(object$posterior_parameters$b_tilde / 
                                                    object$posterior_parameters$a_tilde * 
                                                    (X_VInverse_X + 1.0) ) ),
                  PI_upper = 
                    extraDistr::qlst(1.0 - alpha_pi/2.0,
                                     df = object$posterior_parameters$a_tilde,
                                     mu = .data$`Post Mean`,
                                     sigma = sqrt(object$posterior_parameters$b_tilde / 
                                                    object$posterior_parameters$a_tilde * 
                                                    (X_VInverse_X + 1.0) ) ),
                  CI_lower = 
                    extraDistr::qlst(alpha_ci/2.0,
                                     df = object$posterior_parameters$a_tilde,
                                     mu = .data$`Post Mean`,
                                     sigma = sqrt(object$posterior_parameters$b_tilde / 
                                                    object$posterior_parameters$a_tilde * 
                                                    (X_VInverse_X) ) ),
                  CI_upper = 
                    extraDistr::qlst(1.0 - alpha_ci/2.0,
                                     df = object$posterior_parameters$a_tilde,
                                     mu = .data$`Post Mean`,
                                     sigma = sqrt(object$posterior_parameters$b_tilde / 
                                                    object$posterior_parameters$a_tilde * 
                                                    (X_VInverse_X) ) )
    )
  
  # Get new draws
  if(n_draws > 0){
    new_draws = matrix(0.0,
                       nrow(newdata),
                       n_draws,
                       dimnames = list(NULL,
                                       paste("y_new",1:n_draws,sep="")))
    for(it in 1:n_draws){
      new_draws[,it] = 
        extraDistr::rlst(nrow(newdata),
                         df = object$posterior_parameters$a_tilde,
                         mu = newdata$`Post Mean`,
                         sigma = sqrt(object$posterior_parameters$b_tilde / 
                                        object$posterior_parameters$a_tilde * 
                                        (X_VInverse_X + 1.0) ) )
    }
    
    newdata =
      newdata |> 
      dplyr::bind_cols(new_draws |> 
                         tibble::as_tibble())
  }
  
  return(newdata)
}
