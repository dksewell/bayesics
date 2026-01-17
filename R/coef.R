#' @name coef
#' 
#' @title Coefficient extraction for bayesics objects
#' 
#' @param object bayesics object
#' @param ... optional arguments.
#' 
#' 
#' @returns vector of coefficients
#' 
#' @examples
#' set.seed(2025)
#' N = 500
#' test_data = 
#'   data.frame(x1 = rep(letters[1:5],N/5))
#' test_data$outcome = 
#'   rnorm(N,-1 + 2 * (test_data$x1 %in% c("d","e")) )
#' 
#' # Fit 1-way ANOVA model
#' fit1 <-
#'   aov_b(outcome ~ x1,
#'         test_data,
#'         prior_mean_mu = 2,
#'         prior_mean_nu = 0.5,
#'         prior_var_shape = 0.01,
#'         prior_var_rate = 0.01)
#' coef(fit1)
#' 

#' @rdname coef
#' @export
coef.lm_b = function(object, ...){
  object$summary$`Post Mean`
}

#' @rdname coef
#' @export
coef.aov_b = function(object, ...){
  object$posterior_parameters$mu_g
}

#' @rdname coef
#' @export
coef.np_glm_b = function(object, ...){
  ret = object$summary$`Post Mean`
  names(ret) = object$summary$Variable
  ret
}

#' @rdname coef
#' @export
coef.glm_b = function(object, ...){
  ret = object$summary$`Post Mean`[1:(nrow(object$summary) -
                                        (object$family$family == "negbinom") )]
  names(ret) = object$summary$Variable[1:(nrow(object$summary) -
                                            (object$family$family == "negbinom") )]
  ret
}

#' @rdname coef
#' @export
coef.lm_b_bma = function(object, ...){
  ret = object$summary$`Post Mean`
  names(ret) = object$summary$Variable
  ret[-length(ret)]
}
