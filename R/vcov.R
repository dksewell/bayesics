#' @name vcov
#' 
#' @title Calculate Posterior Variance-Covariance Matrix for a Bayesian Fitted Model Object
#' @description
#' Returns the posterior covariance matrix of the main parameters of a fitted \code{bayesics} object
#' 
#' @param object a fitted model object from \code{bayesics}.
#' @param ... Passed to methods.
#' 
#' @returns A matrix of the covariance matrix for the regression coefficients.  If the posterior 
#' is a multivariate t distribution (or consists of independent t's in the case of heteroscedastic 
#' 1-way ANOVA), the degrees of freedom are returned as the \code{df} attribute of the matrix.  Note 
#' that for \code{lm_b} and \code{aov_b} objects, this function already takes into account the 
#' uncertainty around the residual variance.
#' 
#' @rdname vcov
#' @method vcov aov_b
#' @export
vcov.aov_b = function(object,...){
  covmat = 
    diag(object$posterior_parameters$b_g / 
           object$posterior_parameters$a_g / 
           object$posterior_parameters$nu_g)
  attr(covmat,"df") = object$posterior_parameters$a_g
  
  return(covmat)
}

#' @rdname vcov
#' @method vcov lm_b
#' @export
vcov.lm_b = function(object,...){
  covmat = NULL
  try({
    covmat = 
      chol2inv(chol(object$posterior_parameters$V_tilde))
  }, silent=TRUE)
  if(is.null(covmat)){
    try({
      covmat = 
        qr.solve(object$posterior_parameters$V_tilde)
    }, silent=TRUE)
  }
  if(is.null(covmat)){
    try({
      covmat = 
        solve(object$posterior_parameters$V_tilde)
    }, silent=TRUE)
  }
  if(is.null(covmat)) stop("Hessian is not invertible.")
  attr(covmat,"df") = object$posterior_parameters$a_tilde
  
  return(covmat)
}

#' @rdname vcov
#' @method vcov glm_b
#' @export
vcov.glm_b = function(object,...){
  if("posterior_covariance" %in% names(object)){
    return(object$posterior_covariance)
  }else{
    first_moment =
      object$summary$`Post Mean`
    second_moment = 
      crossprod(object$proposal_draws,
                Diagonal(x = object$importance_sampling_weights) %*% object$proposal_draws)
    
    return(as.matrix(second_moment - tcrossprod(first_moment)))
  }
}



#' @rdname vcov
#' @method vcov np_glm_b
#' @export
vcov.np_glm_b = function(object,...){
  if("posterior_covariance" %in% names(object)){
    return(object$posterior_covariance)
  }else{
    return(cov(na.omit(object$posterior_draws)))
  }
}
