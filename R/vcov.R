#' @name vcov
#' 
#' @title Calculate Posterior Variance-Covariance Matrix for a Bayesian Fitted Model Object
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
#' vcov(fit1)
#' }
#' 
#' @export

#' @rdname vcov
#' @method vcov lm_b
#' @export
vcov.lm_b = function(object,...){
  
  if("posterior_covariance" %in% names(object)){ # Handles lm, glm\IS, np_glm\bootstrapping
    
    covmat = 
      object$posterior_covariance
    attr(covmat,"df") = object$df
    
  }else{
    
    if("importance_sampling_weights" %in% names(object)){ # Handles glm IS
      
      covmat = 
        crossprod(object$proposal_draws,
                  object$importance_sampling_weights * object$proposal_draws) - 
        tcrossprod(colSums(object$importance_sampling_weights * object$proposal_draws))
      attr(covmat,"df") = NA
      
    }else{
      if("posterior_draws" %in% names(object)){ # Handles np_glm bootstrapping, bma_inference
        
        covmat = 
          cov(object$posterior_draws)
        attr(covmat,"df") = NA
        
      }#End: posterior_draws if
    }#End: importance_sampling_weights ifelse
  }#End: posterior_covariance ifelse
  
  
  return(covmat)
}