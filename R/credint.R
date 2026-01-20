#' @name credint
#' 
#' @title Credible Intervals for Model Parameters
#' 
#' @description Computes credible intervals for one or more parameters in a 
#' fitted model.
#' 
#' @param object a fitted model object from \code{bayesics}
#' @param CI_level the credible level required
#' @param which character.  For \code{aov_b} only. Either "means" (for the 
#' group means) or "pairwise" (for pairwise difference in means).
#' @param ... Passed to methods.
#' 
#' @returns Matrix of credible intervals
#' 
#' @examples
#' \donttest{
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
#' credint(fit1)
#' }
#' 
#' 
#' @export
credint = function(object,...){
  UseMethod("credint")
}

#' @rdname credint
#' @exportS3Method credint lm_b 
credint.lm_b = function(object,
                        CI_level = 0.95,
                        ...){
  alpha = 1 - CI_level
  summ = object$summary[,c("Lower","Upper")]
  
  if(object$prior != "improper"){
    summ$Lower = 
      qlst(alpha/2,
           object$posterior_parameters$a_tilde,
           object$posterior_parameters$mu_tilde,
           sqrt(object$posterior_parameters$b_tilde/object$posterior_parameters$a_tilde * 
                  diag(qr.solve(object$posterior_parameters$V_tilde))))
    summ$Upper = 
      qlst(1.0 - alpha/2,
           object$posterior_parameters$a_tilde,
           object$posterior_parameters$mu_tilde,
           sqrt(object$posterior_parameters$b_tilde/object$posterior_parameters$a_tilde * 
                  diag(qr.solve(object$posterior_parameters$V_tilde))))
  }else{
    summ$Lower = 
      qlst(alpha/2,
           nrow(object$data) - length(object$posterior_parameters$mu_tilde),
           object$posterior_parameters$mu_tilde,
           sqrt(diag(object$posterior_parameters$Sigma)))
    summ$Upper = 
      qlst(1.0 - alpha/2,
           nrow(object$data) - length(object$posterior_parameters$mu_tilde),
           object$posterior_parameters$mu_tilde,
           sqrt(diag(object$posterior_parameters$Sigma)))
  }
  
  summ = 
    as.matrix(summ)
  colnames(summ) = 
    c(paste0(100 * 0.5 * alpha,"%"),
      paste0(100 * (1.0 - 0.5 * alpha),"%"))
  rownames(summ) = 
    object$summary$Variable
  
  return(summ)
}



#' @rdname credint
#' @exportS3Method credint aov_b
credint.aov_b = function(object,
                         CI_level = 0.95,
                         which = "means",
                         ...){
  
  alpha = 1 - CI_level
  
  which = 
    c("means","pairwise")[pmatch(tolower(which),
                                 c("means","pairwise"),
                                 duplicates.ok = FALSE)]
  
  if( !(which %in% c("means","pairwise")) )
    stop("The 'which' argument must be either 'means' or 'pairwise'")
  
  
  if(which == "means"){
    summ = object$summary[,c("Lower","Upper")]
    
    summ$Lower = 
      c(extraDistr::qlst(alpha/2, 
                         df = object$posterior_parameters$a_g,
                         mu = object$posterior_parameters$mu_g,
                         sigma = sqrt(object$posterior_parameters$b_g / object$posterior_parameters$nu_g / object$posterior_parameters$a_g)),
        extraDistr::qinvgamma(alpha/2, 
                              alpha = object$posterior_parameters$a_g/2, 
                              beta = object$posterior_parameters$b_g/2))
    summ$Upper = 
      c(extraDistr::qlst(1 - alpha/2, 
                         df = object$posterior_parameters$a_g,
                         mu = object$posterior_parameters$mu_g,
                         sigma = sqrt(object$posterior_parameters$b_g / object$posterior_parameters$nu_g / object$posterior_parameters$a_g)),
        extraDistr::qinvgamma(1 - alpha/2, 
                              alpha = object$posterior_parameters$a_g/2, 
                              beta = object$posterior_parameters$b_g/2))
    
    summ = 
      as.matrix(summ)
    colnames(summ) = 
      c(paste0(100 * 0.5 * alpha,"%"),
        paste0(100 * (1.0 - 0.5 * alpha),"%"))
    rownames(summ) = 
      object$summary$Variable
    
    return(summ)
    
  }else{#End: which == "means"; Start: which == "pairwise"
    pw_summ = 
      object$pairwise_summary[,c("Lower","Upper")] |> 
      as.data.frame()
    
    temp = 
      combn(1:length(levels(object$data$group)),2)
    for(i in 1:nrow(pw_summ)){
      pw_summ[i,c("Lower","Upper")] = 
        stats::quantile(object$posterior_draws[,temp[1,i]] - 
                   object$posterior_draws[,temp[2,i]],
                 probs = c(alpha/2, 
                           1 - alpha/2))
    }
    pw_summ = as.matrix(pw_summ)
    rownames(pw_summ) = 
      object$pairwise_summary$Comparison
    colnames(pw_summ) = 
      c(paste0(100 * 0.5 * alpha,"%"),
        paste0(100 * (1.0 - 0.5 * alpha),"%"))
    
    return(pw_summ)
  }#End: which == "pairwise"
  
}


#' @rdname credint
#' @exportS3Method credint glm_b
credint.glm_b = function(object,
                         CI_level = 0.95,
                         ...){
  alpha = 1 - CI_level
  summ = object$summary[,c("Lower","Upper")]
  if("posterior_covariance" %in% names(object)){
    summ$Lower = 
      qnorm(alpha / 2,
            object$summary$`Post Mean`,
            sd = sqrt(diag(object$posterior_covariance)))
    summ$Upper = 
      qnorm(1 - alpha / 2,
            object$summary$`Post Mean`,
            sd = sqrt(diag(object$posterior_covariance)))
  }else{
    # Get CI bounds
    CI_from_weighted_sample = function(x,w){
      w = cumsum(w[order(x)])
      x = x[order(x)]
      LB = max(which(w <= 0.5 * alpha))
      UB = min(which(w >= 1.0 - 0.5 * alpha))
      return(c(lower = x[LB],
               upper = x[UB]))
    }
    CI_bounds = 
      apply(object$proposal_draws,2,
            CI_from_weighted_sample,
            w = object$importance_sampling_weights)
    summ$Lower = 
      CI_bounds["lower",]
    summ$Upper = 
      CI_bounds["upper", ]
  }
  
  summ = 
    as.matrix(summ)
  colnames(summ) = 
    c(paste0(100 * 0.5 * alpha,"%"),
      paste0(100 * (1.0 - 0.5 * alpha),"%"))
  rownames(summ) = 
    object$summary$Variable
  
  return(summ)
}


#' @rdname credint
#' @exportS3Method credint np_glm_b
credint.np_glm_b = function(object,
                            CI_level = 0.95,
                            ...){
  alpha = 1 - CI_level
  summ = object$summary[,c("Lower","Upper")]
  if("posterior_covariance" %in% names(object)){
    summ$Lower = 
      qnorm(alpha / 2,
            object$summary$`Post Mean`,
            sd = sqrt(diag(as.matrix(object$posterior_covariance))))
    summ$Upper = 
      qnorm(1 - alpha / 2,
            object$summary$`Post Mean`,
            sd = sqrt(diag(as.matrix(object$posterior_covariance))))
  }else{
    summ$Lower = 
      object$posterior_draws |> 
      apply(2,stats::quantile,prob = alpha / 2)
    summ$Upper = 
      object$posterior_draws |> 
      apply(2,stats::quantile,prob = 1.0 - alpha / 2)
  }
  
  summ = 
    as.matrix(summ)
  colnames(summ) = 
    c(paste0(100 * 0.5 * alpha,"%"),
      paste0(100 * (1.0 - 0.5 * alpha),"%"))
  rownames(summ) = 
    object$summary$Variable
  
  return(summ)
}

#' @rdname credint
#' @exportS3Method credint lm_b_bma
credint.lm_b_bma = function(object,
                            CI_level = 0.95,
                            ...){
  alpha = 1 - CI_level
  summ = object$summary[,c("Lower","Upper")]
  summ$Lower = 
    apply(object$posterior_draws,2,stats::quantile,probs = alpha/2)
  summ$Upper =
    apply(object$posterior_draws,2,stats::quantile,probs = 1.0 - alpha/2)
  
  summ = 
    as.matrix(summ)
  colnames(summ) = 
    c(paste0(100 * 0.5 * alpha,"%"),
      paste0(100 * (1.0 - 0.5 * alpha),"%"))
  rownames(summ) = 
    object$summary$Variable
  
  return(summ)
}






















