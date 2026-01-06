#' @name bayes_factors
#'
#' 
#' @title Bayes factors for lm_b, glm_b, and survfit_b
#' 
#' @description
#' Bayes factors for Bayesian regression objects using the Savage-Dickey ratio
#' 
#' 
#' @details
#' 
#' Bayes factors are given in terms of favoring the two-tailed alternative hypothesis 
#' vs. the null hypothesis that the regression coefficient equals zero.   
#' Currently implemented for \code{lm_b} or \code{glm_b} objects.  Note 
#' that for \code{glm_b} objects, if importance sampling was used, 
#' the model will be refit using fixed form variational Bayes to get 
#' the multivariate posterior density.
#' 
#' Interpretation is taken from Kass and Raftery.
#' 
#' @param object lm_b or glm_b object
#' @param by character.  Either "coefficient" or "variable".  If the former, 
#' Bayes factors will be computed for each regression coefficient separately. 
#' If the latter, Bayes factors will be computed for each covariate separately.
#' 
#' @references 
#' 
#' James M. Dickey. "The Weighted Likelihood Ratio, Linear Hypotheses on Normal Location Parameters." Ann. Math. Statist. 42 (1) 204 - 223, February, 1971. https://doi.org/10.1214/aoms/1177693507
#' 
#' Kass, R. E., & Raftery, A. E. (1995). Bayes Factors. Journal of the American Statistical Association, 90(430), 773–795.
#' 
#' 
#' @export

#' @export
bayes_factors = function(object,...){
  UseMethod("bayes_factors")
}


#' @rdname bayes_factors
#' @exportS3Method bayes_factors lm_b
bayes_factors.lm_b = function(object,
                              by = "coefficient",
                              ...){
  if(is.na(object$hyperparameters[1])) 
    stop("Cannot compute Bayes factors with an improper prior.")
  
  by = match.arg(tolower(by),
                 c("coefficient",
                   "variable"))
  
  if(by == "coefficient"){
    log_numerators = 
      extraDistr::dlst(0.0,
                       df = object$hyperparameters$a,
                       mu = object$hyperparameters$mu,
                       sigma = 
                         sqrt(object$hyperparameters$b / 
                                object$hyperparameters$a * 
                                diag(qr.solve(object$hyperparameters$V))),
                       log = TRUE)
    
    log_denominators = 
      extraDistr::dlst(0.0,
                       df = object$posterior_parameters$a,
                       mu = object$posterior_parameters$mu,
                       sigma = 
                         sqrt(object$posterior_parameters$b / 
                                object$posterior_parameters$a * 
                                diag(qr.solve(object$posterior_parameters$V))),
                       log = TRUE)
    
    results = 
      tibble::tibble(Variable = colnames(object$posterior_parameters$V_tilde),
                     `BF favoring alternative` = 
                       exp(log_numerators - log_denominators))
    
  }else{#End: by = "coefficient"
    
    # Find which coefficients correspond to which covariates
    ## Get term labels
    tl = 
      attr(terms(object),"term.labels")
    ## See how these term labels map to the coefficients of the model
    coef_to_tl = 
      attr(model.matrix(terms(object),
                        data = object$data),
           "assign")
    ## See which covariates appear in which terms/coefficients
    covars_per_tl = 
      lapply(tl,
             function(x){
               all.vars(as.formula(paste0("~",x)))
             })
    ## Extract the unique covariates
    unique_covars = 
      unique(unlist(covars_per_tl))
    ## For each covariate, find the corresponding terms
    tl_to_covars = 
      lapply(1:length(unique_covars),
             function(i){
               sapply(covars_per_tl,function(x) unique_covars[i] %in% x) |> 
                 which()
             })
    ## Finally, for each covariate, find the corresponding regression coefficients
    coefs_to_covars = 
      lapply(tl_to_covars,
             function(x){
               which(sapply(coef_to_tl,function(z) z %in% x))
             })
    
    # Now compute Bayes factors via SD ratio
    log_numerators = 
      log_denominators = 
      numeric(length(unique_covars))
    beta_covmat = 
      vcov(object)
    prior_beta_covmat = 
      object$hyperparameters$b / 
      object$hyperparameters$a * 
      qr.solve(object$hyperparameters$V)
    for(j in 1:length(unique_covars)){
      j_index = coefs_to_covars[[j]]
      log_numerators[j] = 
        mvtnorm::dmvt(numeric(length(j_index)),
                      df = object$hyperparameters$a,
                      delta = object$hyperparameters$mu[j_index],
                      sigma = 
                        prior_beta_covmat[j_index,j_index,drop=FALSE],
                      log = TRUE)
      log_denominators[j] = 
        mvtnorm::dmvt(numeric(length(j_index)),
                      df = object$posterior_parameters$a,
                      delta = object$posterior_parameters$mu[j_index],
                      sigma = 
                        beta_covmat[j_index,j_index,drop=FALSE],
                      log = TRUE)
    }
    
    
    
    results = 
      tibble::tibble(Variable = unique_covars,
                     `BF favoring alternative` = 
                       exp(log_numerators - log_denominators))
    
  }
  
  results = 
    results |>
    dplyr::mutate(bf_max = pmax(.data$`BF favoring alternative`,
                                1.0 / .data$`BF favoring alternative`)) |> 
    dplyr::mutate(Interpretation = 
                    ifelse(.data$bf_max <= 3.2,
                           "Not worth more than a bare mention",
                           ifelse(.data$bf_max <= 10,
                                  "Substantial",
                                  ifelse(.data$bf_max <= 100,
                                         "Strong",
                                         "Decisive"))) |> 
                    paste(ifelse(.data$`BF favoring alternative` > 1,
                                 " (in favor of keeping in the model)",
                                 " (in favor of exluding from the model")))
  results$bf_max = NULL
  
  return(results)
}




#' @rdname bayes_factors
#' @exportS3Method bayes_factors glm_b
bayes_factors.glm_b = function(object,
                               by = "coefficient",
                               ...){
  if(is.na(object$hyperparameters[1])) 
    stop("Cannot compute Bayes factors with an improper prior.")
  
  by = match.arg(tolower(by),
                 c("coefficient",
                   "variable"))
  
  # If object was fit via IS, use VB instead.
  if(object$algorithm == "IS"){
    suppressMessages({
      object <- 
        glm_b(formula = object$formula,
              data = object$data,
              family = object$family,
              trials = object$trials,
              prior_beta_mean = object$hyperparameters$prior_beta_mean,
              prior_beta_precision = object$hyperparameters$prior_beta_precision,
              algorithm = "VB")
    })
  }
  
  # Get prior covariance
  prior_cov = NULL
  try({
    prior_cov = 
      chol2inv(chol(object$hyperparameters$prior_beta_precision))
  },silent = TRUE)
  if(is.null(prior_cov)){
    try({
      prior_cov = 
        qr.solve(object$hyperparameters$prior_beta_precision)
    },silent = TRUE)
  }
  if(is.null(prior_cov)){
    try({
      prior_cov = 
        solve(object$hyperparameters$prior_beta_precision)
    },silent = TRUE)
  }
  
  # Get posterior covariance
  post_cov = vcov(object)
  
  # Get dimension
  p = 
    nrow(post_cov) - (object$family$family == "negbinom")
  
  if(by == "coefficient"){
    
    log_numerators = 
      dnorm(0.0,
            object$hyperparameters$prior_beta_mean[1:p],
            sqrt(diag(prior_cov)[1:p]),
            log = TRUE)
    
    log_denominators = 
      dnorm(0.0,
            object$summary$`Post Mean`[1:p],
            sqrt(diag(object$posterior_covariance)[1:p]),
            log = TRUE)
    
    
    results = 
      tibble::tibble(Variable = object$summary$Variable[1:p],
                     `BF favoring alternative` = 
                       exp(log_numerators - log_denominators))
    
  }else{#End: by = "coefficients"
    # Find which coefficients correspond to which covariates
    ## Get term labels
    tl = 
      attr(terms(object),"term.labels")
    ## See how these term labels map to the coefficients of the model
    coef_to_tl = 
      attr(model.matrix(terms(object),
                        data = object$data),
           "assign")
    ## See which covariates appear in which terms/coefficients
    covars_per_tl = 
      lapply(tl,
             function(x){
               all.vars(as.formula(paste0("~",x)))
             })
    ## Extract the unique covariates
    unique_covars = 
      unique(unlist(covars_per_tl))
    ## For each covariate, find the corresponding terms
    tl_to_covars = 
      lapply(1:length(unique_covars),
             function(i){
               sapply(covars_per_tl,function(x) unique_covars[i] %in% x) |> 
                 which()
             })
    ## Finally, for each covariate, find the corresponding regression coefficients
    coefs_to_covars = 
      lapply(tl_to_covars,
             function(x){
               which(sapply(coef_to_tl,function(z) z %in% x))
             })
    
    # Now compute Bayes factors via SD ratio
    log_numerators = 
      log_denominators = 
      numeric(length(unique_covars))
    for(j in 1:length(unique_covars)){
      j_index = coefs_to_covars[[j]]
      log_numerators[j] = 
        mvtnorm::dmvnorm(numeric(length(j_index)),
                         mean = object$hyperparameters$prior_beta_mean[j_index],
                         sigma = 
                           prior_cov[j_index,j_index,drop=FALSE],
                         log = TRUE)
      log_denominators[j] = 
        mvtnorm::dmvnorm(numeric(length(j_index)),
                         mean = coef(object)[j_index],
                         sigma = 
                           post_cov[j_index,j_index,drop=FALSE],
                         log = TRUE)
    }
    
    
    results = 
      tibble::tibble(Variable = unique_covars,
                     `BF favoring alternative` = 
                       exp(log_numerators - log_denominators))
    
  }
  
  
  results = 
    results |>
    dplyr::mutate(bf_max = pmax(.data$`BF favoring alternative`,
                                1.0 / .data$`BF favoring alternative`)) |> 
    dplyr::mutate(Interpretation = 
                    ifelse(.data$bf_max <= 3.2,
                           "Not worth more than a bare mention",
                           ifelse(.data$bf_max <= 10,
                                  "Substantial",
                                  ifelse(.data$bf_max <= 100,
                                         "Strong",
                                         "Decisive"))) |> 
                    paste(ifelse(.data$`BF favoring alternative` > 1,
                                 " (in favor of keeping in the model)",
                                 " (in favor of exluding from the model")))
  results$bf_max = NULL
  
  return(results)
}




#' @rdname bayes_factors
#' @exportS3Method bayes_factors survfit_b
bayes_factors.survfit_b = function(object,
                                   object2,
                                   ...){
  
  # Make quick checks
  if(object$single_group_analysis)
    object$group_names = character(0)
  
  if(object2$single_group_analysis)
    object2$group_names = character(0)
  
  if(isTRUE(all.equal(object$group_names, object2$group_names)))
    stop("Models must be different.")
  
  if(!isTRUE(all.equal(model.response(object$data),
                       model.response(object2$data))))
    stop("Data must be the same to compare via Bayes factors")
  
  G1 = length(object$group_names)
  G2 = length(object2$group_names)
  
  BF = 
    exp( 
      object2$marginal_likelihood - 
        object$marginal_likelihood
    )
  
  bf_max = 
    pmax(BF, 1.0 / BF)
  
  favor_or_against = 
    ifelse(BF > 1,
           " (in favor of the ",
           " (against the ")
  
  Interpretation = 
    ifelse(bf_max <= 3.2,
           "Not worth more than a bare mention",
           ifelse(bf_max <= 10,
                  "Substantial",
                  ifelse(bf_max <= 100,
                         "Strong",
                         "Decisive"))) |> 
    paste(ifelse( G2 > G1,
                  paste0(favor_or_against,"larger model)"),
                  ifelse( G2 < G1,
                          paste0(favor_or_against,"smaller model)"),
                          ifelse(G2 == G1,
                                 paste0(favor_or_against,"the second model)"),
                          ))))
  
  cat("\n----------\n\nSemi-parametric survival curve fitting using Bayesian techniques\n")
  cat("\n----------\n\n")
  cat(paste0("The Bayes factor equaled ",
             format(signif(BF, 3)),
             ".\nInterpretation: ",
             Interpretation))
  cat("\n----------\n\n")
  
  invisible(BF)
}



