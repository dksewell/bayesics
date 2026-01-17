#' @name summary
#' 
#' @title Summary functions for bayesics objects
#' 
#' @param object bayesics object
#' @param CI_level Posterior probability covered by credible interval
#' @param ... optional arguments.
#' 
#' @returns tibble with summary values
#' 
#' @examples
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
#' summary(fit1)
#' 
#' @export

#' @rdname summary
#' @method summary lm_b 
#' @export
summary.lm_b = function(object,
                        CI_level = 0.95,
                        ...){
  alpha = 1 - CI_level
  summ = object$summary
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
  
  if(object$prior != "improper"){
    BF = bayes_factors(object)
    summ  = 
      summ |> 
      left_join(BF,
                by = "Variable")
  }
  
  summ
}

#' @rdname summary
#' @method summary aov_b 
#' @export
summary.aov_b = function(object,
                         CI_level = 0.95,
                         ...){
  alpha = 1 - CI_level
  summ = object$summary
  pw_summ = 
    object$pairwise_summary |> 
    as.data.frame()
  
  if("BF_for_different_vs_same_means" %in% names(object)){
    
    bf_max = 
      max(object$BF_for_different_vs_same_means, 
          1.0 / object$BF_for_different_vs_same_means)
    
    cat("\n---\n") 
    cat(paste0(
      "Bayes factor in favor of the full vs. null model: ",
      format(signif(object$BF_for_different_vs_same_means, 3), 
             scientific = 
               (object$BF_for_different_vs_same_means > 1e3) | 
               (object$BF_for_different_vs_same_means < 1e-3)),
      ";\n      =>Level of evidence: ", 
      ifelse(bf_max <= 3.2,
             "Not worth more than a bare mention",
             ifelse(bf_max <= 10,
                    "Substantial",
                    ifelse(bf_max <= 100,
                           "Strong",
                           "Decisive")))
    )
    )
    
  }
  cat("\n\n\n\n--- Summary of factor level means ---\n")
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
  print(summ)
  
  
  cat("\n\n\n\n--- Summary of pairwise differences ---\n")
  temp = 
    combn(1:length(levels(object$data$group)),2)
  for(i in 1:nrow(pw_summ)){
    pw_summ[i,c("Lower","Upper")] = 
      quantile(object$posterior_draws[,temp[1,i]] - 
                 object$posterior_draws[,temp[2,i]],
               probs = c(alpha/2, 
                         1 - alpha/2))
  }
  pw_summ = as_tibble(pw_summ)
  print(pw_summ)
  cat("\n\n   *Note: EPR (Exceedence in Pairs Rate) for a Comparison of g-h = Pr(Y_(gi) > Y_(hi)|parameters) ")
  
  if(is.null(object$contrasts)){
    invisible(list(summary = object$summary,
                   pairwise = object$pairwise_summary))
  }else{
    
    cat("\n\n\n\n--- Summary of Contrasts ---\n")
    csumm = object$contrasts$summary
    contrast_draws = 
      tcrossprod(object$posterior_draws[,grep("mean_",colnames(object$posterior_draws))],
                 object$contrasts$L)
    csumm$Lower = 
      apply(contrast_draws,2,quantile,probs = alpha/2)
    csumm$Upper = 
      apply(contrast_draws,2,quantile,probs = 1 - alpha/2)
    
    print(csumm)
    
    invisible(list(summary = object$summary,
                   pairwise = object$pairwise_summary,
                   contrasts = 
                     list(L = object$contrasts$L,
                          summary = csumm)))
  }
}

#' @param interpretable_scale ADD description!
#' @rdname summary
#' @method summary np_glm_b 
#' @export
summary.np_glm_b = function(object,
                            CI_level = 0.95,
                            interpretable_scale = TRUE,
                            ...){
  alpha = 1 - CI_level
  summ = object$summary
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
      apply(2,quantile,prob = alpha / 2)
    summ$Upper = 
      object$posterior_draws |> 
      apply(2,quantile,prob = 1.0 - alpha / 2)
  }
  
  # Exponentiate
  if(( (object$family$family == "binomial") & 
       (object$family$link != "logit") ) | 
     ( (object$family$family == "poisson") & 
       (object$family$link != "log") ) | 
     (object$family$family == "gaussian") ){
    interpretable_scale = FALSE
  }
  if(interpretable_scale){
    if("ROPE bounds" %in% colnames(summ)){
      rbounds = 
        sapply(summ$`ROPE bounds`,
               function(x){
                 scan(text = 
                        gsub("[()]",
                             "",
                             x),
                      what = numeric(),
                      sep = ",",
                      quiet = TRUE)
               }) |> 
        t() |> 
        exp()
      for(i in 1:nrow(summ)){
        summ$`ROPE bounds`[i] = 
          paste0("(",
                 round(rbounds[i,1],3),
                 ",",
                 round(rbounds[i,2],3),
                 ")")
      }
    }
    
    if("log(phi)" %in% summ$Variable)
      summ$Variable[which(summ$Variable == "log(phi)")] = "phi"
    
    paste0("\n----------\n\nValues given in terms of ",
           ifelse(object$family$family == "binomial",
                  "odds ratios",
                  "rate ratios")
    ) |> 
      cat()
    cat("\n\n----------\n\n")
    summ = summ[-1,]
    summ[,c("Post Mean","Lower","Upper")] =
      summ[,c("Post Mean","Lower","Upper")] |> 
      exp()
  }
  
  summ
}

#' @rdname summary
#' @method summary lm_b_bma 
#' @export
summary.lm_b_bma = function(object,
                            CI_level = 0.95,
                            ...){
  alpha = 1 - CI_level
  summ = object$summary
  summ$Lower = 
    apply(object$posterior_draws,2,quantile,probs = alpha/2)
  summ$Upper =
    apply(object$posterior_draws,2,quantile,probs = 1.0 - alpha/2)
  
  
  
  summ
}


#' @rdname summary
#' @method summary glm_b 
#' @export
summary.glm_b = function(object,
                         CI_level = 0.95,
                         interpretable_scale = TRUE,
                         ...){
  alpha = 1 - CI_level
  summ = object$summary
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
  
  # Exponentiate
  if(( (object$family$family == "binomial") & 
       (object$family$link != "logit") ) | 
     ( (object$family$family == "poisson") & 
       (object$family$link != "log") )){
    interpretable_scale = FALSE
  }
  
  if(interpretable_scale){
    paste0("\n----------\n\nValues given in terms of ",
           ifelse(object$family$family == "binomial",
                  "odds ratios",
                  "rate ratios")
    ) |> 
      cat()
    cat("\n\n----------\n\n")
    summ = summ[-1,]
    summ[,c("Post Mean","Lower","Upper")] =
      summ[,c("Post Mean","Lower","Upper")] |> 
      exp()
    summ[,"ROPE bounds"] = 
      paste("(",
            round(exp(-object$ROPE[-1]),3),
            ",",
            round(exp(object$ROPE[-1]),3),
            ")",
            sep="")
    if(object$family$family == "negbinom"){
      summ$Variable[nrow(summ)] = "phi"
    }
  }
  
  
  
  if(object$prior != "improper"){
    BF = bayes_factors(object)
    if(object$family$family == "negbinom"){
      BF = 
        dplyr::bind_rows(BF,
                         tibble::tibble(Variable = 
                                          ifelse(interpretable_scale,
                                                 "phi",
                                                 "log(phi)"),
                                        `BF favoring alternative` = NA,
                                        Interpretation = NA))
    }
    summ  = 
      summ |> 
      dplyr::left_join(BF,
                       by = "Variable")
  }
  
  summ
}

#' @rdname summary
#' @method summary mediate_b 
#' @export
summary.mediate_b = function(object,
                             CI_level = 0.95,
                             ...){
  alpha_ci = 1 - CI_level
  summ = object$summary
  nr = nrow(summ)
  
  # Simple case
  if(nr == 4){
    summ$Lower = 
      c(quantile(object$posterior_draws$ACME,
                 probs = 0.5 * alpha_ci),
        quantile(object$posterior_draws$ADE,
                 probs = 0.5 * alpha_ci),
        quantile(object$posterior_draws$`Total Effect`,
                 probs = 0.5 * alpha_ci),
        quantile(object$posterior_draws$ACME / 
                   object$posterior_draws$`Total Effect`,
                 probs = 0.5 * alpha_ci))
    summ$Upper =
      c(quantile(object$posterior_draws$ACME,
                 probs = 1.0 - 0.5 * alpha_ci),
        quantile(object$posterior_draws$ADE,
                 probs = 1.0 - 0.5 * alpha_ci),
        quantile(object$posterior_draws$`Total Effect`,
                 probs = 1.0 - 0.5 * alpha_ci),
        quantile(object$posterior_draws$ACME / 
                   object$posterior_draws$`Total Effect`,
                 probs = 1.0 - 0.5 * alpha_ci))
  }else{#End: simple case
    # Complex case
    summ$Lower = 
      c(quantile(object$posterior_draws$ACME_control,0.5 * alpha_ci),
        quantile(object$posterior_draws$ACME_treat,0.5 * alpha_ci),
        quantile(object$posterior_draws$ADE_control,0.5 * alpha_ci),
        quantile(object$posterior_draws$ADE_treat,0.5 * alpha_ci),
        quantile(object$posterior_draws$Tot_Eff,0.5 * alpha_ci),
        0.5 * quantile(object$posterior_draws$ACME_control + 
                         object$posterior_draws$ACME_treat,0.5 * alpha_ci),
        0.5 * quantile(object$posterior_draws$ADE_control + 
                         object$posterior_draws$ADE_treat,0.5 * alpha_ci),
        quantile( (object$posterior_draws$ACME_control + 
                     object$posterior_draws$ACME_treat) / 
                    (object$posterior_draws$ACME_control + 
                       object$posterior_draws$ACME_treat + 
                       object$posterior_draws$ADE_control + 
                       object$posterior_draws$ADE_treat), 0.5 * alpha_ci )
      )
    summ$Upper = 
      c(quantile(object$posterior_draws$ACME_control,1.0 - 0.5 * alpha_ci),
        quantile(object$posterior_draws$ACME_treat,1.0 - 0.5 * alpha_ci),
        quantile(object$posterior_draws$ADE_control,1.0 - 0.5 * alpha_ci),
        quantile(object$posterior_draws$ADE_treat,1.0 - 0.5 * alpha_ci),
        quantile(object$posterior_draws$Tot_Eff,1.0 - 0.5 * alpha_ci),
        0.5 * quantile(object$posterior_draws$ACME_control + 
                         object$posterior_draws$ACME_treat,1.0 - 0.5 * alpha_ci),
        0.5 * quantile(object$posterior_draws$ADE_control + 
                         object$posterior_draws$ADE_treat,1.0 - 0.5 * alpha_ci),
        quantile( (object$posterior_draws$ACME_control + 
                     object$posterior_draws$ACME_treat) / 
                    (object$posterior_draws$ACME_control + 
                       object$posterior_draws$ACME_treat + 
                       object$posterior_draws$ADE_control + 
                       object$posterior_draws$ADE_treat), 1.0 - 0.5 * alpha_ci )
      )
  }
  
  summ
}
