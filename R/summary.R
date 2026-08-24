#' @name summary
#' 
#' @title Summary Functions for bayesics Objects
#' 
#' @param object bayesics object
#' @param CI_level Posterior probability covered by credible interval
#' @param interpretable_scale If a GLM is fit using 
#' \code{binomial(link="logit")}, \code{poisson(link="log")}, or 
#' \code{negbinom()}, and if \code{interpretable_scale = TRUE} 
#' then the results will be exponentiated.
#' @param print_results logical
#' @param ... optional arguments.
#' 
#' @returns tibble with summary values
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
#' summary(fit1)
#' }
#' 
#' @export

#' @rdname summary
#' @method summary lm_b 
#' @export
summary.lm_b = function(object,
                        CI_level = 0.95,
                        interpretable_scale = TRUE,
                        print_results = TRUE,
                        ...){
  alpha = 1 - CI_level
  p = 
    length(setdiff(object$summary$Variable,
                   "log(phi)"))
  summ = object$summary
  
  if("posterior_covariance" %in% names(object)){ # Handles lm, glm\IS, np_glm\bootstrapping
    
    summ$Lower = 
      qlst(alpha / 2.0,
           object$df,
           object$summary$`Post Mean`,
           sqrt(diag(as.matrix(object$posterior_covariance))))
    summ$Upper = 
      qlst(1.0 - alpha / 2.0,
           object$df,
           object$summary$`Post Mean`,
           sqrt(diag(as.matrix(object$posterior_covariance))))
    
  }else{
    if("importance_sampling_weights" %in% names(object)){ # Handles glm IS
      
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
      
      
    }else{
      if("posterior_draws" %in% names(object)){ # Handles np_glm bootstrapping, bma_inference
        
        summ$Lower = 
          object$posterior_draws[,1:nrow(summ)] |> 
          as.matrix() |> 
          apply(2,quantile,prob = alpha / 2)
        summ$Upper = 
          object$posterior_draws[,1:nrow(summ)] |> 
          as.matrix() |> 
          apply(2,quantile,prob = 1.0 - alpha / 2)
        
      }#End: posterior_draws if
    }#End: importance_sampling_weights ifelse
  }#End: posterior_covariance ifelse
  
  
  # Make interpretable if family != "gaussian"
  if(( (object$family$family == "binomial") & 
       (object$family$link != "logit") ) | 
     ( (object$family$family == "poisson") & 
       (object$family$link != "log") ) |
     ( object$family$family == "gaussian" )){
    interpretable_scale = FALSE
  }
  
  
  # Print analysis high level info
  if("lm_b_fits" %in% names(object)){
    header = 
      "\n----------\n\nBayesian model averaging for linear regression models\n"
  }else{
    header = 
      paste0("\n----------\n\n",
             ifelse(object$family$family == "gaussian",
                    "Linear ",
                    "Generalized linear ")
      ) |> 
      paste0("regression fit using Bayesian techniques",
             ifelse(object$model_type == "nonparametric",
                    " (non-parametric)",""),
             "\n")
  }
  cat(header)
  cat("\n----------\n\n")
  print(object$formula)
  cat("\n----------\n\n")
  
  
  # Print if scale is transformed
  if(interpretable_scale){
    if(print_results){
      paste0("Values given in terms of ",
             ifelse(object$family$family == "binomial",
                    "odds ratios",
                    "rate ratios")
      ) |> 
        cat()
      cat("\n\n----------\n\n")
    }
    
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
  
  # Add sigma^2 if lm_b or lm_b_bma
  if("sigma_sq" %in% names(object)){
    if("posterior_parameters" %in% names(object)){ #Handles lm_b
      summ = 
        bind_rows(
          summ,
          tibble(Variable = "Residual variance",
                 `Post Mean` = object$sigma_sq["Estimate"],
                 Lower = 
                   extraDistr::qinvgamma(alpha / 2.0,
                                         0.5 * object$posterior_parameters$a_tilde,
                                         0.5 * object$posterior_parameters$b_tilde),
                 Upper = 
                   extraDistr::qinvgamma(1.0 - alpha / 2.0,
                                         0.5 * object$posterior_parameters$a_tilde,
                                         0.5 * object$posterior_parameters$b_tilde),
                 `Prob Dir` = NA,
                 ROPE = NA,
                 `ROPE bounds` = "(NA,NA)"
          )
        )
    }
    if("posterior_draws" %in% names(object)){ #Handles lm_b_bma
      summ = 
        bind_rows(
          summ,
          tibble(Variable = "Residual variance",
                 `Post Mean` = object$sigma_sq["Estimate"],
                 Lower = 
                   quantile(unlist(object$posterior_draws[,p+1]),
                            alpha / 2.0),
                 Upper = 
                   quantile(unlist(object$posterior_draws[,p+1]),
                            1.0 - alpha / 2.0),
                 `Prob Dir` = NA,
                 ROPE = NA,
                 `ROPE bounds` = "(NA,NA)"
          )
        )
    }
    
  }#End: add in s^2
  
  # Print results
  if(print_results) print(summ)
  
  # Print CI level
  cat("\n----------\n")
  cat(paste0("(Note: Lower and upper bounds are for the ",
             100 * CI_level,
             "% credible interval.)\n"))
  
  invisible(summ)
}

#' @rdname summary
#' @method summary aov_b 
#' @export
summary.aov_b = function(object,
                         CI_level = 0.95,
                         print_results = TRUE,
                         ...){
  alpha = 1 - CI_level
  
  summary_object = 
    list(
      summary = object$summary,
      pw_summary = 
        object$pairwise_summary |> 
        as.data.frame()
    )
  
  if("BF_for_different_vs_same_means" %in% names(object)){
    
    bf_max = 
      max(object$BF_for_different_vs_same_means, 
          1.0 / object$BF_for_different_vs_same_means)
    
    summary_object$BF = 
      list(BF = 
             object$BF_for_different_vs_same_means,
           interpretation = 
             paste0(
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
    
    if(print_results){
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
    
  }
  
  if(print_results) cat("\n\n\n\n--- Summary of factor level means ---\n")
  summary_object$summary$Lower = 
    c(extraDistr::qlst(alpha/2, 
                       df = object$posterior_parameters$a_g,
                       mu = object$posterior_parameters$mu_g,
                       sigma = sqrt(object$posterior_parameters$b_g / object$posterior_parameters$nu_g / object$posterior_parameters$a_g)),
      extraDistr::qinvgamma(alpha/2, 
                            alpha = object$posterior_parameters$a_g/2, 
                            beta = object$posterior_parameters$b_g/2))
  summary_object$summary$Upper = 
    c(extraDistr::qlst(1 - alpha/2, 
                       df = object$posterior_parameters$a_g,
                       mu = object$posterior_parameters$mu_g,
                       sigma = sqrt(object$posterior_parameters$b_g / object$posterior_parameters$nu_g / object$posterior_parameters$a_g)),
      extraDistr::qinvgamma(1 - alpha/2, 
                            alpha = object$posterior_parameters$a_g/2, 
                            beta = object$posterior_parameters$b_g/2))
  if(print_results) print(summary_object$summary)
  
  
  if(print_results) cat("\n\n\n\n--- Summary of pairwise differences ---\n")
  temp = 
    combn(1:length(levels(object$data$group)),2)
  for(i in 1:nrow(summary_object$pw_summary)){
    summary_object$pw_summary[i,c("Lower","Upper")] = 
      quantile(object$posterior_draws[,temp[1,i]] - 
                 object$posterior_draws[,temp[2,i]],
               probs = c(alpha/2, 
                         1 - alpha/2))
  }
  summary_object$pw_summary = as_tibble(summary_object$pw_summary)
  if(print_results) print(summary_object$pw_summary)
  if(print_results) cat("\n\n   *Note: EPR (Exceedence in Pairs Rate) for a Comparison of g-h = Pr(Y_(gi) > Y_(hi)|parameters) ")
  
  if(is.null(object$contrasts)){
    invisible(summary_object)
  }else{
    
    summary_object$contrasts = 
      list(L = object$contrasts$L)
    
    if(print_results) cat("\n\n\n\n--- Summary of Contrasts ---\n")
    summary_object$contrasts$summary = object$contrasts$summary
    contrast_draws = 
      tcrossprod(object$posterior_draws[,grep("mean_",colnames(object$posterior_draws))],
                 object$contrasts$L)
    summary_object$contrasts$summary$Lower = 
      apply(contrast_draws,2,quantile,probs = alpha/2)
    summary_object$contrasts$summary$Upper = 
      apply(contrast_draws,2,quantile,probs = 1 - alpha/2)
    
    if(print_results) print(summary_object$contrasts$summary)
    
    invisible(summary_object)
  }
}



#' @rdname summary
#' @method summary mediate_b 
#' @export
summary.mediate_b = function(object,
                             CI_level = 0.95,
                             print_results = TRUE,
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
  
  if(print_results) print(summ)
  invisible(summ)
}
