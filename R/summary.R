#' @name summary
#' 
#' @title Summarizing \code{bayesics} Fits
#' 
#' @description
#' \code{summary} method for objects of class \code{lm_b}, 
#' \code{aov_b}, \code{mediate_b}, \code{b_procedure}, and 
#' \code{survfit_b}.
#' 
#' 
#' @param object \code{bayesics} object
#' @param CI_level Posterior probability covered by credible interval.  
#' Unused for \code{b_procedure} objects.
#' @param interpretable_scale If a GLM is fit using 
#' \code{binomial(link="logit")}, \code{poisson(link="log")}, or 
#' \code{negbinom()}, and if \code{interpretable_scale = TRUE} 
#' then the results will be exponentiated.
#' @param print_results logical
#' @param ... optional arguments.  for \code{print.survfit_b}, 
#' this goes into `tibble::print.tbl_df`.
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
  if(print_results){
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
  }
  
  
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
  if(print_results){
    cat("\n----------\n")
    cat(paste0("(Note: Lower and upper bounds are for the ",
               100 * CI_level,
               "% credible interval.)\n"))
  }
  
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




#' @rdname summary
#' @method summary b_procedure
#' @export
summary.b_procedure = function(object,
                               ...){
  cat(paste0("\n----------\n\n",
             object$name,
             " using Bayesian techniques\n\n----------\n\n"))
  
  # Data
  if(object$print_data){
    cat("Data: \n")
    print(object$data)
    cat("\n")
  }
  
  
  # Prior
  if(is.list(object$prior)){
    
    cat("\n\n")
    cat(object$prior$description)
    cat("\n")
    format(signif(object$prior$prior, 3), 
           scientific = FALSE) |> 
      noquote() |> 
      print()
    
  }else{
    
    cat(object$prior)
    
  }
  
  
  # Results
  ## Estimate, CI, ROPE, pdir
  if(isTRUE(object$display_as_matrices)){ # This is for chisq_test_b
    
    ## Get row and column numbers
    results = 
      object$results |> 
      mutate(row = 
               as.integer(stringr::str_extract(.data$Quantity, "(?<=Row )\\d+")),
             col = 
               as.integer(stringr::str_extract(.data$Quantity, "(?<=Col )\\d+"))
      )
    nR = max(results$row)
    nC = max(results$col)
    
    ## Get the type of probability being modeled
    prob_type = 
      dplyr::case_when(
        is.null(object$sampling_design) ~ "",
        object$sampling_design == "multinomial" ~ "P_(row,col)",
        object$sampling_design == "fixed columns" ~ "P_(row|col)",
        object$sampling_design == "multinomial" ~ "P(col|row)"
      )
    
    ## Print posterior mean
    cat(paste0(
      "\n\nEstimated probabilities ",
      prob_type,
      ":\n"))
    x_matrix = 
      matrix(0.0,nR,nC,
             dimnames = dimnames(object$data))
    x_matrix[cbind(results$row,
                   results$col)] = 
      results$`Post Mean`
    x_matrix |> 
      signif(3) |> 
      format(scientific = FALSE) |> 
      noquote() |> 
      print()
    
    
    ## Print CIs
    cat(paste0("\n\n",
               100 * object$CI_level, "% credible intervals: \n"))
    
    ci_lower = ci_upper = x_matrix
    ci_lower[cbind(results$row,
                   results$col)] = 
      results$Lower
    ci_upper[cbind(results$row,
                   results$col)] = 
      results$Upper
    credints = matrix("",nR,nC,
                      dimnames = dimnames(object$data))
    for(i in 1:nR){
      for(j in 1:nC){
        credints[i, j] = paste0("(", format(signif(ci_lower[i,j], 3),
                                            scientific = FALSE),
                                ", ",
                                format(signif(ci_upper[i,j], 3),
                                       scientific = FALSE),
                                ")")
      }
    }
    credints |> 
      noquote() |> 
      print()
    
    
    ## Print ROPE
    if(!is.null(object$ROPE)){
      cat(
        paste0("\n\n",
               object$ROPE$description,
               " is in the ROPE (i.e., between ",
               format(signif(object$ROPE$ROPE_lower_bound, 3), 
                      scientific = FALSE),
               " and ",
               format(signif(object$ROPE$ROPE_upper_bound, 3), 
                      scientific = FALSE),
               "): \n")
      )
      
      x_matrix[cbind(results$row,
                     results$col)] = 
        results$Pr_in_ROPE
      x_matrix |> 
        signif(3) |> 
        format(scientific = FALSE) |> 
        noquote() |> 
        print()
      
    }
    
    
    ## Print pdir
    if(!is.null(object$pdir)){
      cat(paste0("\n\n",
                 object$pdir$description,
                 ": \n"))
      x_matrix[cbind(results$row,
                     results$col)] = 
        object$pdir$pdir
      x_matrix |> 
        signif(3) |> 
        format(scientific = FALSE) |> 
        noquote() |> 
        print()
    }
    
    
  }else{ #End: if(isTRUE(object$display_as_matrices))
    
    cat("\n\nPosterior Results:\n")
    
    for(j in 1:nrow(object$results)){
      cat(paste0("\n---",
                 object$results$Quantity[j],
                 "\n"))
      cat(
        paste0("      Estimate: ",
               format(signif(object$results$`Post Mean`[j], 3), 
                      scientific = FALSE),
               "\n      ",
               object$CI_level*100,
               "% CI: (",
               format(signif(object$results$Lower[j], 3), 
                      scientific = FALSE),
               ",",
               format(signif(object$results$Upper[j], 3), 
                      scientific = FALSE),
               ")")
      )
      if( ("ROPE_lower_bound" %in% colnames(object$results)) &&
          (!is.na(object$results$Pr_in_ROPE[j])) ){
        cat(
          paste0("\n      Probability that ",
                 tolower(object$results$Quantity[j]),
                 " is between ",
                 format(signif(object$results$ROPE_lower_bound[j], 3), 
                        scientific = FALSE),
                 " and ",
                 format(signif(object$results$ROPE_upper_bound[j], 3), 
                        scientific = FALSE),
                 ": ",
                 format(signif(object$results$Pr_in_ROPE[j], 3), 
                        scientific = FALSE))
        )
      }
    }
    
    ## PDir
    if(!is.null(object$pdir)){
      cat(paste0("\n\n",
                 object$pdir$description,
                 ": ",
                 format(signif(object$pdir$pdir, 3),
                        scientific = FALSE)))
    }
    
  }#End: if(!isTRUE(object$display_as_matrices))
  
  
  ## Overall ROPE (see chisq_test)
  if(!is.null(object$overall_ROPE)){
    cat(paste0("\n\n",
               object$overall_ROPE$description,
               ": ",
               format(signif(object$overall_ROPE$Pr_in_ROPE, 3), 
                      scientific = FALSE)))
  }
  
  
  ## Bayes factor
  if(!is.null(object$BF)){
    cat(paste0("\n\n",
               object$BF$description,
               ": ",
               format(signif(object$BF$BF, 3), 
                      scientific = FALSE),
               "\n    =>Level of evidence: ",
               object$BF$interpretation))
  }
  
  
  
  cat("\n\n----------\n\n")
  
  if(!is.null(object$notes)){
    for(j in 1:length(object$notes)){
      message(paste0(paste(rep("*",j),collapse=""),
                     "Note: ",
                     object$notes[j]))
    }
  }
  
  invisible(object)
}





#' @rdname summary
#' @method summary survfit_b
#' @export
summary.survfit_b = function(object, ...){
  cat("\n----------\n\nSemi-parametric survival curve fitting using Bayesian techniques\n")
  cat("\n----------\n\n")
  
  if(object$single_group_analysis){
    
    cat(paste0("Number of intervals: ",
               nrow(object$intervals),
               "\nSurvival curve fitted up to: ",
               max(object$intervals),
               "\n\n"))
    
    summary_object = 
      tibble::tibble(Interval = 
                       object$intervals |> 
                       apply(1,function(x) paste0("(",
                                                  format(signif(x[1], 3)),
                                                  ",",
                                                  format(signif(x[2], 3)),
                                                  ")")),
                     `Estimated rate` = 
                       object$posterior_parameters[,1] / 
                       object$posterior_parameters[,2],
                     `2.5%` = 
                       qgamma(0.025,
                              object$posterior_parameters[,1],
                              object$posterior_parameters[,2]),
                     `97.5%` =
                       qgamma(0.975,
                              object$posterior_parameters[,1],
                              object$posterior_parameters[,2]),
                     Shape = 
                       format(signif(object$posterior_parameters[,1], 3)),
                     Rate = 
                       format(signif(object$posterior_parameters[,2], 3))
      )
    
    print(summary_object, ...)
    
  }else{
    
    cat(paste0("\nNumber of intervals: ",
               nrow(object[[1]]$intervals),
               "\n\nSurvival curve fitted up to: ",
               max(object[[1]]$intervals),
               "\n\n"))
    
    
    for(g in object$group_names){
      cat(g)
      cat("\n\n")
      
      temp = 
        tibble::tibble(Interval = 
                         object[[g]]$intervals |> 
                         apply(1,function(x) paste0("(",
                                                    format(signif(x[1], 3)),
                                                    ",",
                                                    format(signif(x[2], 3)),
                                                    ")")),
                       `Estimated rate` = 
                         object[[g]]$posterior_parameters[,1] / 
                         object[[g]]$posterior_parameters[,2],
                       `2.5%` = 
                         qgamma(0.025,
                                object[[g]]$posterior_parameters[,1],
                                object[[g]]$posterior_parameters[,2]),
                       `97.5%` =
                         qgamma(0.975,
                                object[[g]]$posterior_parameters[,1],
                                object[[g]]$posterior_parameters[,2]),
                       Shape = 
                         format(signif(object[[g]]$posterior_parameters[,1], 3)),
                       Rate = 
                         format(signif(object[[g]]$posterior_parameters[,2], 3))
        )
      
      print(temp, ...)
      
      if(g == object$group_names[1]){
        summary_object = 
          temp |> 
          dplyr::mutate(Group = g)
      }else{
        summary_object = 
          bind_rows(
            summary_object,
            temp |> 
              dplyr::mutate(Group = g)
          )
      }
      summary_object = 
        summary_object |> 
        dplyr::relocate(Group)
      
      cat("\n----------\n\n")
      
    }
    
  }
  
  cat("Note: The time-to-event data follow a piecewise exponential model.  Each interval follows an exponential distribution, whose rate has a posterior of Gamma(<Shape>,<Rate>).\n")
  
  invisible(summary_object)
}

