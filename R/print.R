#' @name print
#' 
#' @title Print \code{bayesics} Objects.
#' 
#' @param x an object used to select a method.
#' @param ... optional arguments.
#' 
#' @returns None
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
#' print(fit1)
#' }
#' 

#' @rdname print
#' @method print aov_b 
#' @export
print.aov_b = function(x, ...){
  cat("\n----------\n\nAnalysis of Variance fit using Bayesian techniques\n")
  cat("\n----------\n\n")
  print(x$formula)
  cat("\n----------\n\n") 
  if("BF_for_different_vs_same_means" %in% names(x)){
    cat(paste0(
      "Bayes factor in favor of the full vs. null model: ",
      format(signif(x$BF_for_different_vs_same_means, 3), 
             scientific = 
               (x$BF_for_different_vs_same_means > 1e3) | 
               (x$BF_for_different_vs_same_means < 1e-3))))
    
    cat("\n\n----------\n\n")
  }
  print(x$summary)
  cat("\n----------\n")
  cat(paste0("(Note: Lower and upper bounds are for the ",
             100 * x$CI_level,
             "% credible interval.)"))
}

#' @rdname print
#' @method print lm_b 
#' @export
print.lm_b = function(x, ...){
  
  if("lm_b_fits" %in% names(x)){
    header = 
      "\n----------\n\nBayesian model averaging for linear regression models\n"
  }else{
    header = 
      paste0("\n----------\n\n",
             ifelse(x$family$family == "gaussian",
                    "Linear ",
                    "Generalized linear ")
      ) |> 
      paste0("regression fit using Bayesian techniques",
             ifelse(x$model_type == "nonparametric",
                    " (non-parametric)",""),
             "\n")
  }
  
  cat(header)
  cat("\n----------\n\n")
  print(x$formula)
  cat("\n----------\n\n")
  print(x$summary)
  cat("\n----------\n")
  cat(paste0("(Note: Lower and upper bounds are for the ",
             100 * x$CI_level,
             "% credible interval.)"))
}


#' @rdname print
#' @method print mediate_b
#' @export
print.mediate_b = function(x, ...){
  cat("\n----------\n\nMediation analysis using Bayesian techniques\n")
  cat("\n----------\n\n")
  cat("Mediator model:\n")
  print(x$model_m$formula)
  cat("\nOutcome model:\n")
  print(x$model_y$formula)
  cat("\n----------\n\n")
  print(x$summary)
  cat("\n----------\n")
  cat(paste0("(Note: Lower and upper bounds are for the ",
               100 * x$CI_level,
               "% credible interval.)"))
}




#' @rdname print
#' @method print survfit_b 
#' @export
print.survfit_b = function(x, ...){
  cat("\n----------\n\nSemi-parametric survival curve fitting using Bayesian techniques\n")
  cat("\n----------\n\n")
  
  if(x$single_group_analysis){
  
    tibble::tibble(Interval = 
                     x$intervals |> 
                     apply(1,function(x) paste0("(",
                                                format(signif(x[1], 3)),
                                                ",",
                                                format(signif(x[2], 3)),
                                                ")")),
                   `Estimated rate` = 
                     x$posterior_parameters[,1] / 
                     x$posterior_parameters[,2],
                   `2.5%` = 
                     qgamma(0.025,
                            x$posterior_parameters[,1],
                            x$posterior_parameters[,2]),
                   `97.5%` =
                     qgamma(0.975,
                            x$posterior_parameters[,1],
                            x$posterior_parameters[,2]),
                   Shape = 
                     format(signif(x$posterior_parameters[,1], 3)),
                   Rate = 
                     format(signif(x$posterior_parameters[,2], 3))
    ) |> 
      print()
    
  }else{
    
    for(g in x$group_names){
      cat(g)
      cat("\n\n")
      
      tibble::tibble(Interval = 
                       x[[g]]$intervals |> 
                       apply(1,function(x) paste0("(",
                                                  format(signif(x[1], 3)),
                                                  ",",
                                                  format(signif(x[2], 3)),
                                                  ")")),
                     `Estimated rate` = 
                       x[[g]]$posterior_parameters[,1] / 
                       x[[g]]$posterior_parameters[,2],
                     `2.5%` = 
                       qgamma(0.025,
                              x[[g]]$posterior_parameters[,1],
                              x[[g]]$posterior_parameters[,2]),
                     `97.5%` =
                       qgamma(0.975,
                              x[[g]]$posterior_parameters[,1],
                              x[[g]]$posterior_parameters[,2]),
                     Shape = 
                       format(signif(x[[g]]$posterior_parameters[,1], 3)),
                     Rate = 
                       format(signif(x[[g]]$posterior_parameters[,2], 3))
      ) |> 
        print()
      
      cat("\n----------\n\n")
      
    }
    
  }
  
  cat("Note: The time-to-event data follows a piecewise exponential model.  Each interval follows an exponential distribution, whose rate has a posterior of Gamma(<Shape>,<Rate>).")
}


#' @rdname print
#' @method print b_procedure
#' @export
print.b_procedure = function(x, ...){
  cat(paste0("\n----------\n\n",
             x$name,
             " using Bayesian techniques\n\n----------\n\n"))
  
  # Data
  if(x$print_data){
    cat("Data: \n")
    print(x$data)
    cat("\n")
  }
  
  
  # Prior
  if(is.list(x$prior)){
    cat("\n\n")
    cat(x$prior$description)
    cat("\n")
    format(signif(x$prior$prior, 3), 
             scientific = FALSE) |> 
      noquote() |> 
      print()
  }else{
    
    cat(x$prior)
  }
  
  
  # Results
  ## Estimate, CI, ROPE, pdir
  if(isTRUE(x$display_as_matrices)){ # This is for chisq_test_b
    
    ## Get row and column numbers
    results = 
      x$results |> 
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
        is.null(x$sampling_design) ~ "",
        x$sampling_design == "multinomial" ~ "P_(row,col)",
        x$sampling_design == "fixed columns" ~ "P_(row|col)",
        x$sampling_design == "multinomial" ~ "P(col|row)"
      )
    
    ## Print posterior mean
    cat(paste0(
      "\n\nEstimated probabilities ",
      prob_type,
      ":\n"))
    x_matrix = 
      matrix(0.0,nR,nC,
             dimnames = dimnames(x$data))
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
               100 * x$CI_level, "% credible intervals: \n"))
    
    ci_lower = ci_upper = x_matrix
    ci_lower[cbind(results$row,
                   results$col)] = 
      results$Lower
    ci_upper[cbind(results$row,
                   results$col)] = 
      results$Upper
    credints = matrix("",nR,nC,
                     dimnames = dimnames(x$data))
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
    if(!is.null(x$ROPE)){
      cat(
        paste0("\n\n",
               x$ROPE$description,
               " is in the ROPE (i.e., between ",
               format(signif(x$ROPE$ROPE_lower_bound, 3), 
                      scientific = FALSE),
               " and ",
               format(signif(x$ROPE$ROPE_upper_bound, 3), 
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
    if(!is.null(x$pdir)){
      cat(paste0("\n\n",
                 x$pdir$description,
                 ": \n"))
      x_matrix[cbind(results$row,
                     results$col)] = 
        x$pdir$pdir
      x_matrix |> 
        signif(3) |> 
        format(scientific = FALSE) |> 
        noquote() |> 
        print()
    }
    
    
  }else{ #End: if(isTRUE(x$display_as_matrices))
    
    cat("\n\nPosterior Results:\n")
    
    for(j in 1:nrow(x$results)){
        cat(paste0("\n---",
                   x$results$Quantity[j],
                   "\n"))
        cat(
          paste0("      Estimate: ",
                 format(signif(x$results$`Post Mean`[j], 3), 
                        scientific = FALSE),
                 "\n      ",
                 x$CI_level*100,
                 "% CI: (",
                 format(signif(x$results$Lower[j], 3), 
                        scientific = FALSE),
                 ",",
                 format(signif(x$results$Upper[j], 3), 
                        scientific = FALSE),
                 ")")
        )
        if( ("ROPE_lower_bound" %in% colnames(x$results)) &&
             (!is.na(x$results$Pr_in_ROPE[j])) ){
          cat(
            paste0("\n      Probability that ",
                   tolower(x$results$Quantity[j]),
                   " is between ",
                   format(signif(x$results$ROPE_lower_bound[j], 3), 
                          scientific = FALSE),
                   " and ",
                   format(signif(x$results$ROPE_upper_bound[j], 3), 
                          scientific = FALSE),
                   ": ",
                   format(signif(x$results$Pr_in_ROPE[j], 3), 
                          scientific = FALSE))
          )
        }
    }
    
    ## PDir
    if(!is.null(x$pdir)){
      cat(paste0("\n\n",
                 x$pdir$description,
                 ": ",
                 format(signif(x$pdir$pdir, 3),
                        scientific = FALSE)))
    }
    
  }#End: if(!isTRUE(x$display_as_matrices))
  
  
  ## Overall ROPE (see chisq_test)
  if(!is.null(x$overall_ROPE)){
    cat(paste0("\n\n",
               x$overall_ROPE$description,
               ": ",
               format(signif(x$overall_ROPE$Pr_in_ROPE, 3), 
                      scientific = FALSE)))
  }
  
  
  ## Bayes factor
  if(!is.null(x$BF)){
    cat(paste0("\n\n",
               x$BF$description,
               ": ",
               format(signif(x$BF$BF, 3), 
                      scientific = FALSE),
               "\n    =>Level of evidence: ",
               x$BF$interpretation))
  }
  
  
  
  cat("\n\n----------\n\n")
  
  if(!is.null(x$notes)){
    for(j in 1:length(x$notes)){
      message(paste0(paste(rep("*",j),collapse=""),
                     "Note: ",
                     x$notes[j]))
    }
  }
  
  invisible(x)
}
