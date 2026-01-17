#' @name print
#' 
#' @title Print bayesics objects.
#' 
#' @param x an object used to select a method.
#' @param ... optional arguments.
#' 
#' @returns None
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
#' print(fit1)
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
  cat("\n----------\n\nLinear regression fit using Bayesian techniques\n")
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
#' @method print np_glm_b
#' @export
print.np_glm_b = function(x, ...){
  cat("\n----------\n\nNon-parametric linear regression fit using Bayesian techniques\n")
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
#' @method print lm_b_bma
#' @export
print.lm_b_bma = function(x, ...){
  cat("\n----------\n\nBayesian model averaging for linear regression models\n")
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
#' @method print glm_b
#' @export
print.glm_b = function(x, ...){
  cat("\n----------\n\nGeneralized linear regression fit using Bayesian techniques\n")
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
