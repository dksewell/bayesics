#' Bayesian Wilcoxon Rank Sum (aka Mann-Whitney U) and Signed Rank Analyses
#' 
#' 
#' @details
#' 
#' \strong{Bayesian Wilcoxon signed rank analysis}
#' For a single input vector or paired data, the Bayesian signed rank 
#' analysis will be performed.  The estimand is the proportion of (differenced) 
#' values that are positive.  For more information, see \link[DFBA]{dfba_wilcoxon} 
#' and vignette("dfba_wilcoxon",package = "DFBA").
#' 
#' \strong{Bayesian Wilcoxon rank sum/Mann-Whitney analysis}
#' For unpaired x and y inputs, the Bayesian rank sum analysis will be performed. 
#' The estimand is \eqn{\Omega_x:=\lim_{n\to\infty} \frac{U_x}{U_x + U_y}}, where 
#' \eqn{U_x} is the number of pairs \eqn{(i,j)} such that \eqn{x_i > y_j}, and 
#' vice versa for \eqn{U_y}.  That is, it is the population proportion of all 
#' untied pairs for which \eqn{x > y}.  Larger values imply that \eqn{x} is 
#' stochastically larger than \eqn{y}. For more information, see \link[DFBA]{dfba_mann_whitney} 
#' and vignette("dfba_mann_whitney",package = "DFBA").
#' 
#' 
#' @param x numeric vector of data values. Non-finite (e.g., infinite or 
#' missing) values will be omitted.
#' @param y an optional numeric vector of data values: as with x non-finite values will be omitted.
#' @param paired if \code{TRUE} and \code{y} is supplied, x-y will be the input of the Bayesian 
#' Wilcoxon signed rank test.
#' @param p numeric.
#' \itemize{
#'  \item Signed rank: \code{wilcox_test_b} will return the 
#' posterior probability that the population proportion of positive values 
#' (i.e., \eqn{x>y}) is greater than this value.
#'  \item Rank sum/Mann-Whitney U: \code{wilcox_test_b} will return the 
#' posterior probability that the \eqn{\Omega_x} (see details) is greater than 
#' this value.
#' }
#' @param ROPE If a single number, ROPE will be \code{p}\eqn{\pm}\code{ROPE}. 
#'  If a vector of length 2, these will serve as the ROPE bounds. Defaults to 
#'  \eqn{\pm 0.05}.
#' @param prior Prior used on the probability that x > y.  Either 
#' "uniform" (Beta(1,1)), or "centered" (Beta(2,2)). 
#' This is ignored if prior_shapes is provided.
#' @param prior_shapes Vector of length two, giving the shape parameters for the 
#' beta distribution that will act as the prior on the population proportions.
#' @param CI_level The posterior probability to be contained in the credible interval.
#' @param plot logical.  Should a plot be shown?
#' @param seed Always set your seed! (Unused for \eqn{\geq} 20 observations.)
#' 
#' @returns (returned invisible) If signed rank analysis is implemented, a list with the following:
#' \itemize{
#'  \item Posterior mean
#'  \item Credible interval
#'  \item Probability proportion of values that are positive is less than the reference p
#'  \item Probability proportion of values that are positive is in the ROPE
#'  \item Underlying DFBA object
#'  \item Prior and posterior plot
#'  \item Posterior parameters
#'  \item Bayes factor
#' }
#' If rank sum analysis is implemented, a list with the following:
#' \itemize{
#'  \item Posterior mean
#'  \item Credible interval
#'  \item Probability \eqn{\Omega_x} is less than the reference p
#'  \item Probability \eqn{\Omega_x} is in the ROPE
#'  \item Underlying DFBA object
#'  \item Prior and posterior plot
#'  \item Posterior parameters
#'  \item Bayes factor
##' }
#' 
#' @references 
#' Chechile, R.A. (2020). Bayesian Statistics for Experimental Scientists: A General Introduction to Distribution-Free Methods. Cambridge: MIT Press.
#' 
#' Chechile, R. A. (2018) A Bayesian analysis for the Wilcoxon signed-rank statistic. Communications in Statistics - Theory and Methods, https://doi.org/10.1080/03610926.2017.1388402
#' 
#' Chechile, R.A. (2020). A Bayesian analysis for the Mann-Whitney statistic. Communications in Statistics – Theory and Methods 49(3): 670-696. https://doi.org/10.1080/03610926.2018.1549247.
#' 
#' 
#' 
#' @examples
#' \donttest{
#' # Signed rank analysis
#' ## Generate data
#' N = 150
#' set.seed(2025)
#' test_data = 
#'   data.frame(x = rbeta(N,2,10),
#'              y = rbeta(N,5,10))
#' 
#' ## input differenced data
#' wilcoxon_test_b(test_data$x - test_data$y)
#' ## input paired data vectors individually
#' wilcoxon_test_b(test_data$x,
#'                 test_data$y,
#'                 paired = TRUE)
#' 
#' ## Use different priors
#' wilcoxon_test_b(test_data$x - test_data$y,
#'                 prior = "uniform")
#' wilcoxon_test_b(test_data$x - test_data$y,
#'                 prior_shapes = c(5,5))
#' 
#' ## Change ROPE bounds
#' wilcoxon_test_b(test_data$x - test_data$y,
#'                 ROPE = 0.1)
#' 
#' # Rank sum analysis
#' ## Generate data
#' set.seed(2025)
#' N = 150
#' x = rbeta(N,2,10)
#' y = rbeta(N + 1,5,10)
#' 
#' ## Perform analysis
#' wilcoxon_test_b(x,y)
#' }
#' 
#' 
#' @export

wilcoxon_test_b = function(x,
                           y,
                           paired = FALSE,
                           p = 0.5,
                           ROPE,
                           prior = "centered",
                           prior_shapes,
                           CI_level = 0.95,
                           plot = TRUE,
                           seed = 1){
  if( !("numeric" %in% class(x)) ) 
    stop("x must be numeric")
  if( !missing(y) && !("numeric" %in% class(y)) )
    stop("y must be numeric")
  
  
  ## Get ROPE
  if(missing(ROPE)){
    ROPE = 
      min(c(0.05,
            p / 2,
            (1.0 - p) / 2))
    ROPE_bounds = p + c(-1,1) * ROPE
  }else{
    if(!(length(ROPE) %in% 1:2))
      stop("ROPE must be of length 1 or 2")
    if(length(ROPE) == 1){
      if(ROPE < 0)
        stop("ROPE must be positive")
      if(ROPE > min(c(p,1.0 - p) / 2))
        stop("ROPE extends from p beyond 0 or 1")
      ROPE_bounds = p + c(-1,1) * ROPE
    }else{
      if(diff(ROPE) <= 0)
        stop("ROPE lower bound must be less than the upper bound")
      if(ROPE[1] < 0)
        stop("Lower bound of ROPE cannot be less than 0")
      if(ROPE[2] > 1)
        stop("Upper bound of ROPE cannot be greater than 1")
      ROPE_bounds = ROPE
    }
  }
  
  ## Prior distribution
  if(missing(prior_shapes)){
    prior = c("uniform",
              "centered")[pmatch(tolower(prior),
                                 c("uniform",
                                   "centered"))]
    
    if(prior == "uniform"){
      message("Prior shape parameters were not supplied.\nA uniform prior will be used.")
      prior_shapes = rep(1.0,2)
    }
    if(prior == "centered"){
      message("Prior shape parameters were not supplied.\nBeta(2,2) prior will be used.")
      prior_shapes = rep(2.0,2)
    }
  }else{
    if(any(prior_shapes <= 0))
      stop("Prior shape parameters must be positive.")
  }
  
  
  
  
  if(paired | missing(y)){
    
    # Do signed rank test
    
    ## Get outcome vector of differences
    if(!missing(y)){
      x = na.omit(x); y = na.omit(y)
      if(length(x) != length(y))
        stop("If paired = TRUE, then x and y must be of the same length.")
      x = x - y
    }
    
    ## Use DFBA package to do analysis
    dfba_object = 
      DFBA::dfba_wilcoxon(Y1 = x,
                          Y2 = numeric(length(x)),
                          a0 = prior_shapes[1],
                          b0 = prior_shapes[2],
                          prob_interval = CI_level)
    
    # Compute results
    ## If small, do:
    if("phipost" %in% names(dfba_object)){
      results = 
        list(posterior_mean = 
               sum(dfba_object$phiv * dfba_object$phipost),
             CI = 
               c(dfba_object$hdi_lower,
                 dfba_object$hdi_upper),
             Pr_less_than_p = 
               dfba_object$cumulative_phi[max(which(dfba_object$phiv <= p))],
             Pr_in_ROPE = 
               dfba_object$cumulative_phi[min(which(dfba_object$phiv >= ROPE_bounds[2]))] - 
               dfba_object$cumulative_phi[max(which(dfba_object$phiv <= ROPE_bounds[1]))],
             posterior_distribution = 
               tibble::tibble(`Pr(x>y)` = dfba_object$phiv,
                              `Posterior value` = dfba_object$phipost),
             dfba_wilcoxon_object = dfba_object
        )
      
      # Plot (if requested)
      if(plot){
        
        results$prob_plot = 
          tibble::tibble(x = seq(0.001,0.999,#seq(.Machine$double.eps,1.0 - .Machine$double.eps,
                                 l = 50)) |> 
          ggplot(aes(x=x)) +
          stat_function(fun = 
                          function(x){
                            dbeta(x,
                                  prior_shapes[1],
                                  prior_shapes[2])
                          },
                        aes(color = "Prior"),
                        linewidth = 2) + 
          geom_smooth(data = results$posterior_distribution |> 
                        dplyr::mutate(`Posterior value` = 
                                        .data$`Posterior value` / 
                                        max(.data$`Posterior value`) * 
                                        max(dbeta(seq(0.001,0.999,l = 50),
                                                  prior_shapes[1],
                                                  prior_shapes[2]))),
                      aes(x = .data$`Pr(x>y)`,
                          y = .data$`Posterior value`,
                          color = "Posterior"),
                      linewidth = 2,
                      method = "loess", 
                      span = 0.3,
                      se = FALSE) +
          scale_color_manual(values = c("Prior" = "#440154FF", 
                                        "Posterior" = "#FDE725FF")) +
          theme_classic(base_size = 15) +
          xlab("") + 
          ylab("") + 
          labs(color = "Distribution") + 
          ggtitle(ifelse(missing(y),"Pr(x > 0)","Pr(x > y)")) + 
          ylim(0,1.025 * max(dbeta(seq(0.001,0.999,l = 50),
                                   prior_shapes[1],
                                   prior_shapes[2])))
        
        suppressWarnings({
          print(results$prob_plot)
        })
        
      }
      
      
    }else{
      ## If large, do:
      results = 
        list(posterior_mean = 
               dfba_object$a_post / 
               sum( dfba_object$a_post + dfba_object$b_post),
             CI = 
               c(dfba_object$eti_lower,
                 dfba_object$eti_upper),
             Pr_less_than_p = 
               pbeta(p,
                     dfba_object$a_post,
                     dfba_object$b_post),
             Pr_in_ROPE = 
               pbeta(ROPE_bounds[2],
                     dfba_object$a_post,
                     dfba_object$b_post) - 
               pbeta(ROPE_bounds[1],
                     dfba_object$a_post,
                     dfba_object$b_post),
             dfba_wilcoxon_object = dfba_object
        )
      
      
      # Plot (if requested)
      if(plot){
        
        results$prob_plot = 
          tibble::tibble(x = seq(0.001,0.999,#seq(.Machine$double.eps,1.0 - .Machine$double.eps,
                                 l = 50)) |> 
          ggplot(aes(x=x)) +
          stat_function(fun = 
                          function(x){
                            dbeta(x,
                                  prior_shapes[1],
                                  prior_shapes[2])
                          },
                        aes(color = "Prior"),
                        linewidth = 2) + 
          
          stat_function(fun = 
                          function(x){
                            dbeta(x,
                                  dfba_object$a_post,
                                  dfba_object$b_post)
                          },
                        aes(color = "Posterior"),
                        linewidth = 2) + 
          
          scale_color_manual(values = c("Prior" = "#440154FF", 
                                        "Posterior" = "#FDE725FF")) +
          theme_classic(base_size = 15) +
          xlab("") + 
          ylab("") + 
          labs(color = "Distribution") + 
          ggtitle(ifelse(missing(y),"Pr(x > 0)","Pr(x > y)"))
        
        
        print(results$prob_plot)
        
      }
      
      # Add posterior parameters to returned object
      results$posterior_parameters = 
        c(shape_1 = dfba_object$a_post,
          shape_2 = dfba_object$b_post)
      
    }
    
    results$BF_for_phi_gr_onehalf_vs_phi_less_onehalf = 
      dfba_object$BF10
    bf_max = max(results$BF_for_phi_gr_onehalf_vs_phi_less_onehalf,
                 1.0 / results$BF_for_phi_gr_onehalf_vs_phi_less_onehalf)
    BF_evidence =
      ifelse(bf_max <= 3.2,
             "Not worth more than a bare mention",
             ifelse(bf_max <= 10,
                    "Substantial",
                    ifelse(bf_max <= 100,
                           "Strong",
                           "Decisive")))
    
    
    # Print results
    cat("\n----------\n\nWilcoxon signed-rank analysis using Bayesian techniques\n")
    cat("\n----------\n\n")
    cat(paste0("Prior used: Beta(", 
               format(signif(prior_shapes[1], 3), 
                      scientific = FALSE),
               ",",
               format(signif(prior_shapes[2], 3), 
                      scientific = FALSE),
               ")\n\n"))
    cat(paste0("Posterior mean: ", 
               format(signif(results$posterior_mean, 3), 
                      scientific = FALSE),
               "\n\n"))
    cat(paste0(100 * CI_level,
               "% credible interval: (", 
               format(signif(results$CI[1], 3), 
                      scientific = FALSE),
               ", ",
               format(signif(results$CI[2], 3), 
                      scientific = FALSE),
               ")\n\n"))
    cat(paste0("Probability that Pr(x > y) > ",
               format(signif(p, 3), 
                      scientific = FALSE),
               " = ",
               format(signif(1.0 - results$Pr_less_than_p, 3), 
                      scientific = FALSE),
               "\n\n"))
    cat(paste0("Probability that Pr(x > y) is in the ROPE, defined to be (",
               format(signif(ROPE_bounds[1], 3), 
                      scientific = FALSE),
               ",",
               format(signif(ROPE_bounds[2], 3), 
                      scientific = FALSE),
               ") = ",
               format(signif(results$Pr_in_ROPE, 3), 
                      scientific = FALSE),
               "\n\n")) 
    if(p == 0.5){
      cat(paste0("Bayes factor in favor of phi>0.5 vs. phi<=0.5: ",
                 format(signif(results$BF_for_phi_gr_onehalf_vs_phi_less_onehalf, 3), 
                        scientific = FALSE),
                 ";\n      =>Level of evidence: ", 
                 BF_evidence,
                 "\n\n")) 
      
    }
    cat("\n----------\n\n")
    
    
    
    
    invisible(results)
    
  }else{#End: Wilcoxon signed rank analysis
    
    # Do rank sum test
    
    ## Use DFBA package to do analysis
    dfba_object = 
      DFBA::dfba_mann_whitney(E = x,
                              C = y,
                              a0 = prior_shapes[1],
                              b0 = prior_shapes[2],
                              prob_interval = CI_level)
    
    # Compute results
    ## If small, do:
    if("omegapost" %in% names(dfba_object)){
      results = 
        list(posterior_mean = 
               sum(dfba_object$omega_E * dfba_object$omegapost),
             CI = 
               c(dfba_object$eti_lower,
                 dfba_object$eti_upper),
             Pr_less_than_p = 
               dfba_object$cumulative_omega[max(which(dfba_object$omega_E <= p))],
             Pr_in_ROPE = 
               dfba_object$cumulative_omega[min(which(dfba_object$omega_E >= ROPE_bounds[2]))] - 
               dfba_object$cumulative_omega[max(which(dfba_object$omega_E <= ROPE_bounds[1]))],
             posterior_distribution = 
               tibble::tibble(`Pr(x>y)` = dfba_object$omega_E,
                              `Posterior value` = dfba_object$omegapost),
             dfba_wilcoxon_object = dfba_object
        )
      
      # Plot (if requested)
      if(plot){
        
        results$prob_plot = 
          tibble::tibble(x = seq(0.001,0.999,#seq(.Machine$double.eps,1.0 - .Machine$double.eps,
                                 l = 50)) |> 
          ggplot(aes(x=x)) +
          stat_function(fun = 
                          function(x){
                            dbeta(x,
                                  prior_shapes[1],
                                  prior_shapes[2])
                          },
                        aes(color = "Prior"),
                        linewidth = 2) + 
          geom_smooth(data = results$posterior_distribution |> 
                        mutate(`Posterior value` = 
                                 .data$`Posterior value` / 
                                 max(.data$`Posterior value`) * 
                                 max(dbeta(seq(0.001,0.999,l = 50),
                                           prior_shapes[1],
                                           prior_shapes[2]))),
                      aes(x = .data$`Pr(x>y)`,
                          y = .data$`Posterior value`,
                          color = "Posterior"),
                      linewidth = 2,
                      method = "loess", 
                      span = 0.3,
                      se = FALSE) +
          scale_color_manual(values = c("Prior" = "#440154FF", 
                                        "Posterior" = "#FDE725FF")) +
          theme_classic(base_size = 15) +
          xlab("") + 
          ylab("") + 
          labs(color = "Distribution") + 
          ggtitle(expression(Omega[x])) + 
          ylim(0,1.025 * max(dbeta(seq(0.001,0.999,l = 50),
                                   prior_shapes[1],
                                   prior_shapes[2])))
        
        suppressWarnings({
          print(results$prob_plot)
        })
        
      }
      
      
    }else{
      ## If large, do:
      results = 
        list(posterior_mean = 
               dfba_object$a_post / 
               sum( dfba_object$a_post + dfba_object$b_post),
             CI = 
               c(dfba_object$eti_lower,
                 dfba_object$eti_upper),
             Pr_less_than_p = 
               pbeta(p,
                     dfba_object$a_post,
                     dfba_object$b_post),
             Pr_in_ROPE = 
               pbeta(ROPE_bounds[2],
                     dfba_object$a_post,
                     dfba_object$b_post) - 
               pbeta(ROPE_bounds[1],
                     dfba_object$a_post,
                     dfba_object$b_post),
             dfba_wilcoxon_object = dfba_object
        )
      
      
      # Plot (if requested)
      if(plot){
        
        results$prob_plot = 
          tibble::tibble(x = seq(0.001,0.999,#seq(.Machine$double.eps,1.0 - .Machine$double.eps,
                                 l = 50)) |> 
          ggplot(aes(x=x)) +
          stat_function(fun = 
                          function(x){
                            dbeta(x,
                                  prior_shapes[1],
                                  prior_shapes[2])
                          },
                        aes(color = "Prior"),
                        linewidth = 2) + 
          
          stat_function(fun = 
                          function(x){
                            dbeta(x,
                                  dfba_object$a_post,
                                  dfba_object$b_post)
                          },
                        aes(color = "Posterior"),
                        linewidth = 2) + 
          
          scale_color_manual(values = c("Prior" = "#440154FF", 
                                        "Posterior" = "#FDE725FF")) +
          theme_classic(base_size = 15) +
          xlab("") + 
          ylab("") + 
          labs(color = "Distribution") + 
          ggtitle(expression(Omega[x])) 
        
        print(results$prob_plot)
        
      }
      
      # Add posterior parameters to returned object
      results$posterior_parameters = 
        c(shape_1 = dfba_object$a_post,
          shape_2 = dfba_object$b_post)
      
    }
    
    results$BF_for_Omegax_gr_onehalf_vs_Omegax_less_onehalf = 
      dfba_object$BF10
    bf_max = max(results$BF_for_Omegax_gr_onehalf_vs_Omegax_less_onehalf,
                 1.0 / results$BF_for_Omegax_gr_onehalf_vs_Omegax_less_onehalf)
    BF_evidence =
      ifelse(bf_max <= 3.2,
             "Not worth more than a bare mention",
             ifelse(bf_max <= 10,
                    "Substantial",
                    ifelse(bf_max <= 100,
                           "Strong",
                           "Decisive")))
    
    # Print results
    cat("\n----------\n\nWilcoxon rank sum analysis using Bayesian techniques\n")
    cat("\n----------\n\n")
    cat(
      "NOTE: Estimand is Omega_x := Proportion of (non-tied) pairs where x is bigger than y"
    )
    cat("\n\n----------\n\n")
    cat(paste0("Prior used: Beta(", 
               format(signif(prior_shapes[1], 3), 
                      scientific = FALSE),
               ",",
               format(signif(prior_shapes[2], 3), 
                      scientific = FALSE),
               ")\n\n"))
    cat(paste0("Posterior mean: ", 
               format(signif(results$posterior_mean, 3), 
                      scientific = FALSE),
               "\n\n"))
    cat(paste0(100 * CI_level,
               "% credible interval: (", 
               format(signif(results$CI[1], 3), 
                      scientific = FALSE),
               ", ",
               format(signif(results$CI[2], 3), 
                      scientific = FALSE),
               ")\n\n"))
    cat(paste0("Probability that Omega_x > ",
               format(signif(p, 3), 
                      scientific = FALSE),
               " = ",
               format(signif(1.0 - results$Pr_less_than_p, 3), 
                      scientific = FALSE),
               "\n\n"))
    cat(paste0("Probability that Omega_x is in the ROPE, defined to be (",
               format(signif(ROPE_bounds[1], 3), 
                      scientific = FALSE),
               ",",
               format(signif(ROPE_bounds[2], 3), 
                      scientific = FALSE),
               ") = ",
               format(signif(results$Pr_in_ROPE, 3), 
                      scientific = FALSE),
               "\n\n"))
    if(p == 0.5){
      cat(paste0("Bayes factor in favor of phi>0.5 vs. phi<=0.5: ",
                 format(signif(results$BF_for_Omegax_gr_onehalf_vs_Omegax_less_onehalf, 3), 
                        scientific = FALSE),
                 ";\n      =>Level of evidence: ", 
                 BF_evidence,
                 "\n\n")) 
      
    }
    cat("\n----------\n\n")
    
    
    
    
    invisible(results)
    
    
  }
  
}
