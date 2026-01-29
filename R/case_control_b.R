#' Case-Control Analysis
#' 
#' Bayesian analysis of a case-control study (without covariates).
#' 
#' @details
#' If \code{large_sample_approx = TRUE} (the default if left missing and all 
#' cell counts are at least 5), then the likelihood is
#' \deqn{
#'  \log(\hat\omega) \sim N\left(\log(\omega),\frac{1}{n_{11}} + \frac{1}{n_{12}} + 
#'  \frac{1}{n_{21}} + \frac{1}{n_{22}} \right),
#' }
#' where \eqn{\omega} is the odds ratio, \eqn{\hat\omega} is the 
#' empirical odds ratio, \eqn{n_{ij}}, \eqn{i,j = 1,2} are the cells of the 
#' 2x2 contingency table. The prior on \eqn{\log\omega} is
#' \deqn{
#'  \log\omega \sim N(a,b^2).
#' }
#' 
#' If the large sample approximation is not used, then inference is made on 
#' the odds ratio by instead putting uniform priors on \eqn{\Pr(exposure|outcome)}.
#' 
#' 
#' 
#' 
#' @param cases vector of length 2, giving the numbers at risk and not at risk,
#' respectively, for cases
#' @param controls vector of length 2, giving the numbers at risk and not at risk,
#' respectively, for controls
#' @param x 2x2 contingency table.  The rows should depict the at risk status 
#' (first row is at risk, second row is not at risk), and the columns should 
#' depict the case control status (first column is case, second column is control).
#' @param large_sample_approx If all cell counts of \code{x} are not too low 
#' (\eqn{\geq 5}) then use the approximation that the empirical log odds are 
#' normally distributed.  (See details for more.)  If missing, this will be 
#' set to \code{TRUE} iff all cell counts are greater than or equal to 5.
#' @param ROPE ROPE for odds ratio. Provide either a single value or a vector 
#' of length two.  If the former, the ROPE will be taken as (1/ROPE,ROPE).  
#' If the latter, these will be the bounds of the ROPE.
#' @param prior_mean numeric.  The prior mean on the log odds ratio
#' @param prior_sd numeric.  The prior sd on the log odds ratio. See details 
#' for default values.
#' @param plot logical.  Should a plot be shown?
#' @param CI_level The posterior probability to be contained in the 
#' credible interval.
#' @param seed integer.  Always set your seed!!! (ignored if \code{large_sample_approx = TRUE}.)
#' @param mc_error The relative monte carlo error of the quantiles of the CIs. 
#' (ignored if \code{large_sample_approx = TRUE}.)
#' 
#' @returns (returned invisible) list including the following:
#' \itemize{
#'  \item \code{data}: data
#'  \item \code{posterior_mean}: posterior mean of the odds ratio (cases vs. controls)
#'  \item \code{CI}: Credible interval
#'  \item \code{Pr_oddsratio_in_ROPE}: Probability the odds ratio (cases vs. controls) is in the ROPE
#'  \item \code{posterior_draws}: posterior draws of the odds ratio (cases vs. controls)
#'  \item \code{or_plot}: odds ratio (cases vs. controls) posterior plot
#' }
#' 
#' @examples
#' case_control_b(matrix(c(8,47,1,26),2,2))
#' 
#' case_control_b(c(8,47),
#'                c(1,26))
#' 
#' 
#' 
#' @export



case_control_b = function(cases,
                          controls,
                          x,
                          large_sample_approx,
                          ROPE,
                          prior_mean = 0.0,
                          prior_sd = log(10.0) / 1.96,
                          plot = TRUE,
                          CI_level = 0.95,
                          seed = 1,
                          mc_error = 0.005){
  
  alpha_ci = 1.0 - CI_level
  
  # Get 2x2 table, and do checks along the way
  if(missing(x)){
    if(missing(cases))
      stop("Either x or cases must be provided")
    if("table" %in% class(cases)){
      x = 
        cases |> 
        matrix(2,2)
    }
    if("matrix" %in% class(cases)){
      x = cases
    }
    if( ("integer" %in% class(cases)) | ("numeric" %in% class(cases)) ){
      if(length(cases) != 2) 
        stop("Length of cases must be 2")
      if(missing(controls))
        stop("Must supply both cases and controls")
      if(length(controls) != 2) 
        stop("Length of controls must be 2")
      
      x = cbind(cases,controls)
    }
  }else{
    if( !("table" %in% class(x)) & !("matrix" %in% class(x)))
      stop("x must be either a table or a 2x2 matrix")
    if(!all(near(dim(x),c(2,2))))
      stop("x must be 2x2")
  }
  
  # Get estimate of odds ratio
  or_hat =
    x[1,1] * x[2,2] / 
    x[1,2] / x[2,1]
  log_or_hat = log(or_hat)
  
  # Get ROPE
  if(missing(ROPE)){
    ROPE = c(1.0 / 1.125, 1.125)
    # From Kruchke (2018) on rate ratios from FDA <1.25. (Use half of small effect size for ROPE, hence 0.25/2) 
    #   Use the same thing for odds ratios.
  }else{
    if(length(ROPE) > 2) stop("ROPE must be given as an upper bound, or given as both lower and upper bounds.")
    if((length(ROPE) > 1) & (ROPE[1] >= ROPE[2])) stop("ROPE lower bound must be smaller than ROPE upper bound")
    if(length(ROPE) == 1) ROPE = c(1.0 / ROPE, ROPE)
  }
  
  
  # Determine if large sample approximation should be used
  if(missing(large_sample_approx)){
    large_sample_approx = 
      all(c(x) >= 5)
  }
  
  
  # Perform inference
  colnames(x) = c("Cases","Controls")
  rownames(x) = c("At risk","Not at risk")
  results = list(data=x)
  
  if(large_sample_approx){
    
    ## Get posterior parameters
    s2 = sum(1.0 / c(x))
    results$posterior_parameters = 
      c((prior_sd^2 * log_or_hat + s2 * prior_mean) / 
          (prior_sd^2 + s2),
        sqrt(prior_sd^2 * s2 / (prior_sd^2 + s2)))
    names(results$posterior_parameters) = c("mean","sd")
    
    ## Get point estimates
    results$posterior_mean = 
      exp(results$posterior_parameters["mean"] + 
            0.5 * results$posterior_parameters["sd"]^2)
    results$posterior_median = 
      exp(results$posterior_parameters["mean"])
    
    ## Get interval estimates
    results$CI = 
      exp(
        qnorm(c(0.5 * alpha_ci,1.0 - 0.5 * alpha_ci),
              results$posterior_parameters["mean"],
              results$posterior_parameters["sd"])
      )
    
    ## Get ROPE
    results$Pr_oddsratio_in_ROPE = 
      pnorm(log(ROPE[2]),
            results$posterior_parameters["mean"],
            results$posterior_parameters["sd"]) -
      pnorm(log(ROPE[1]),
            results$posterior_parameters["mean"],
            results$posterior_parameters["sd"])
    
    # Plot (if requested)
    if(plot){
      results$or_plot = 
        tibble::tibble(x = seq(qlnorm(0.005,
                                      results$posterior_parameters["mean"],
                                      results$posterior_parameters["sd"]),
                               qlnorm(0.995,
                                      results$posterior_parameters["mean"],
                                      results$posterior_parameters["sd"]),
                               l = 50)) |> 
        ggplot(aes(x=x)) +
        stat_function(fun = 
                        function(x){
                          dlnorm(x,
                                 prior_mean,
                                 prior_sd)
                        },
                      aes(color = "Prior"),
                      linewidth = 2) + 
        stat_function(fun = 
                        function(x){
                          dlnorm(x,
                                 results$posterior_parameters["mean"],
                                 results$posterior_parameters["sd"])
                        },
                      aes(color = "Posterior"),
                      linewidth = 2) + 
        scale_color_manual(values = c("Prior" = "#440154FF", 
                                      "Posterior" = "#FDE725FF")) +
        theme_classic(base_size = 15) +
        xlab("") + 
        ylab("") + 
        labs(color = "Distribution") + 
        ggtitle("Population proportion")
      
      print(results$or_plot)
      
    }
    
    
    
  }else{#End: large sample approx
    set.seed(seed)
    message("Cell sizes were too small for large sample approximation.\nInstead, setting uniform prior on Pr(exposure|outcome) and making exact finite sample inference.")
    
    # Get posterior parameters
    post_shapes = 
      t(x) + 1.0
    
    # Get posterior draws
    ## Get preliminary draws
    p1_draws = 
      rbeta(500,
            post_shapes[1,1],
            post_shapes[1,2])
    p2_draws = 
      rbeta(500,
            post_shapes[2,1],
            post_shapes[2,2])
    odds_ratios = 
      p1_draws / (1.0 - p1_draws) * (1.0 - p2_draws) / p2_draws
    fhat = 
      density(odds_ratios,
              from = 0.0 + .Machine$double.eps)
    n_draws = 
      0.5 * alpha_ci * (1.0 - 0.5 * alpha_ci) *
      (
        qnorm(0.5 * (1.0 - 0.99)) / 
          mc_error /
          fhat$y[which.min(abs(fhat$x - 
                                 quantile(odds_ratios, 0.5 * alpha_ci)))]
      )^2 |> 
      round()
    
    ## Finish posterior draws
    p1_draws = 
      c(p1_draws,
        rbeta(n_draws - length(p1_draws),
              post_shapes[1,1],
              post_shapes[1,2]))
    p2_draws = 
      c(p2_draws,
        rbeta(n_draws - length(p2_draws),
              post_shapes[2,1],
              post_shapes[2,2]))
    
    odds_ratios = 
      p1_draws / (1.0 - p1_draws) * (1.0 - p2_draws) / p2_draws
    
    ## Get point estimates
    results$posterior_mean = 
      mean(odds_ratios)
    
    ## Get interval estimates
    results$CI = 
      quantile(odds_ratios,c(0.5 * alpha_ci,
                             1.0 - 0.5 * alpha_ci))
    
    ## Get ROPE
    results$Pr_oddsratio_in_ROPE = 
      mean( (odds_ratios <= ROPE[2]) & 
              (odds_ratios >= ROPE[1]) )
    
    ## Save posterior draws
    results$posterior_draws = 
      odds_ratios
    
    
    ## Plot if requested
    if(plot){
      results$or_plot = 
        data.frame(or = 
                     odds_ratios[which(odds_ratios <= quantile(odds_ratios,0.99))]) |> 
        ggplot(aes(x = .data$or)) + 
        geom_density(bounds = c(0,Inf),
                     linewidth = 2) + 
        theme_classic(base_size = 15) +
        xlab("") + 
        ylab("") + 
        ggtitle("Posterior of odds ratio (cases vs. controls)")
      
      print(results$or_plot)
    }
    
  }#End: small sample inference
  
  
  # Print results
  message("\n----------\n\nCase-control analysis using Bayesian techniques\n")
  message("\n----------\n\n")
  message("Data: \n")
  message(paste(capture.output(print(x)),
                collapse = "\n"))
  message("\n\n")
  if(large_sample_approx){
    message(paste0("Prior used on log odds ratio: N(", 
               format(signif(prior_mean, 3), 
                      scientific = FALSE),
               ",",
               format(signif(prior_sd, 3), 
                      scientific = FALSE),
               ")\n\n"))
  }
  message(paste0("Posterior mean of the odds ratio: ", 
             format(signif(results$posterior_mean, 3), 
                    scientific = FALSE),
             "; Population 2 = ",
             "\n\n"))
  if(large_sample_approx){
    message(paste0("Posterior median of the odds ratio: ", 
               format(signif(results$posterior_median, 3), 
                      scientific = FALSE),
               "; Population 2 = ",
               "\n\n"))
  }
  message(paste0(100 * CI_level,
             "% credible interval: (", 
             format(signif(results$CI[1], 3), 
                    scientific = FALSE),
             ", ",
             format(signif(results$CI[2], 3), 
                    scientific = FALSE),
             ")\n\n"))
  message(paste0("Probability that the odds ratio is in the ROPE, defined to be (",
             format(signif(ROPE[1], 3), 
                    scientific = FALSE),
             ",",
             format(signif(ROPE[2], 3), 
                    scientific = FALSE),
             ") = ",
             format(signif(results$Pr_oddsratio_in_ROPE, 3), 
                    scientific = FALSE),
             "\n\n")) 
  message("\n----------\n\n")
  
  
  invisible(results)
}