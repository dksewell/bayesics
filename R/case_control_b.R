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
#'  \log\omega \sim N(\texttt{prior\_mean},\texttt{prior\_sd}^2).
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
#' @param prior_mean numeric.  The prior mean on the log odds ratio.  Defaults 
#' to 0 (i.e., odds ratio of 1).
#' @param prior_sd numeric.  The prior sd on the log odds ratio. Defaults to 
#' place 95% prior probability that the odds ratio is between 0.1 and 10.
#' @param plot logical.  Should a plot be shown?
#' @param CI_level The posterior probability to be contained in the 
#' credible interval.
#' @param seed integer.  Always set your seed!!! (ignored if \code{large_sample_approx = TRUE}.)
#' @param mc_error The relative monte carlo error of the quantiles of the CIs. 
#' (ignored if \code{large_sample_approx = TRUE}.)
#' 
#' 
#' @returns An object of class \code{\link{b_procedure-class}}.
#' 
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
  
  ## ---- construct 2x2 table x ----
  
  if (missing(x)) {
    
    if (missing(cases))
      stop("Either `x` or `cases` must be supplied", call. = FALSE)
    
    # cases supplied; interpret its type
    if (is.matrix(cases) || inherits(cases, "table")) {
      x <- as.matrix(cases)
      
    } else if (is.numeric(cases)) {
      
      if (length(cases) != 2)
        stop("`cases` must be a numeric vector of length 2", call. = FALSE)
      
      if (missing(controls))
        stop("Both `cases` and `controls` must be supplied", call. = FALSE)
      
      if (!is.numeric(controls) || length(controls) != 2)
        stop("`controls` must be a numeric vector of length 2", call. = FALSE)
      
      x <- cbind(cases, controls)
      
    } else {
      stop(
        "`cases` must be a numeric vector of length 2, a matrix, or a table",
        call. = FALSE
      )
    }
    
  } else {
    
    # x supplied explicitly
    if (!(is.matrix(x) || inherits(x, "table")))
      stop("`x` must be a 2x2 matrix or table", call. = FALSE)
    
    x <- as.matrix(x)
  }
  
  
  ## ---- validate canonical x ----
  
  if (!is.numeric(x))
    stop("`x` must be numeric", call. = FALSE)
  
  if (!all(dim(x) == c(2L, 2L)))
    stop("`x` must be a 2x2 matrix", call. = FALSE)
  
  if (any(x < 0))
    stop("`x` must contain non-negative counts", call. = FALSE)
  
  if (any(x %% 1 != 0))
    stop("`x` must contain integer counts", call. = FALSE)
  
  
  ## ---- large_sample_approx ----
  
  if (!missing(large_sample_approx)) {
    if (!is.logical(large_sample_approx) ||
        length(large_sample_approx) != 1)
      stop("`large_sample_approx` must be TRUE or FALSE",
           call. = FALSE)
  }
  
  ## ---- ROPE ----
  
  if (!missing(ROPE)) {
    
    if (!is.numeric(ROPE))
      stop("`ROPE` must be numeric",
           call. = FALSE)
    
    if (!(length(ROPE) %in% c(1L, 2L)))
      stop("`ROPE` must be a numeric value or a numeric vector of length 2",
           call. = FALSE)
    
    if (any(ROPE <= 0))
      stop("All values of `ROPE` must be positive",
           call. = FALSE)
    
    if ( (length(ROPE) == 2) &&
         (ROPE[1] >= ROPE[2]) )
      stop("The first element of `ROPE` must be smaller than its second element.",
           call. = FALSE)
  }
  
  ## ---- prior_mean ----
  
  if (!is.numeric(prior_mean) || length(prior_mean) != 1)
    stop("`prior_mean` must be a numeric scalar",
         call. = FALSE)
  
  ## ---- prior_sd ----
  
  if (!is.numeric(prior_sd) || length(prior_sd) != 1 || prior_sd <= 0)
    stop("`prior_sd` must be a positive numeric scalar",
         call. = FALSE)
  
  ## ---- plot ----
  
  if (!is.logical(plot) || length(plot) != 1)
    stop("`plot` must be TRUE or FALSE",
         call. = FALSE)
  
  ## ---- CI_level ----
  
  if (!is.numeric(CI_level) || length(CI_level) != 1)
    stop("`CI_level` must be a numeric scalar",
         call. = FALSE)
  
  if (CI_level <= 0 || CI_level >= 1)
    stop("`CI_level` must be between 0 and 1",
         call. = FALSE)
  
  ## ---- seed ----
  
  if (!is.numeric(seed) || length(seed) != 1 || seed %% 1 != 0)
    stop("`seed` must be a single integer value",
         call. = FALSE)
  
  ## ---- mc_error ----
  
  if (!is.numeric(mc_error) || length(mc_error) != 1 || mc_error <= 0)
    stop("`mc_error` must be a positive numeric scalar",
         call. = FALSE)
  
  
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
  results = 
    list(name = "Case-control analysis",
         data = x,
         print_data = TRUE,
         CI_level = CI_level)
  
  if(large_sample_approx){
    
    ## Set prior
    results$prior =
      paste0("Prior on log odds is: N(",
             format(signif(prior_mean,3),scientific = FALSE),
             ", ",
             format(signif(prior_sd,3),scientific = FALSE),
             "^2)")
    
    ## Get posterior results
    s2 = sum(1.0 / c(x))
    posterior_parameters = 
      c((prior_sd^2 * log_or_hat + s2 * prior_mean) / 
          (prior_sd^2 + s2),
        sqrt(prior_sd^2 * s2 / (prior_sd^2 + s2)))
    names(posterior_parameters) = c("mean","sd")
    
    ## Get point and interval estimates
    results$results = 
      tibble(Quantity = 
               "Odds ratio (at risk vs. not at risk)",
             `Post Mean` = 
               exp(posterior_parameters["mean"] + 
                     0.5 * posterior_parameters["sd"]^2),
             Lower = 
               exp(
                 qnorm(0.5 * alpha_ci,
                       posterior_parameters["mean"],
                       posterior_parameters["sd"])
               ),
             Upper = 
               exp(
                 qnorm(1.0 - 0.5 * alpha_ci,
                       posterior_parameters["mean"],
                       posterior_parameters["sd"])
               ))
    
    
    ## Get ROPE
    results$results = 
      results$results |> 
      mutate(Pr_in_ROPE = 
               pnorm(log(ROPE[2]),
                     posterior_parameters["mean"],
                     posterior_parameters["sd"]) -
               pnorm(log(ROPE[1]),
                     posterior_parameters["mean"],
                     posterior_parameters["sd"]),
             ROPE_lower_bound = ROPE[1],
             ROPE_upper_bound = ROPE[2]
      )
    
    ## Get PDir
    results$pdir = 
      list(pdir = 
             pnorm(0.0,
                   posterior_parameters["mean"],
                   posterior_parameters["sd"])
      )
    results$pdir$description = 
      paste0("Probability that the odds ratio is ",
             ifelse(results$pdir$pdir > 0.5,
                    "less",
                    "greater"),
             " than 1")
    results$pdir$pdir = 
      max(results$pdir$pdir,
          1.0 - results$pdir$pdir)
    
    
    
    # Plot (if requested)
    if(plot){
      results$plot = 
        tibble::tibble(x = seq(qlnorm(0.005,
                                      posterior_parameters["mean"],
                                      posterior_parameters["sd"]),
                               qlnorm(0.995,
                                      posterior_parameters["mean"],
                                      posterior_parameters["sd"]),
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
                                 posterior_parameters["mean"],
                                 posterior_parameters["sd"])
                        },
                      aes(color = "Posterior"),
                      linewidth = 2) + 
        scale_color_manual(values = c("Prior" = "#440154FF", 
                                      "Posterior" = "#FDE725FF")) +
        theme_classic(base_size = 15) +
        xlab("") + 
        ylab("") + 
        labs(color = "Distribution") + 
        ggtitle("Posterior of odds ratio (at risk vs. not at risk)")
    }
    
    
    
  }else{#End: large sample approx
    set.seed(seed)
    message("Cell sizes were too small for large sample approximation.\nInstead, setting uniform prior on Pr(exposure|outcome) and making exact finite sample inference.")
    
    ## Set prior
    results$prior = 
      "Prior on probability of exposure given outcome is: Unif(0,1)"
    
    
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
    
    ## Get point and interval estimates
    results$results = 
      tibble(Quantity = 
               "Odds ratio (at risk vs. not at risk)",
             `Post Mean` = 
               mean(odds_ratios),
             Lower = 
               quantile(odds_ratios,0.5 * alpha_ci),
             Upper = 
               quantile(odds_ratios,1.0 - 0.5 * alpha_ci)
      )
    
    
    ## Get ROPE
    results$results = 
      results$results |> 
      mutate(Pr_in_ROPE = 
               mean( (odds_ratios <= ROPE[2]) & 
                       (odds_ratios >= ROPE[1]) ),
             ROPE_lower_bound = ROPE[1],
             ROPE_upper_bound = ROPE[2]
      )
    
    
    ## Get PDir
    results$pdir = 
      list(pdir = 
             mean(odds_ratios < 1.0)
      )
    results$pdir$description = 
      paste0("Probability that the odds ratio is ",
             ifelse(results$pdir > 0.5,
                    "less",
                    "greater"),
             " than 1")
    results$pdir$pdir = 
      max(results$pdir$pdir,
          1.0 - results$pdir$pdir)
    
    
    ## Plot if requested
    if(plot){
      results$plot = 
        data.frame(or = 
                     odds_ratios[which(odds_ratios <= quantile(odds_ratios,0.99))]) |> 
        ggplot(aes(x = .data$or)) + 
        geom_density(bounds = c(0,Inf),
                     linewidth = 2) + 
        theme_classic(base_size = 15) +
        xlab("") + 
        ylab("") + 
        ggtitle("Posterior of odds ratio (at risk vs. not at risk)")
    }
    
  }#End: small sample inference
  
  results = 
    structure(results,
              class = "b_procedure")
  
  return(results)
}