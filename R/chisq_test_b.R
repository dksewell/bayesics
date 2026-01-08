#' @name chisq_test_b
#' @aliases independence_b
#' 
#' @title Test of independence for 2-way contingency tables
#' 
#' @param x Either a table or a matrix of counts
#' @param sampling_design Either "multinomial", "fixed rows", or "fixed columns"
#' @param prior Either "jeffreys" (Dirichlet(1/2)) or "uniform" (Dirichlet(1)).  
#' This is ignored if prior_shapes is provided.
#' @param prior_shapes Either a single positive scalar, in which case a 
#' symmetric Dirichlet is used, or else a matrix matching the dimensions 
#' of x or a vector of length \code{prod(dim(x))}.
#' @param CI_level The posterior probability to be contained in the credible 
#' interval.
#' @param ROPE vector of positive values giving ROPE boundaries for each regression. 
#' @param seed Always set your seed!
#' @param mc_error This is the error in probability from the posterior CDF 
#' evaluated at the ROPE bounds. Note that if it is estimated that these
#' probabilities are between 0.11 and 0.89, the more relaxed value of 0.01 is 
#' used.
#' 
#' @details
#' For a 2-way contingency table with R rows and C columns, evaluate 
#' the probability that 
#' \itemize{
#'  \item the joint probabilities \eqn{p_{ij}} are all 
#' within the ROPE of \eqn{p_{i\cdot}\times p_{\cdot j}} for \code{sampling_design = "multinomial"}
#'  \item the probabilities \eqn{p_{j|i}} are all 
#' within the ROPE of \eqn{p_{\cdot j}} if \code{sampling_design = "fixed rows"} or 
#' \code{"fixed columns"}
#' }
#' 
#' @references 
#' 
#' Gunel, Erdogan &  Dickey, James (1974). Bayes factors for independence in contingency tables, Biometrika, 61(3), Pages 545–557, https://doi.org/10.1093/biomet/61.3.545
#' 
#' Kass, R. E., & Raftery, A. E. (1995). Bayes Factors. Journal of the American Statistical Association, 90(430), 773–795.
#' 
#' 
#' @examples
#' # Generate data
#' set.seed(2025)
#' N = 500
#' nR = 5
#' nC = 3
#' dep_probs = 
#'   extraDistr::rdirichlet(1,rep(2,nR*nC)) |> 
#'   matrix(nR,nC)
#' 
#' # Multinomial sampling
#' ## Test independence
#' independence_b(round(N * dep_probs))
#' 
#' ## Use other priors
#' independence_b(round(N * dep_probs),
#'                prior = "uniform")
#' independence_b(round(N * dep_probs),
#'                prior_shapes = 2)
#' independence_b(round(N * dep_probs),
#'                prior_shapes = matrix(1:(nR*nC),nR,nC))
#' 
#' # Fixed marginals
#' independence_b(round(N * dep_probs),
#'                sampling_design = "rows")
#' independence_b(round(N * dep_probs),
#'                sampling_design = "cols")
#' 
#' 
#' @export
independence_b = function(x,
                          sampling_design = "multinomial",
                          ROPE,
                          prior = "jeffreys",
                          prior_shapes,
                          CI_level = 0.95,
                          seed = 1,
                          mc_error = 0.002){
  set.seed(seed)
  alpha_ci = 1.0 - CI_level
  
  sampling_design = 
    c("multinomial",
      "multinomial",
      "fixed rows",
      "fixed rows",
      "fixed columns",
      "fixed columns",
      "fixed columns")[pmatch(tolower(sampling_design),
                              c("poisson",
                                "multinomial",
                                "fixed rows",
                                "rows",
                                "fixed columns",
                                "columns",
                                "cols"))]
  
  if( !("matrix" %in% class(x)) & 
      !("table" %in% class(x)) )
    stop("x must be a table or a matrix.")
  nR = nrow(x)
  nC = ncol(x)
  x = matrix(x,nR,nC)
  if(is.null(colnames(x)))
    colnames(x) = paste("column",1:ncol(x),sep="_")
  if(is.null(rownames(x)))
    rownames(x) = paste("row",1:nrow(x),sep="_")
  
  
  # Multinomial sampling design
  if(sampling_design == "multinomial"){
    
    ## Get prior hyperparameters
    if(missing(prior_shapes)){
      prior = c("uniform",
                "jeffreys")[pmatch(tolower(prior),
                                   c("uniform",
                                     "jeffreys"))]
      
      if(prior == "uniform"){
        message("Prior shape parameters were not supplied.\nA uniform prior will be used.")
        prior_shapes = rep(1.0,nR * nC)
      }
      if(prior == "jeffreys"){
        message("Prior shape parameters were not supplied.\nJeffrey's prior will be used.")
        prior_shapes = rep(0.5,nR * nC)
      }
    }else{
      if(any(prior_shapes <= 0))
        stop("Prior shape parameters must be positive.")
      prior_shapes = c(prior_shapes)
      if("matrix" %in% class(prior_shapes))
        prior_shapes = c(prior_shapes)
      if( !(length(prior_shapes) %in% c(1,nR * nC)) )
        stop("The length of prior_shapes must equal either 1 or the product of the dimension of the table x.")
      if(length(prior_shapes) == 1)
        prior_shapes = rep(prior_shapes,nR * nC)
    }
    
    ## Get ROPE
    if(missing(ROPE)){
      ROPE = c(1.0 / 1.125, 1.125)
      # From Kruchke (2018) on rate ratios from FDA <1.25. (Use half of small effect size for ROPE, hence 0.25/2) 
      #   Use the same thing for odds ratios.
    }else{
      if(length(ROPE) > 2) stop("ROPE must be given as an upper bound, or given as both lower and upper bounds.")
      if((length(ROPE) > 1) & (ROPE[1] >= ROPE[2])) stop("ROPE lower bound must be smaller than ROPE upper bound")
      if(length(ROPE) == 1) ROPE = c(1.0 / ROPE, ROPE)
    }
    
    
    
    ## Get posterior summary statistics
    results = list()
    results$posterior_shapes = 
      x + matrix(prior_shapes,nR,nC)
    
    results$posterior_mean = 
      results$posterior_shapes / sum(results$posterior_shapes)
    
    results$lower_bound = 
      matrix(qbeta(0.5 * alpha_ci,
                   c(results$posterior_shapes),
                   sum(results$posterior_shapes) - c(results$posterior_shapes)),
             nR,nC,
             dimnames = dimnames(x))
    
    results$upper_bound = 
      matrix(qbeta(1.0 - 0.5 * alpha_ci,
                   c(results$posterior_shapes),
                   sum(results$posterior_shapes) - c(results$posterior_shapes)),
             nR,nC,
             dimnames = dimnames(x))
    
    ## Compute CIs and ROPE for independence
    ### Get preliminary samples
    p_draws = 
      extraDistr::rdirichlet(500,c(results$posterior_shapes)) |> 
      array(c(500,nR,nC))
    p_product = 
      future.apply::future_lapply(1:500,
                                  function(n){
                                    tcrossprod(rowSums(p_draws[n,,]),colSums(p_draws[n,,]))
                                  })
    get_odds = function(p) p / (1.0 - p)
    odds_ratios = 
      future.apply::future_lapply(1:500,
                                  function(n){
                                    get_odds(p_draws[n,,]) / 
                                      get_odds(p_product[[n]])
                                  }) |> 
      unlist() |> 
      array(c(nR,nC,500))
    ### Compute the number of draws needed
    #     Use Raftery and Lewis, focusing just on the ROPE region, rather 
    #     than quantiles of a CI.
    #     However, we only need a high degree of precision for values close to
    #     0 or 1.  Adaptively lower precision needed for when quantiles are 
    #     between 0.1 and 0.9.
    get_adaptive_mc_error = function(qq,
                                     mc_error_baseline = mc_error,
                                     mc_error_max = 0.01){
      qtilde = sapply(qq,function(z) min(z, 1.0 - z))
      
      ifelse(qtilde <= 0.1,
             mc_error_baseline,
             sapply(qtilde,function(z)min(mc_error_max,
                                          max(mc_error_baseline,
                                              z - 0.1)
             )
             )
      )
    }
    n_draws = 
      odds_ratios |> 
      future.apply::future_apply(1:2,
                                 function(x){
                                   v1 = mean(x < ROPE[1])
                                   v1 = v1 * (1.0 - v1)
                                   v2 = mean(x < ROPE[2])
                                   v2 = v2 * (1.0 - v2)
                                   adaptive_mc_error = 
                                     get_adaptive_mc_error(c(v1,v2))
                                   qnorm(0.5 * (1.99))^2 *
                                     max( c(v1,v2) / adaptive_mc_error)^2
                                 }) |> 
      c() |> 
      max() |> 
      round()
    #   Save this code for if we end up wanting to use mc_error on CIs 
    # fhats = 
    #   lapply(1:(nR*nC),
    #          function(i){
    #            density(odds_ratios[1 + (i - 1) %% nR,
    #                                1 + (i - 1) %/% nR,],
    #                    from = 0.0 + .Machine$double.eps)
    #          })
    # n_draws = 
    #   future_sapply(1:length(fhats),
    #                 function(i){
    #                   0.5 * alpha_ci * (1.0 - 0.5 * alpha_ci) *
    #                     (
    #                       qnorm(0.5 * (1.0 - 0.99)) / 
    #                         mc_error /
    #                         fhats[[i]]$y[which.min(abs(fhats[[i]]$x - 
    #                                                      quantile(odds_ratios[1 + (i - 1) %% nR,
    #                                                                           1 + (i - 1) %/% nR,],
    #                                                               0.5 * alpha_ci)))]
    #                     )^2
    #                 }) |> 
    #   max() |> 
    #   round()
    
    ### Get all posterior draws needed
    p_draws = 
      extraDistr::rdirichlet(n_draws,c(results$posterior_shapes)) |> 
      array(c(n_draws,nR,nC))
    p_product = 
      future.apply::future_lapply(1:n_draws,
                                  function(n){
                                    tcrossprod(rowSums(p_draws[n,,]),colSums(p_draws[n,,]))
                                  })
    get_odds = function(p) p / (1.0 - p)
    odds_ratios = 
      future.apply::future_lapply(1:n_draws,
                                  function(n){
                                    get_odds(p_draws[n,,]) / 
                                      get_odds(p_product[[n]])
                                  }) |> 
      unlist() |> 
      array(c(nR,nC,n_draws))
    
    ### Compute ROPE
    results$individual_ROPE = 
      odds_ratios |> 
      apply(1:2,
            function(x){
              mean(x <= ROPE[2]) - 
                mean(x <= ROPE[1])
            })
    dimnames(results$individual_ROPE) = 
      dimnames(x)
    
    odds_ratios_binary_ROPE = 
      (odds_ratios <= ROPE[2]) & 
      (odds_ratios >= ROPE[1])
    results$overall_ROPE = 
      apply(odds_ratios_binary_ROPE,3,all) |> 
      mean()
    
    
    ### Compute PDir
    results$prob_pij_less_than_p_i_times_p_j = 
      odds_ratios |> 
      apply(1:2,
            function(x){
              mean(x <= 1)
            })
    results$prob_direction = 
      apply(results$prob_pij_less_than_p_i_times_p_j,1:2,
            function(x) max(x, 1.0 - x)
      )
    dimnames(results$prob_pij_less_than_p_i_times_p_j) = 
      dimnames(results$prob_direction) = 
      dimnames(x)
    
    
    ### Bayes factor using intrinsic prior
    #   The following commented out code is from https://doi.org/10.1198/jasa.2009.tm08106
    #   Unless I have a mistake in my code or misunderstand their supplementary 
    #   material, this does not work at all in practice.
    # n = sum(x)
    # BF10 = 
    #   lgamma(1.0 + nR * nC) -
    #     lgamma(1.0 + n + nR * nC) +
    #     lgamma(n + nR) + lgamma(1 + nC) -
    #     lgamma(1.0 + nR) - lgamma(1.0 + nC)
    # for(i in 1:nR){
    #   for(j in 1:nC){
    #     x_plus_y = x
    #     x_plus_y[i,j] = x_plus_y[i,j] + 1
    #     BF10 = 
    #       BF10 + 
    #       sum(lgamma(1.0 + x_plus_y)) - 
    #       sum(lgamma(1.0 + rowSums(x_plus_y))) - 
    #       sum(lgamma(1.0 + colSums(x_plus_y)))
    #   }
    # }
    #   Instead, we can use results from Gunel and Dickey, which actually match 
    #   with the rest of the inference
    prior_shapes = matrix(prior_shapes,nR,nC)
    lbeta1 = function(x) sum(lgamma(x)) - lgamma(sum(x)) 
    BF01 = 
      lbeta1( rowSums(x) + rowSums(prior_shapes) - (ncol(x) - 1.0) ) - 
      lbeta1( rowSums(prior_shapes) - (ncol(x) - 1.0) ) +
      lbeta1( colSums(x) + colSums(prior_shapes) - (nrow(x) - 1.0) ) - 
      lbeta1( colSums(prior_shapes) - (nrow(x) - 1.0) ) - 
      lbeta1( c(x) + c(prior_shapes) ) +
      lbeta1( c(prior_shapes) )
    
    
    BF10 = exp(-BF01)
    results$BF_for_dependence_vs_independence = BF10
    bf_max = max(BF10,
                 1.0 / BF10)
    results$BF_evidence =
      ifelse(bf_max <= 3.2,
             "Not worth more than a bare mention",
             ifelse(bf_max <= 10,
                    "Substantial",
                    ifelse(bf_max <= 100,
                           "Strong",
                           "Decisive")))
    
    # Print results
    cat("\n----------\n\n2-way table test for independence using Bayesian techniques\n")
    cat("\n----------\n\n")
    
    cat("Prior used: Dirichlet with shape parameters = \n")
    prior_shapes = 
      matrix(prior_shapes,
             nR,nC,
             dimnames = dimnames(x))
    format(signif(prior_shapes, 3), 
           scientific = FALSE) |> 
      noquote() |> 
      print()
    cat("\n\n")
    
    cat("Posterior mean:\n")
    format(signif(results$posterior_mean, 3), 
           scientific = FALSE) |> 
      noquote() |> 
      print()
    cat("\n\n")
    
    cat(paste0(100 * CI_level,
               "% (marginal) credible intervals: \n"))
    credints = 
      matrix("",nR,nC,dimnames = dimnames(x))
    for(i in 1:nR){
      for(j in 1:nC){
        credints[i,j] = 
          paste0("(",
                 format(signif(results$lower_bound[i,j],3),
                        scientific = FALSE),
                 ", ",
                 format(signif(results$upper_bound[i,j],3),
                        scientific = FALSE),
                 ")")
      }
    }
    credints |> 
      noquote() |> 
      print()
    cat("\n\n")
    
    cat("Probability that p_ij < p_(i.) x p_(.j):\n")
    format(signif(results$prob_pij_less_than_p_i_times_p_j, 3), 
           scientific = FALSE) |> 
      noquote() |> 
      print()
    cat("\n\n")
    
    cat("Probability of direction:\n")
    format(signif(results$prob_direction, 3), 
           scientific = FALSE) |> 
      noquote() |> 
      print()
    cat("\n\n")
    
    cat(paste0("Probability that all odds ratios (unrestricted vs. independence) are in the ROPE, defined to be (",
               format(signif(ROPE[1], 3), 
                      scientific = FALSE),
               ",",
               format(signif(ROPE[2], 3), 
                      scientific = FALSE),
               ") = ",
               format(signif(results$overall_ROPE, 3), 
                      scientific = FALSE),
               "\n\n")) 
    
    cat("The marginal probabilities that each odds ratio is in the ROPE:\n")
    format(signif(results$individual_ROPE, 3), 
           scientific = FALSE) |> 
      noquote() |> 
      print()
    cat("\n\n")
    
    cat(paste0("Bayes factor in favor of dependence: ",
               format(signif(BF10, 3), 
                      scientific = FALSE),
               ";\n      =>Level of evidence: ", 
               results$BF_evidence,
               "\n\n")) 
    
    cat("\n----------\n\n")
    
  }
  if(grepl("fixed",sampling_design)){
    
    ## Get prior hyperparameters
    if(missing(prior_shapes)){
      prior = c("uniform",
                "jeffreys")[pmatch(tolower(prior),
                                   c("uniform",
                                     "jeffreys"))]
      
      if(prior == "uniform"){
        message("Prior shape parameters were not supplied.\nA uniform prior will be used.")
        prior_shapes = matrix(1.0,nR, nC)
      }
      if(prior == "jeffreys"){
        message("Prior shape parameters were not supplied.\nJeffrey's prior will be used.")
        prior_shapes = matrix(0.5,nR, nC)
      }
    }else{
      if(any(prior_shapes <= 0))
        stop("Prior shape parameters must be positive.")
      prior_shapes = c(prior_shapes)
      if(!("matrix" %in% class(prior_shapes)))
        prior_shapes = matrix(prior_shapes,nR,nC)
      if( !(length(prior_shapes) %in% c(1,nR * nC)) )
        stop("The length of prior_shapes must equal either 1 or the product of the dimension of the table x.")
      if(length(prior_shapes) == 1)
        prior_shapes = matrix(prior_shapes,nR, nC)
    }
    
    
    ## Procedure is symmetric, so only code for fixed row sums.  Correct at end.
    if(sampling_design == "fixed columns"){
      x = t(x)
      prior_shapes = t(prior_shapes)
      flipped = TRUE
      nR = nrow(x)
      nC = ncol(x)
    }else{
      flipped = FALSE
    }
    
    ## Get ROPE
    if(missing(ROPE)){
      ROPE = c(1.0 / 1.125, 1.125)
      # From Kruchke (2018) on rate ratios from FDA <1.25. (Use half of small effect size for ROPE, hence 0.25/2) 
      #   Use the same thing for odds ratios.
    }else{
      if(length(ROPE) > 2) stop("ROPE must be given as an upper bound, or given as both lower and upper bounds.")
      if((length(ROPE) > 1) & (ROPE[1] >= ROPE[2])) stop("ROPE lower bound must be smaller than ROPE upper bound")
      if(length(ROPE) == 1) ROPE = c(1.0 / ROPE, ROPE)
    }
    
    
    ## Get posterior summary statistics
    results = list()
    results$posterior_shapes = 
      x + matrix(prior_shapes,nR,nC)
    
    results$posterior_mean = 
      results$posterior_shapes |> 
      apply(1,function(x) x / sum(x)) |> 
      t()
    
    results$lower_bound = 
      matrix(qbeta(0.5 * alpha_ci,
                   c(results$posterior_shapes),
                   c(matrix(rowSums(results$posterior_shapes),nR,nC) - 
                       results$posterior_shapes)),
             nR,nC,
             dimnames = dimnames(x))
    
    results$upper_bound = 
      matrix(qbeta(1.0 - 0.5 * alpha_ci,
                   c(results$posterior_shapes),
                   c(matrix(rowSums(results$posterior_shapes),nR,nC) - 
                       results$posterior_shapes)),
             nR,nC,
             dimnames = dimnames(x))
    
    ## Compute CIs and ROPE for independence
    
    ### Get predetermined row marginals
    p_i_dot = rowSums(x) / sum(x)
    
    ### Get preliminary samples
    #### Get posterior draws of p_j_given_i for each row
    p_j_given_i = array(0.0,c(500,nR,nC))
    for(i in 1:nR){
      p_j_given_i[,i,] = 
        extraDistr::rdirichlet(500,results$posterior_shapes[i,])
    }
    #### From this, compute p_j:
    p_j = 
      future.apply::future_sapply(1:500,
                                  function(n){
                                    colSums(p_i_dot * p_j_given_i[n,,])
                                  })
    #### Then compute odds ratios and go from there.
    get_odds = function(p) p / (1.0 - p)
    odds_ratios = 
      future.apply::future_lapply(1:500,
                                  function(n){
                                    get_odds(p_j_given_i[n,,]) / 
                                      get_odds(matrix(p_j[,n],nR,nC,byrow = TRUE))
                                  }) |> 
      unlist() |> 
      array(c(nR,nC,500))
    #### Compute the number of draws needed
    #     Use Raftery and Lewis, focusing just on the ROPE region, rather 
    #     than quantiles of a CI.
    #     However, we only need a high degree of precision for values close to
    #     0 or 1.  Adaptively lower precision needed for when quantiles are 
    #     between 0.1 and 0.9.
    get_adaptive_mc_error = function(qq,
                                     mc_error_baseline = mc_error,
                                     mc_error_max = 0.01){
      qtilde = sapply(qq,function(z) min(z, 1.0 - z))
      
      ifelse(qtilde <= 0.1,
             mc_error_baseline,
             sapply(qtilde,function(z)min(mc_error_max,
                                          max(mc_error_baseline,
                                              z - 0.1)
             )
             )
      )
    }
    n_draws = 
      odds_ratios |> 
      future.apply::future_apply(1:2,
                                 function(x){
                                   v1 = mean(x < ROPE[1])
                                   v1 = v1 * (1.0 - v1)
                                   v2 = mean(x < ROPE[2])
                                   v2 = v2 * (1.0 - v2)
                                   adaptive_mc_error = 
                                     get_adaptive_mc_error(c(v1,v2))
                                   qnorm(0.5 * (1.99))^2 *
                                     max( c(v1,v2) / adaptive_mc_error)^2
                                 }) |> 
      c() |> 
      max() |> 
      round()
    
    ### Get all posterior draws needed
    #### Get posterior draws of p_j_given_i for each row
    p_j_given_i = array(0.0,c(n_draws,nR,nC))
    for(i in 1:nR){
      p_j_given_i[,i,] = 
        extraDistr::rdirichlet(n_draws,results$posterior_shapes[i,])
    }
    #### From this, compute p_j:
    p_j = 
      future.apply::future_sapply(1:n_draws,
                                  function(n){
                                    colSums(p_i_dot * p_j_given_i[n,,])
                                  })
    #### Then compute odds ratios and go from there.
    get_odds = function(p) p / (1.0 - p)
    odds_ratios = 
      future.apply::future_lapply(1:n_draws,
                                  function(n){
                                    get_odds(p_j_given_i[n,,]) / 
                                      get_odds(matrix(p_j[,n],nR,nC,byrow = TRUE))
                                  }) |> 
      unlist() |> 
      array(c(nR,nC,n_draws))
    
    
    ### Compute ROPE
    results$individual_ROPE = 
      odds_ratios |> 
      apply(1:2,
            function(x){
              mean(x <= ROPE[2]) - 
                mean(x <= ROPE[1])
            })
    dimnames(results$individual_ROPE) = 
      dimnames(x)
    
    odds_ratios_binary_ROPE = 
      (odds_ratios <= ROPE[2]) & 
      (odds_ratios >= ROPE[1])
    results$overall_ROPE = 
      apply(odds_ratios_binary_ROPE,3,all) |> 
      mean()
    
    
    ### Compute PDir
    results$prob_p_j_given_i_less_than_p_j = 
      odds_ratios |> 
      apply(1:2,
            function(x){
              mean(x <= 1)
            })
    results$prob_direction = 
      apply(results$prob_p_j_given_i_less_than_p_j,1:2,
            function(x) max(x, 1.0 - x)
      )
    dimnames(results$prob_p_j_given_i_less_than_p_j) = 
      dimnames(results$prob_direction) = 
      dimnames(x)
    
    
    ### Bayes factor 
    prior_shapes = matrix(prior_shapes,nR,nC)
    lbeta1 = function(x) sum(lgamma(x)) - lgamma(sum(x)) 
    BF01 = 
      lbeta1(colSums(x) + colSums(prior_shapes) - (nrow(x) - 1.0)) - 
      lbeta1(colSums(prior_shapes) - (nrow(x) - 1.0)) -
      lbeta1(rowSums(prior_shapes)) + 
      lbeta1(rowSums(x) + rowSums(prior_shapes)) -
      lbeta1(c(x) + c(prior_shapes)) + 
      lbeta1(c(prior_shapes))
    
    
    BF10 = exp(-BF01)
    results$BF_for_dependence_vs_independence = BF10
    bf_max = max(BF10,
                 1.0 / BF10)
    results$BF_evidence =
      ifelse(bf_max <= 3.2,
             "Not worth more than a bare mention",
             ifelse(bf_max <= 10,
                    "Substantial",
                    ifelse(bf_max <= 100,
                           "Strong",
                           "Decisive")))
    
    ## Flip back if needed
    if(flipped){
      for(j in names(results))
        results[[j]] = t(results[[j]])
      nR = ncol(x)
      nC = nrow(x)
      prior_shapes = t(prior_shapes)
      x = t(x)
    }
    
    
    ## Print results
    cat("\n----------\n\n2-way table test for independence using Bayesian techniques\n")
    cat("\n----------\n\n")
    
    cat("Prior used: Dirichlet with shape parameters = \n")
    prior_shapes = 
      matrix(prior_shapes,
             nR,nC,
             dimnames = dimnames(x))
    format(signif(prior_shapes, 3), 
           scientific = FALSE) |> 
      noquote() |> 
      print()
    cat("\n\n")
    
    cat("Posterior mean:\n")
    format(signif(results$posterior_mean, 3), 
           scientific = FALSE) |> 
      noquote() |> 
      print()
    cat("\n\n")
    
    cat(paste0(100 * CI_level,
               "% (marginal) credible intervals: \n"))
    credints = 
      matrix("",nR,nC,dimnames = dimnames(x))
    for(i in 1:nR){
      for(j in 1:nC){
        credints[i,j] = 
          paste0("(",
                 format(signif(results$lower_bound[i,j],3),
                        scientific = FALSE),
                 ", ",
                 format(signif(results$upper_bound[i,j],3),
                        scientific = FALSE),
                 ")")
      }
    }
    credints |> 
      noquote() |> 
      print()
    cat("\n\n")
    
    cat("Probability that p_ij < p_(i.) x p_(.j):\n")
    format(signif(results$prob_p_j_given_i_less_than_p_j, 3), 
           scientific = FALSE) |> 
      noquote() |> 
      print()
    cat("\n\n")
    
    cat("Probability of direction:\n")
    format(signif(results$prob_direction, 3), 
           scientific = FALSE) |> 
      noquote() |> 
      print()
    cat("\n\n")
    
    cat(paste0("Probability that all odds ratios (unrestricted vs. independence) are in the ROPE, defined to be (",
               format(signif(ROPE[1], 3), 
                      scientific = FALSE),
               ",",
               format(signif(ROPE[2], 3), 
                      scientific = FALSE),
               ") = ",
               format(signif(results$overall_ROPE, 3), 
                      scientific = FALSE),
               "\n\n")) 
    
    cat("The marginal probabilities that each odds ratio is in the ROPE:\n")
    format(signif(results$individual_ROPE, 3), 
           scientific = FALSE) |> 
      noquote() |> 
      print()
    cat("\n\n")
    
    cat(paste0("Bayes factor in favor of dependence: ",
               format(signif(BF10, 3), 
                      scientific = FALSE),
               ";\n      =>Level of evidence: ", 
               results$BF_evidence,
               "\n\n")) 
    
    
    cat("\n----------\n\n")
    
    
  }
  
  invisible(results)
}


# #' @export
# homogeneity_b = function(x,
#                          y,
#                          ROPE,
#                          prior = "uniform",
#                          prior_shapes,
#                          CI_level = 0.95,
#                          seed = 1,
#                          mc_error = 0.01){
#   
#   # Get data
#   x = as.matrix(x)
#   if(dim(x)[2] == 1){
#     x = as.vector(x)
#     if(missing(y))
#       stop("if x is a vector, y must also be supplied")
#     if (length(x) != length(y)) 
#       stop("'x' and 'y' must have the same length")
#     x = rbind(x,c(y))
#   }
#   if(nrow(x) != 2)
#     stop("x must have two rows.")
#   
#   J = ncol(x)
#   alpha_ci = 1.0 - CI_level
#   
#   
#   # Prior distribution
#   if(missing(prior_shapes)){
#     prior = c("uniform",
#               "jeffreys")[pmatch(tolower(prior),
#                                  c("uniform",
#                                    "jeffreys"))]
#     
#     if(prior == "uniform"){
#       message("Prior shape parameters were not supplied.\nA uniform prior will be used.")
#       prior_shapes = rep(1.0,J)
#     }
#     if(prior == "jeffreys"){
#       message("Prior shape parameters were not supplied.\nJeffrey's prior will be used.")
#       prior_shapes = rep(0.5,J)
#     }
#   }else{
#     if(any(prior_shapes <= 0))
#       stop("Prior shape parameters must be positive.")
#     if(length(prior_shapes) == 1)
#       prior_shapes == rep(prior_shapes,J)
#     if(length(prior_shapes) != J)
#       stop("Length of prior_shapes should either match the number of categories or be of length 1")
#   }
#   
#   # Get ROPE
#   if(missing(ROPE)){
#     ROPE = c(1.0 / 1.125, 1.125)
#     # From Kruchke (2018) on rate ratios from FDA <1.25. (Use half of small effect size for ROPE, hence 0.25/2) 
#     #   Use the same thing for odds ratios.
#   }else{
#     if(length(ROPE) > 2) stop("ROPE must be given as an upper bound, or given as both lower and upper bounds.")
#     if((length(ROPE) > 1) & (ROPE[1] >= ROPE[2])) stop("ROPE lower bound must be smaller than ROPE upper bound")
#     if(length(ROPE) == 1) ROPE = c(1.0 / ROPE, ROPE)
#   }
#   
#   # Get posterior parameters
#   post_shapes = 
#     tcrossprod(matrix(1.0,2,1), prior_shapes) + x
#   
#   # Get posterior draws
#   set.seed(seed)
#   ## Get preliminary draws
#   p1_draws = 
#     rdirichlet(500,
#                post_shapes[1,])
#   p2_draws = 
#     rdirichlet(500,
#                post_shapes[2,])
#   ## Use CLT for empirical quantiles:
#   #     A Central Limit Theorem For Empirical Quantiles in the Markov Chain Setting. Peter W. Glynn and Shane G. Henderson
#   #     With prob 0.99 we will be within mc_relative_error of the alpha_ci/2 quantile
#   fhat = 
#     lapply(1:J,
#            function(j){
#              density(p1_draws[,j] - p2_draws[,j],
#                      from = -1.0 + .Machine$double.eps,
#                      to = 1.0 - .Machine$double.eps)
#            })
#     
#   n_draws = 
#     sapply(1:J,
#            function(j){
#              0.5 * alpha_ci * (1.0 - 0.5 * alpha_ci) *
#                (
#                  qnorm(0.5 * (1.0 - 0.99)) / 
#                    mc_relative_error /
#                    quantile(p1_draws[,j] - p2_draws[,j], 0.5 * alpha_ci) /
#                    fhat[[j]]$y[which.min(abs(fhat[[j]]$x - 
#                                           quantile(p1_draws[,j] - p2_draws[,j], 0.5 * alpha_ci)))]
#                )^2 |> 
#                round()
#            })
#     
#   
#   
#   
# }
# 
# 
