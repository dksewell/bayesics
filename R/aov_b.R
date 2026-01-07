#' Analysis of Variance using Bayesian methods
#' 
#' @details
#' 
#' \strong{MODEL:}
#' The likelihood model is given by 
#' \deqn{
#'  y_{gi} \overset{iid}{\sim} N(\mu_g,\sigma^2_g),
#' }
#' (although if \code{heteroscedastic} is set to FALSE, \eqn{\sigma^2_g=\sigma^2_h} 
#' \eqn{\forall g,h}). 
#' 
#' The prior is given by 
#' \deqn{
#'  \mu_g|\sigma^2_g \overset{iid}{\sim} N\left(\mu,\frac{\sigma^2_g}{\nu}\right), \\
#'  \sigma^2_g \overset{iid}{\sim} \Gamma^{-1}(a/2,b/2),
#' }
#' where \eqn{mu} is set by \code{prior_mean_mu}, \eqn{nu} is set by 
#' \code{prior_mean_nu}, \eqn{a} is set by \code{prior_var_shape}, and 
#' \eqn{b} is set by \code{prior_var_rate}.
#' 
#' The posterior is
#' \deqn{
#'  \mu_g|y,\sigma^2_g \overset{iid}{\sim} N\left(\hat\mu_g,\frac{\sigma^2_g}{\nu_g}\right), \\
#'  \sigma^2_g|y \overset{iid}{\sim} \Gamma^{-1}(a_g/2,b_g/2),
#' }
#' where \eqn{\hat\mu_g}, \eqn{\nu_g}, \eqn{a_g}, and \eqn{b_g} are all returned 
#' by \code{aov_b} in the named element \code{posterior_parameters}.  
#' 
#' \strong{ROPE:}
#' 
#' If missing, the ROPE bounds will be given under the principle of "half of a 
#' small effect size."  Using Cohen's D of 0.2 as a small effect size, the ROPE 
#' is defined in terms of \eqn{-0.1 <} Cohen's D \eqn{ < 0.1}.
#' 
#' 
#' @param formula A formula specifying the model.
#' @param data A data frame in which the variables specified in the formula 
#' will be found. If missing, the variables are searched for in the standard way.
#' @param heteroscedastic logical.  Set to FALSE to assume all groups have 
#' equal variance.
#' @param prior_mean_mu numeric. Hyperparameter for the a priori mean of the 
#' group means.
#' @param prior_mean_nu numeric. Hyperparameter which scales the precision of 
#' the group means.
#' @param prior_var_shape numeric. Twice the shape parameter for the inverse gamma prior on
#' the residual variance(s).  I.e., \eqn{\sigma^2\sim IG}(prior_var_shape/2,prior_var_rate/2).
#' @param prior_var_rate  numeric. Twice the rate parameter for the inverse gamma prior on
#' the residual variance(s).  I.e., \eqn{\sigma^2\sim IG}(prior_var_shape/2,prior_var_rate/2).
#' @param CI_level numeric. Credible interval level.
#' @param ROPE numeric.  Used to compute posterior probability that Cohen's D +/- ROPE
#' @param contrasts numeric/matrix. Either vector of length equal to the number of 
#' levels in the grouping variable, or else a matrix where each row is a separate 
#' contrast, and the number of columns match the number of levels in the grouping 
#' variable.
#' @param improper logical.  Should we use an improper prior that is proportional 
#' to the inverse of the variance?
#' @param seed integer.  Always set your seed!!!
#' @param mc_error The number of posterior draws will ensure that with 99% 
#' probability the bounds of the credible intervals will be within \eqn{\pm} 
#' \code{mc_error}\eqn{\times 4s_y}, that is, within 100\code{mc_error}% of the 
#' trimmed range of y.
#' 
#' @return Object of class "aov_b" with the following elements:
#' \itemize{
#'  \item summary - tibble giving the summary of the model parameters
#'  \item BF_for_different_vs_same_means - Bayes factor in favor of the full 
#'  model (each group has their own mean) vs. the null model (all groups have 
#'  the same mean).
#'  \item pairwise summary - tibble giving the summary comparing all 
#'  factor level means
#'  \item contrasts (if provided) - list with named elements L (the contrasts provided 
#'  by the user) and summary.
#'  \item posterior_draws - mcmc object (see coda package) giving the posterior draws
#'  \item posterior_parameters - 
#'  \itemize{
#'    \item mu_g - the post. means of the group means
#'    \item nu_g - the post. scalars of the precision
#'    \item a_g - (twice) the post. shape of the inv. gamma for the group variances
#'    \item b_g - (twice) the post. rate of the inv. gamma for the group variances.
#'  }
#'  \item hyperparameters - 
#'  \itemize{
#'    \item mu - the prior mean of the group means
#'    \item nu - the prior scalar of the precision
#'    \item a - (twice) the prior shape of the inv. gamma for the group variances
#'    \item b - (twice) the prior rate of the inv. gamma for the group variances.
#'  }
#'  \item \code{mc_error} - absolute errors used to determine number of 
#'  posterior draws for accurate interval estimation
#'  \item other inputs by user
#' }
#' 
#' @references 
#' Charles R. Doss, James M. Flegal, Galin L. Jones, Ronald C. Neath "Markov chain Monte Carlo estimation of quantiles," Electronic Journal of Statistics, Electron. J. Statist. 8(2), 2448-2478, (2014)
#' 
#' @export

aov_b = function(formula,
                 data,
                 heteroscedastic = TRUE,
                 prior_mean_mu,
                 prior_mean_nu = 0.001,
                 prior_var_shape = 0.001,
                 prior_var_rate = 0.001,
                 CI_level = 0.95,
                 ROPE = 0.1,
                 contrasts,
                 improper = FALSE,
                 seed = 1,
                 mc_error = 0.002,
                 compute_bayes_factor = TRUE){
  
  # Set alpha lv
  a = 1 - CI_level
  
  # Extract variable names
  mf = model.frame(formula,data)
  data$y = 
    model.response(mf)
  data$group = mf[[2]]
  data = data[,c("y","group")]
  if(!is.factor(data$group))
    data$group = factor(data$group)
  
  # variables = all.vars(formula)
  # if(!("factor" %in% class(data[[variables[2]]]))){
  #   data[[variables[2]]] = 
  #     factor(data[[variables[2]]])
  # }
  # data =
  #   data |> 
  #   dplyr::rename(y = variables[1],
  #                 group = variables[2])
  
  # Determine if data are too big to compute BF exactly
  if(nrow(data) > 11180){ # This requires ~ 1 GB to create the covariance matrices
    if(!isTRUE(compute_bayes_factor == "force"))
      compute_bayes_factor = FALSE
  }else{
    compute_bayes_factor = 
      ifelse(isTRUE(compute_bayes_factor == "force"),
             TRUE,compute_bayes_factor)
  }
  
  # Check if improper prior \propto 1/\sigma^2 is requested
  if(improper){
    prior_mean_mu = 0.0
    prior_mean_nu = 0.0
    prior_var_shape = -1.0
    prior_var_rate = 0.0
  }else{
    if(missing(prior_mean_mu))
      prior_mean_mu = mean(data$y)
  }
  
  # Get summary stats
  data_quants = 
    data |> 
    dplyr::group_by(.data$group) |> 
    dplyr::summarize(n = dplyr::n(),
                     ybar = mean(.data$y),
                     y2 = sum(.data$y^2),
                     sample_var = var(.data$y)) |> 
    dplyr::mutate(s2 = (n - 1) / n * .data$sample_var)
  
  
  if(heteroscedastic){
    
    # Get posterior parameters
    nu_g = 
      prior_mean_nu + data_quants$n
    mu_g =
      (prior_mean_nu * prior_mean_mu + data_quants$n * data_quants$ybar) /
      nu_g
    a_g =
      prior_var_shape + data_quants$n
    b_g =
      prior_var_rate +
      data_quants$n * data_quants$s2 +
      prior_mean_nu * data_quants$n / (nu_g + data_quants$n) * (prior_mean_mu - data_quants$ybar)^2
    G = length(nu_g)
    
    # Return a summary including the posterior mean, credible intervals, probability of direction, and BF
    ret = list()
    ret$summary = 
      tibble::tibble(Variable = 
                       paste(rep(c("Mean","Var"),each = G),
                             rep(all.vars(formula)[2],2*G),
                             rep(levels(data$group),2),
                             sep = " : "),
                     `Post Mean` = c(mu_g, b_g/2 / (a_g/2 - 1.0)),
                     Lower = c(extraDistr::qlst(a/2, 
                                                df = a_g,
                                                mu = mu_g,
                                                sigma = sqrt(b_g / nu_g / a_g)),
                               extraDistr::qinvgamma(a/2, alpha = a_g/2, beta = b_g/2)),
                     Upper = c(extraDistr::qlst(1 - a/2, 
                                                df = a_g,
                                                mu = mu_g,
                                                sigma = sqrt(b_g / nu_g / a_g)),
                               extraDistr::qinvgamma(1 - a/2, alpha = a_g/2, beta = b_g/2)),
                     `Prob Dir` = c(extraDistr::plst(0, 
                                                  df = a_g,
                                                  mu = mu_g,
                                                  sigma = sqrt(b_g / nu_g / a_g)),
                                 rep(NA,G)))
    ret$summary$`Prob Dir` = 
      sapply(ret$summary$`Prob Dir`, function(x) max(x,1-x))
    
    if( (!improper) & compute_bayes_factor){
      X = model.matrix(y ~ group - 1,
                       data = data)
      ret$BF_for_different_vs_same_means = 
        exp(
          mvtnorm::dmvt(data$y,
                        df = prior_var_shape,
                        delta = rep(prior_mean_mu,nrow(data)),
                        sigma = 
                          prior_var_rate /
                          prior_var_shape * 
                          (diag(nrow(data)) + 
                             1.0 / prior_mean_nu * tcrossprod(X)),
                        log = TRUE) - 
            mvtnorm::dmvt(data$y,
                          df = prior_var_shape,
                          delta = rep(prior_mean_mu,nrow(data)),
                          sigma = 
                            prior_var_rate /
                            prior_var_shape * 
                            (diag(nrow(data)) + 
                               matrix(1.0 / prior_mean_nu,
                                      nrow(data),nrow(data))),
                          log = TRUE)
        )
    }
    
    
    
    
    # Get posterior samples 
    ## Get preliminary draws
    s2_g_draws = 
      future.apply::future_sapply(1:G,
                                  function(g){
                                    extraDistr::rinvgamma(500,
                                                          alpha = a_g[g]/2,
                                                          beta = b_g[g]/2)
                                  },
                                  future.seed = seed)
    mu_g_draws = 
      future.apply::future_sapply(1:G,
                                  function(g){
                                    rnorm(500,
                                          mean = mu_g[g],
                                          sd = sqrt(s2_g_draws[,g] / nu_g[g]))
                                  },
                                  future.seed = seed)
    mu_g_draws = 
      cbind(mu_g_draws,
            matrix(0.0,500,choose(ncol(mu_g_draws),2)))
    dummy = ncol(s2_g_draws) + 1
    for(i in 1:(ncol(s2_g_draws) - 1)){
      for(j in (i + 1):ncol(s2_g_draws)){
        mu_g_draws[,dummy] = 
          mu_g_draws[,i] - mu_g_draws[,j]
        dummy = dummy + 1
      }
    }
    fhats = 
      future.apply::future_lapply(1:ncol(mu_g_draws),
                                  function(i){
                                    stats::density(mu_g_draws[,i])
                                  })
    epsilon = mc_error * 4 * sqrt(data_quants$s2)
    n_draws = 
      future.apply::future_sapply(1:ncol(mu_g_draws),
                                  function(i){
                                    0.5 * a * (1.0 - 0.5 * a) *
                                      (
                                        qnorm(0.5 * (1.0 - 0.99)) / 
                                          epsilon /
                                          fhats[[i]]$y[which.min(abs(fhats[[i]]$x - 
                                                                       quantile(mu_g_draws[,i], 0.5 * a)))]
                                      )^2
                                  }) |> 
      max() |> 
      round()
    
    ## Get all required draws
    s2_g_draws = 
      future.apply::future_sapply(1:G,
                                  function(g){
                                    extraDistr::rinvgamma(n_draws,
                                                          alpha = a_g[g]/2,
                                                          beta = b_g[g]/2)
                                  },
                                  future.seed = seed)
    colnames(s2_g_draws) =
      paste("variance",levels(data$group),sep="_")
    
    mu_g_draws = 
      future.apply::future_sapply(1:G,
                                  function(g){
                                    rnorm(n_draws,
                                          mean = mu_g[g],
                                          sd = sqrt(s2_g_draws[,g] / nu_g[g]))
                                  },
                                  future.seed = seed)
    colnames(mu_g_draws) = 
      paste("mean",levels(data$group),sep="_")
    
    ret$posterior_draws = 
      cbind(mu_g_draws,
            s2_g_draws)
    
    # Get pairwise comparisons
    temp = 
      combn(1:length(levels(data$group)),2)
    ret$pairwise_summary = 
      data.frame(Comparison = 
                   apply(combn(levels(data$group),2),
                         2,
                         function(x) paste(x[1],x[2],sep="-")),
                 Estimate = mu_g[temp[1,]] - mu_g[temp[2,]],
                 Lower = 0.0,
                 Upper = 0.0,
                 pdir = 0.0,
                 ROPE = 0.0,
                 EPR = 0.0,
                 EPR_Lower = 0.0,
                 EPR_Upper = 0.0)
    for(i in 1:nrow(ret$pairwise_summary)){
      ## Get CI for D(g,h)
      ret$pairwise_summary[i,c("Lower","Upper")] = 
        quantile(ret$posterior_draws[,temp[1,i]] - 
                   ret$posterior_draws[,temp[2,i]],
                 probs = c(a/2, 1 - a/2))
      
      # Get Pdir
      ret$pairwise_summary$pdir[i] = 
        mean(ret$posterior_draws[,temp[1,i]] - 
               ret$posterior_draws[,temp[2,i]] < 0)
      ret$pairwise_summary$pdir[i] = 
        max(ret$pairwise_summary$pdir[i],
            1.0 - ret$pairwise_summary$pdir[i])
      
      
      # Get ROPE for D(g,h) based on Cohen's D
      ret$pairwise_summary$ROPE[i] = 
        mean(
          abs(mu_g_draws[,temp[1,i]] - 
                mu_g_draws[,temp[2,i]]) / 
            sqrt(
              (data_quants$n[temp[1,i]] * s2_g_draws[,temp[1,i]] + 
                 data_quants$n[temp[2,i]] * s2_g_draws[,temp[2,i]]) /
                sum(data_quants$n[temp[,i]])
            ) <
            ROPE
        )
      
      # Get EPR(g,h)
      epr_temp = 
        pnorm((mu_g_draws[,temp[1,i]] - mu_g_draws[,temp[2,i]]) / 
                sqrt(s2_g_draws[,temp[1,i]] + s2_g_draws[,temp[2,i]]))
      ret$pairwise_summary$EPR[i] = mean(epr_temp)
      ret$pairwise_summary$EPR_Lower[i] = quantile(epr_temp,a/2)
      ret$pairwise_summary$EPR_Upper[i] = quantile(epr_temp,1.0 - a/2)
      
    }
    colnames(ret$pairwise_summary)[5] = 
      "Prob Dir"
    colnames(ret$pairwise_summary)[6] = 
      paste0("ROPE (",ROPE,")")
    colnames(ret$pairwise_summary)[8:9] = 
      paste("EPR",c("Lower","Upper"))
    ret$pairwise_summary = 
      ret$pairwise_summary |> 
      tibble::as_tibble() |> 
      dplyr::rename(`Post Mean` = dplyr::all_of("Estimate"))
    
    
    
    # Compute contrasts if requested
    if(!missing(contrasts)){
      
      if("matrix" %in% class(contrasts)){
        L = tcrossprod(mu_g_draws, contrasts)
        rownames(contrasts) = paste("contrast",1:nrow(contrasts),sep="_")
        colnames(contrasts) = levels(data$group)
      }else{
        L = mu_g_draws %*% contrasts
        contrasts = 
          matrix(contrasts,nrow=1,
                 dimnames = list("contrasts_1",
                                 levels(data$group)))
      }
      
      if(any(!dplyr::near(rowSums(contrasts),0.0)))
        stop("Contrasts must sum to 0.")
      
      ret$contrasts = 
        list(L = contrasts,
             summary = tibble::tibble(contrast = 1:nrow(contrasts),
                                      `Post Mean` = colMeans(L),
                                      Lower = apply(L,2,quantile,probs = a/2),
                                      Upper = apply(L,2,quantile,probs = 1 - a/2)))
      
    }
    
    ret$formula = formula
    ret$data = 
      data |> 
      dplyr::rename(!!all.vars(formula)[1] := dplyr::all_of("y"))
    ret$mc_error = epsilon
    ret$posterior_parameters = 
      list(mu_g = mu_g,
           nu_g = nu_g,
           a_g = a_g,
           b_g = b_g)
    if(improper){
      ret$hyperparameters = NA
    }else{
      ret$hyperparameters = 
        list(mu = prior_mean_mu,
             nu = prior_mean_nu,
             a = prior_var_shape,
             b = prior_var_rate)
    }
    
    # Get fitted values
    temp = 
      dplyr::left_join(data,
                       tibble::tibble(group = levels(data$group),
                                      fitted = mu_g,
                                      sd = sqrt(b_g / (a_g - 1))),
                       by = "group")
    ret$fitted = temp$fitted
    # Get residuals
    ret$residuals = data$y - ret$fitted
    ret$standardized_residuals = 
      ret$residuals / temp$sd
    
    return(structure(ret,
                     class = "aov_b"))
    
  }else{# start homoscedastic approach
    
    # Get posterior parameters
    nu_g = 
      prior_mean_nu + data_quants$n
    mu_g =
      (prior_mean_nu * prior_mean_mu + data_quants$n * data_quants$ybar) /
      nu_g
    a_G =
      prior_var_shape + sum(data_quants$n)
    b_G =
      prior_var_rate +
      sum(
        data_quants$n * data_quants$s2 +
          prior_mean_nu * data_quants$n / (nu_g + data_quants$n) * (prior_mean_mu - data_quants$ybar)^2
      )
    G = length(nu_g)
    
    # Return a summary including the posterior mean, credible intervals, probability of direction, and BF
    ret = list()
    ret$summary = 
      tibble::tibble(Variable = 
                       c(paste(rep("Mean",G),
                               rep(all.vars(formula)[2],G),
                               levels(data$group),
                               sep = " : "),
                         "Var"),
                     `Post Mean` = c(mu_g, b_G/2 / (a_G/2 - 1.0)),
                     Lower = c(extraDistr::qlst(a/2, 
                                                df = a_G,
                                                mu = mu_g,
                                                sigma = sqrt(b_G / nu_g / a_G)),
                               extraDistr::qinvgamma(a/2, alpha = a_G/2, beta = b_G/2)),
                     Upper = c(extraDistr::qlst(1 - a/2, 
                                                df = a_G,
                                                mu = mu_g,
                                                sigma = sqrt(b_G / nu_g / a_G)),
                               extraDistr::qinvgamma(1 - a/2, alpha = a_G/2, beta = b_G/2)),
                     `Prob Dir` = c(extraDistr::plst(0, 
                                                  df = a_G,
                                                  mu = mu_g,
                                                  sigma = sqrt(b_G / nu_g / a_G)),
                                 NA))
    ret$summary$`Prob Dir` = 
      sapply(ret$summary$`Prob Dir`, function(x) max(x,1-x))
    
    if( (!improper) & compute_bayes_factor){
      X = model.matrix(y ~ group - 1,
                       data = data)
      ret$BF_for_different_vs_same_means = 
        exp(
          mvtnorm::dmvt(data$y,
                        df = prior_var_shape,
                        delta = rep(prior_mean_mu,nrow(data)),
                        sigma = 
                          prior_var_rate /
                          prior_var_shape * 
                          (diag(nrow(data)) + 
                             1.0 / prior_mean_nu * tcrossprod(X)),
                        log = TRUE) - 
            mvtnorm::dmvt(data$y,
                          df = prior_var_shape,
                          delta = rep(prior_mean_mu,nrow(data)),
                          sigma = 
                            prior_var_rate /
                            prior_var_shape * 
                            (diag(nrow(data)) + 
                               matrix(1.0 / prior_mean_nu,
                                      nrow(data),nrow(data))),
                          log = TRUE)
        )
    }
    
    # Get posterior samples
    ## Get preliminary draws
    s2_G_draws =
      extraDistr::rinvgamma(500,
                            alpha = a_G/2,
                            beta = b_G/2)
    mu_g_draws = 
      future.apply::future_sapply(1:G,
                                  function(g){
                                    rnorm(500,
                                          mean = mu_g[g],
                                          sd = sqrt(s2_G_draws / nu_g[g]))
                                  },
                                  future.seed = seed)
    mu_g_draws = 
      cbind(mu_g_draws,
            matrix(0.0,500,choose(ncol(mu_g_draws),2)))
    dummy = G + 1
    for(i in 1:(G - 1)){
      for(j in (i + 1):G){
        mu_g_draws[,dummy] = 
          mu_g_draws[,i] - mu_g_draws[,j]
        dummy = dummy + 1
      }
    }
    fhats = 
      future.apply::future_lapply(1:ncol(mu_g_draws),
                                  function(i){
                                    stats::density(mu_g_draws[,i])
                                  })
    epsilon = mc_error * 4 * sqrt(data_quants$s2)
    n_draws = 
      future.apply::future_sapply(1:ncol(mu_g_draws),
                                  function(i){
                                    0.5 * a * (1.0 - 0.5 * a) *
                                      (
                                        qnorm(0.5 * (1.0 - 0.99)) / 
                                          epsilon /
                                          fhats[[i]]$y[which.min(abs(fhats[[i]]$x - 
                                                                       quantile(mu_g_draws[,i], 0.5 * a)))]
                                      )^2
                                  }) |> 
      max() |> 
      round()
    
    ## Get all required draws
    s2_G_draws =
      extraDistr::rinvgamma(n_draws,
                            alpha = a_G/2,
                            beta = b_G/2)
    
    mu_g_draws = 
      future.apply::future_sapply(1:G,
                                  function(g){
                                    rnorm(n_draws,
                                          mean = mu_g[g],
                                          sd = sqrt(s2_G_draws / nu_g[g]))
                                  },
                                  future.seed = seed)
    colnames(mu_g_draws) = 
      paste("mean",levels(data$group),sep="_")
    
    
    ret$posterior_draws = 
      cbind(mu_g_draws,
            Var = s2_G_draws)
    
    
    # Get pairwise comparisons
    temp = 
      combn(1:length(levels(data$group)),2)
    ret$pairwise_summary = 
      data.frame(Comparison = 
                   apply(combn(levels(data$group),2),
                         2,
                         function(x) paste(x[1],x[2],sep="-")),
                 Estimate = mu_g[temp[1,]] - mu_g[temp[2,]],
                 Lower = 
                   extraDistr::qlst(0.5 * a,
                                    df = a_G,
                                    mu = mu_g[temp[1,]] - mu_g[temp[2,]],
                                    sigma = 
                                      sqrt(b_G / a_G * 
                                             (1.0 / nu_g[temp[1,]] + 
                                                1.0 / nu_g[temp[2,]])
                                      )
                   ),
                 Upper = 
                   extraDistr::qlst(1.0 - 0.5 * a,
                                    df = a_G,
                                    mu = mu_g[temp[1,]] - mu_g[temp[2,]],
                                    sigma = 
                                      sqrt(b_G / a_G * 
                                             (1.0 / nu_g[temp[1,]] + 
                                                1.0 / nu_g[temp[2,]])
                                      )
                   ), 
                 pdir = 
                   pmax(
                     extraDistr::plst(0.0,
                                      df = a_G,
                                      mu = mu_g[temp[1,]] - mu_g[temp[2,]],
                                      sigma = 
                                        sqrt(b_G / a_G * 
                                               (1.0 / nu_g[temp[1,]] + 
                                                  1.0 / nu_g[temp[2,]])
                                        )
                     ),
                     extraDistr::plst(0.0,
                                      df = a_G,
                                      mu = mu_g[temp[1,]] - mu_g[temp[2,]],
                                      sigma = 
                                        sqrt(b_G / a_G * 
                                               (1.0 / nu_g[temp[1,]] + 
                                                  1.0 / nu_g[temp[2,]])
                                        ),
                                      lower.tail = FALSE
                     )
                   ),
                 ROPE = 0.0,
                 EPR = 0.0,
                 EPR_Lower = 0.0,
                 EPR_Upper = 0.0)
    for(i in 1:nrow(ret$pairwise_summary)){
      # Get ROPE for D(g,h) based on Cohen's D
      ret$pairwise_summary$ROPE[i] = 
        mean(
          abs(mu_g_draws[,temp[1,i]] - 
                mu_g_draws[,temp[2,i]]) / 
            sqrt(s2_G_draws) <
            ROPE
        )
      
      # Get EPR(g,h)
      epr_temp = 
        pnorm((mu_g_draws[,temp[1,i]] - mu_g_draws[,temp[2,i]]) / 
                sqrt(2 * s2_G_draws))
      ret$pairwise_summary$EPR[i] = mean(epr_temp)
      ret$pairwise_summary$EPR_Lower[i] = quantile(epr_temp,a/2)
      ret$pairwise_summary$EPR_Upper[i] = quantile(epr_temp,1.0 - a/2)
      
    }
    
    colnames(ret$pairwise_summary)[5] = 
      "Prob Dir"
    colnames(ret$pairwise_summary)[6] = 
      paste0("ROPE (",ROPE,")")
    colnames(ret$pairwise_summary)[8:9] = 
      paste("EPR",c("Lower","Upper"))
    ret$pairwise_summary = 
      ret$pairwise_summary |> 
      tibble::as_tibble() |> 
      dplyr::rename(`Post Mean` =  dplyr::all_of("Estimate"))
    
    # Compute contrasts if requested
    if(!missing(contrasts)){
      
      if("matrix" %in% class(contrasts)){
        L = tcrossprod(mu_g_draws, contrasts)
        rownames(contrasts) = paste("contrast",1:nrow(contrasts),sep="_")
      }else{
        L = mu_g_draws %*% contrasts
        contrasts = 
          matrix(contrasts,nrow=1,dimnames = list("contrasts_1",NULL))
      }
      
      if(any(!dplyr::near(rowSums(contrasts),0.0)))
        stop("Contrasts must sum to 0.")
      
      ret$contrasts = 
        list(L = contrasts,
             summary = tibble::tibble(contrast = 1:nrow(contrasts),
                                      `Post Mean` = colMeans(L),
                                      Lower = apply(L,2,quantile,probs = a/2),
                                      Upper = apply(L,2,quantile,probs = 1 - a/2)))
      
    }
    
    ret$formula = formula
    ret$data = 
      data |> 
      dplyr::rename(!!all.vars(formula)[1] :=  dplyr::all_of("y"))
    ret$mc_error = epsilon
    ret$posterior_parameters = 
      list(mu_g = mu_g,
           nu_g = nu_g,
           a_g = a_G,
           b_g = b_G)
    ret$CI_level = CI_level
    if(improper){
      ret$hyperparameters = NA
    }else{
      ret$hyperparameters = 
        list(mu = prior_mean_mu,
             nu = prior_mean_nu,
             a = prior_var_shape,
             b = prior_var_rate)
    }
    
    # Get fitted values
    temp = 
      dplyr::left_join(data,
                       tibble::tibble(group = levels(data$group),
                                      fitted = mu_g),
                       by = "group")
    ret$fitted = drop(temp$fitted)
    # Get residuals
    ret$residuals = drop(data$y - ret$fitted)
    ret$standardized_residuals = 
      ret$residuals / sqrt(b_G / (a_G - 1))
    
    return(structure(ret,
                     class = "aov_b"))
    
  }
  
}
