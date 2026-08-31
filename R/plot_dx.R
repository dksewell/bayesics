#' @rdname plot_dx
#' 
#' @title Diagnostic Plots for Bayesian Regression Objects
#' 
#' @param x object of class \code{aov_b}, \code{lm_b}, \code{glm_b}, or
#' \code{mediate_b}
#' @param statistic, Statistic used to compute Bayesian p-value.
#' If missing, the default statistic will either be the Shapiro-Wilk 
#' test statistic if the family is \code{gaussian} or else the deviance.  
#' User specified functions are allowed, and must take in the response 
#' variable, its expected value, and if applicable to the family, 
#' dispersion (residual variance for \code{gaussian} and \eqn{\phi} 
#' for \code{negbinom}, where \eqn{Var(y) = \mu + \mu^2/\phi}).
#' If x is of class \code{mediate_b}, \code{statistic} should be 
#' a named list with names equal to "m" and "y" for the mediator 
#' and the outcome models respectively.
#' @param mc_error The number of posterior draws will ensure that with 
#' 99% probability the estimated Bayesian p-value will be within 
#' \eqn{\pm} \code{mc_error} of the actual Bayesian p-value. 
#' NOTE: Cross-platform reproducibility may not be 
#' guaranteed due to numerical differences in stats::density() between 
#' operating systems, leading to slightly different numbers of posterior 
#' samples obtained.
#' @param seed integer.
#' @param return_as_list logical.  If TRUE, a list of ggplots will be returned, 
#' rather than a single plot produced by the patchwork package.
#' @param seed integer.
#' @param ... arguments passed on to plot_dx

#' @export
plot_dx = function(x,
                   statistic,
                   mc_error,
                   seed,
                   return_as_list,
                   ...){
  UseMethod("plot_dx")
}


#' @rdname plot_dx
#' @exportS3Method plot_dx lm_b
plot_dx.lm_b = function(x,
                        statistic,
                        mc_error = 0.005,
                        seed = 1,
                        return_as_list = TRUE,
                        ...){
  
  plot_list = list()
  
  # Check normality, heteroscedasticity, etc.
  if(x$family$family == "gaussian"){
    plot_list[["fitted_vs_residuals"]] =
      tibble::tibble(yhat = x$fitted,
                     epsilon = x$standardized_residuals) |>
      ggplot(aes(y = .data$epsilon,x = .data$yhat)) +
      geom_hline(yintercept = 0,
                 linetype = 2,
                 color = "gray35") +
      geom_point(alpha = 0.6) +
      xlab(expression(hat(y))) +
      ylab(expression(hat(epsilon))) +
      theme_classic() +
      ggtitle("Fitted vs. Residuals")
    
    plot_list[["qqnorm"]] =
      tibble::tibble(yhat = x$fitted,
                     epsilon = x$standardized_residuals) |>
      ggplot(aes(sample = .data$epsilon)) +
      geom_qq(alpha = 0.3) +
      geom_qq_line() +
      xlab("Theoretical quantiles") +
      ylab("Empirical quantiles") +
      theme_classic() +
      ggtitle("QQ norm plot")
  }
  
  if(missing(statistic)){
    bpval = 
      bayes_pvalue(x,
                   mc_error = mc_error,
                   seed = seed)
  }else{
    bpval = 
      bayes_pvalue(x,
                   statistic = statistic,
                   mc_error = mc_error,
                   seed = seed)
  }
  
  
  plot_list$bpvals = 
    bpval$statistic_posterior_draws |> 
    dplyr::mutate(obs_gr_pred = .data$T_y_observed  > .data$T_y_predicted) |> 
    ggplot(aes(x = .data$T_y_predicted,
               y = .data$T_y_observed ,
               color = .data$obs_gr_pred)) + 
    geom_point(alpha = 0.05) + 
    geom_abline(intercept = 0,
                slope = 1) + 
    xlab(bquote(T(y[pred] * "," * beta))) +
    ylab(bquote(T(y[obs] * "," * beta))) +
    theme_classic() +
    scale_color_viridis_d() +
    ggtitle(paste0("Bayesian p-value = ",
                   round(bpval$bpvalue,3))) + 
    theme(legend.position = "none")
  
  if(return_as_list){
    return(plot_list)
  }else{
    return(
      patchwork::wrap_plots(plot_list)
    )
  }
}

#' @rdname plot_dx
#' @exportS3Method plot_dx aov_b
plot_dx.aov_b = function(x,
                         statistic,
                         mc_error = 0.005,
                         seed = 1,
                         return_as_list = TRUE,
                         ...){
  
  plot_list = list()
  
  if(x$family$family == "gaussian"){
    plot_list[["residuals_by_group"]] =
      tibble::tibble(group = x$data$group,
                     yhat = x$fitted,
                     epsilon = x$residuals) |>
      ggplot(aes(y = .data$epsilon,x = .data$group)) +
      geom_hline(yintercept = 0,
                 linetype = 2,
                 color = "gray35") +
      geom_violin(alpha = 0.6) +
      xlab(all.vars(x$formula)[2]) +
      ylab(expression(hat(epsilon))) +
      theme_classic() +
      ggtitle("Residual plot by group")
    
    plot_list[["qqnorm"]] =
      tibble::tibble(group = x$data$group,
                     yhat = x$fitted,
                     epsilon = x$residuals) |>
      ggplot(aes(sample = .data$epsilon)) +
      geom_qq(alpha = 0.3) +
      geom_qq_line() +
      xlab("Theoretical quantiles") +
      ylab("Empirical quantiles") +
      theme_classic() +
      ggtitle("QQ norm plot")
  }
  
  if(missing(statistic)){
    bpval = 
      bayes_pvalue(x,
                   mc_error = mc_error,
                   seed = seed)
  }else{
    bpval = 
      bayes_pvalue(x,
                   statistic = statistic,
                   mc_error = mc_error,
                   seed = seed)
  }
  
  
  plot_list$bpvals = 
    bpval$statistic_posterior_draws |> 
    dplyr::mutate(obs_gr_pred = .data$T_y_observed  > .data$T_y_predicted) |> 
    ggplot(aes(x = .data$T_y_predicted,
               y = .data$T_y_observed ,
               color = .data$obs_gr_pred)) + 
    geom_point(alpha = 0.05) + 
    geom_abline(intercept = 0,
                slope = 1) + 
    xlab(bquote(T(y[pred] * "," * beta))) +
    ylab(bquote(T(y[obs] * "," * beta))) +
    theme_classic() +
    scale_color_viridis_d() +
    ggtitle(paste0("Bayesian p-value = ",
                   round(bpval$bpvalue,3))) + 
    theme(legend.position = "none")
  
  if(return_as_list){
    return(plot_list)
  }else{
    return(
      patchwork::wrap_plots(plot_list)
    )
  }
}


#' @rdname plot_dx
#' @exportS3Method plot_dx mediate_b
plot_dx.mediate_b = function(x,
                             statistic = list(m = NULL,
                                              y = NULL),
                             mc_error = 0.005,
                             seed = 1,
                             return_as_list = TRUE,
                             ...){
  
  plot_list = list()
  
  if(!is.null(statistic$m) && !("function" %in% class(statistic$m)))
    stop("Is statistic for the mediator model is provided, it must be a function that takes in y, E(y), and if applicable a dispersion parameter, in that order.")
  
  if(!is.null(statistic$y) && !("function" %in% class(statistic$y)))
    stop("Is statistic for the outcome model is provided, it must be a function that takes in y, E(y), and if applicable a dispersion parameter, in that order.")
  
  if(is.null(statistic$m)){
    statistic$m <- function(y, mu, dispersion = NULL) {
      switch(x$model_m$family$family,
             gaussian   = shapiro.test((y - mu) / sqrt(dispersion))$statistic,
             binomial   = -2.0 * sum(dbinom(y,x$model_m$trials,mu/x$model_m$trials,log=T)),
             poisson    = -2.0 * sum(dpois(y,mu,log=T)),
             negbinom   = -2.0 * sum(dnbinom(y,mu = mu,size = dispersion,log=T)),
             stop("Unsupported family")
      )
    }
  }
  
  if(is.null(statistic$y)){
    statistic$y <- function(y, mu, dispersion = NULL) {
      switch(x$model_y$family$family,
             gaussian   = shapiro.test((y - mu) / sqrt(dispersion))$statistic,
             binomial   = -2.0 * sum(dbinom(y,x$model_y$trials,mu/x$model_y$trials,log=T)),
             poisson    = -2.0 * sum(dpois(y,mu,log=T)),
             negbinom   = -2.0 * sum(dnbinom(y,mu = mu,size = dispersion,log=T)),
             stop("Unsupported family")
      )
    }
  }
  
  # Mediator model
  plot_list[[1]] = 
    plot_dx(x$model_m,
            statistic = statistic$m,
            mc_error = mc_error,
            seed = seed,
            return_as_list = TRUE)
  for(j in names(plot_list[[1]])){
    plot_list[[1]][[j]] = 
      plot_list[[1]][[j]] +
      ggtitle(paste0(plot_list[[1]][[j]]$labels$title,
                     " (Mediator model)"))
  }
  
  # Outcome model
  plot_list[[2]] = 
    plot_dx(x$model_y,
            statistic = statistic$y,
            mc_error = mc_error,
            seed = seed,
            return_as_list = TRUE)
  for(j in names(plot_list[[2]])){
    plot_list[[2]][[j]] = 
      plot_list[[2]][[j]] +
      ggtitle(paste0(plot_list[[2]][[j]]$labels$title,
                     " (Outcome model)"))
  }
  
  
  
  plot_list = do.call(c,plot_list)
  
  if(return_as_list){
    return(plot_list)
  }else{
    return(
      wrap_plots(plot_list)
    )
  }
}
