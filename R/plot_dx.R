#' @rdname plot_dx
#' 
#' @title Diagnostic plots for Bayesian regression objects
#' 
#' @param x object of class \code{aov_b}, \code{lm_b}, or \code{glm_b}
#' @param statistic, statistic_m, statistic_y Statistic used to compute 
#' Bayesian p-value (\code{statistic_m} and \code{statistic_y} used for the 
#' mediator and outcome model for a \code{mediate_b} object). 
#' Either "deviance", or else a function taking in data, expected value, and 
#' if applicable to the family, disperion (residual variance for \code{gaussian},
#' and \eqn{\phi} for \code{negbinom} where \eqn{Var(y) = \mu + \mu^2/\phi}).
#' @param mc_error The number of posterior draws will ensure that with 
#' 99% probability the estimated Bayesian p-value will be within 
#' \eqn{\pm} \code{mc_error} of the actual Bayesian p-value.
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
                        statistic = "deviance",
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
  
  bpval = 
    bayes_pvalue(x,
                 statistic = statistic,
                 mc_error = mc_error,
                 seed = seed)
  
  
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
    ggtitle(paste0("Bayesian p-value based on deviance = ",
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
                         statistic = "deviance",
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
  
  bpval = 
    bayes_pvalue(x,
                 statistic = statistic,
                 mc_error = mc_error,
                 seed = seed)
  
  
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
    ggtitle(paste0("Bayesian p-value based on deviance = ",
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
                             statistic_m = "deviance",
                             statistic_y = "deviance",
                             mc_error = 0.005,
                             seed = 1,
                             return_as_list = TRUE,
                             ...){
  
  plot_list = list()
  
  # Mediator model
  plot_list[[1]] = 
    plot_dx(x$model_m,
            statistic = statistic_m,
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
            statistic = statistic_y,
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
