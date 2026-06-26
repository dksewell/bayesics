#' @name plot
#' 
#' @title Plots \code{bayesics} Objects.
#' 
#' @param x A \code{bayesics} object
#' @param type character. Select any of "diagnostics", 
#' "cred band", and/or "pred band".  If plotting a 
#' \code{mediate_b} object, the valid values for \code{type} 
#' are "diagnostics" (or "dx"), "acme", or "ade".
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
#' @param variable character. If type = "pdp" , which variable should be plotted?
#' @param exemplar_covariates data.frame or tibble with exactly one row.  
#' Used to fix other covariates while varying the variable of interest for the plot.
#' @param combine_pred_cred logical. If type includes both "cred band" and "pred band", 
#' should the credible band be superimposed on the prediction band or 
#' plotted separately?
#' @param variable_seq_length integer. Number of points used to draw pdp.
#' @param return_as_list logical.  If TRUE, a list of ggplots will be returned, 
#' rather than a single plot produced by the patchwork package.
#' @param CI_level Posterior probability covered by credible interval
#' @param PI_level Posterior probability covered by prediction interval
#' @param backtransformation function.  If a transformation of 
#' the response variable was used, \code{backtransformation} 
#' should be the inverse of this transformation function.  E.g., 
#' if you fit lm_b(log(y) ~ x), then set \code{backtransformation=exp}. 
#' @param n_draws integer.  Number of posterior draws used for visualization 
#' of survival curves.  Ignored if \code{x} is not a \code{survfit_b} object.
#' @param ... optional arguments.
#' 
#' @returns If \code{return_as_list=TRUE}, a list of requested ggplots.
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
#' plot(fit1)
#' }
#' 
#' 
#' @rdname plot
#' @method plot lm_b
#' @export
plot.lm_b = function(x,
                     type = c("diagnostics",
                              "cred band",
                              "pred band"),
                     statistic,
                     mc_error = 0.005,
                     seed = 1,
                     variable,
                     exemplar_covariates,
                     combine_pred_cred = TRUE,
                     variable_seq_length = 30,
                     return_as_list = FALSE,
                     CI_level = 0.95,
                     PI_level = 0.95,
                     backtransformation = function(x){x},
                     ...){
  
  type = match.arg(type,several.ok = TRUE)
  
  plot_list = list()
  
  if(x$model_type == "nonparametric"){
    type = 
      setdiff(type,c("diagnostics",
                     "pred band"))
  }else{
    if( (x$family$family == "binomial") & 
        ("pred band" %in% type) ){
      type = "cred band"
    }
  }
  
  if(length(type) == 0)
    stop("No valid plotting type given.")
  
  if("diagnostics" %in% type){
    if(missing(statistic)){
      plot_list[[1]] = 
        plot_dx(x = x,
                mc_error = mc_error,
                seed = seed,
                return_as_list = TRUE)
    }else{
      plot_list[[1]] = 
        plot_dx(x = x,
                statistic = statistic,
                mc_error = mc_error,
                seed = seed,
                return_as_list = TRUE)
    }
  }
  
  if(length(intersect(type,
                      c("cred band",
                        "pred band"))) != 0){
    
    if(missing(variable) & missing(exemplar_covariates)){
      plot_list[[2]] =
        plot_bands(x,
                   type = setdiff(type,"diagnostics"),
                   combine_pred_cred = combine_pred_cred,
                   variable_seq_length = variable_seq_length,
                   CI_level = CI_level,
                   PI_level = PI_level,
                   backtransformation = backtransformation,
                   return_as_list = TRUE)
    }
    if(missing(variable) & !missing(exemplar_covariates)){
      plot_list[[2]] =
        plot_bands(x,
                   type = setdiff(type,"diagnostics"),
                   exemplar_covariates = exemplar_covariates,
                   combine_pred_cred = combine_pred_cred,
                   variable_seq_length = variable_seq_length,
                   CI_level = CI_level,
                   PI_level = PI_level,
                   backtransformation = backtransformation,
                   return_as_list = TRUE)
    }
    if(!missing(variable) & missing(exemplar_covariates)){
      plot_list[[2]] =
        plot_bands(x,
                   type = setdiff(type,"diagnostics"),
                   variable = variable,
                   combine_pred_cred = combine_pred_cred,
                   variable_seq_length = variable_seq_length,
                   CI_level = CI_level,
                   PI_level = PI_level,
                   backtransformation = backtransformation,
                   return_as_list = TRUE)
    }
    if(!missing(variable) & !missing(exemplar_covariates)){
      plot_list[[2]] =
        plot_bands(x,
                   type = setdiff(type,"diagnostics"),
                   variable = variable,
                   exemplar_covariates = exemplar_covariates,
                   combine_pred_cred = combine_pred_cred,
                   variable_seq_length = variable_seq_length,
                   CI_level = CI_level,
                   PI_level = PI_level,
                   backtransformation = backtransformation,
                   return_as_list = TRUE)
    }
    
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





#' @rdname plot
#' @method plot mediate_b
#' @export
plot.mediate_b = function(x,
                          type = c("diagnostics","acme","ade"),
                          statistic = list(m = NULL,
                                           y = NULL),
                          return_as_list = FALSE,
                          seed = 1,
                          mc_error = 0.005,
                          ...){
  
  type = match.arg(type,several.ok = TRUE)
  
  
  # Start diagnostic plots
  if("diagnostics" %in% type){
    plot_list = 
      plot_dx(x,
              statistic = statistic,
              mc_error = mc_error,
              seed = seed,
      )
  }else{
    plot_list = list()
  }
  
  
  # Start ACME plots
  if("acme" %in% type){
    
    ## Simple case
    if(nrow(x$summary) == 4){
      
      plot_list$acme = 
        x$posterior_draws |> 
        ggplot(aes(x = .data$ACME)) +
        geom_histogram(alpha = 0.5) + 
        theme_classic() +
        ggtitle("Avgerage Causal Mediation Effect")
      
    }else{#End: simple case
      ## Complex case
      plot_list$acme = 
        x$posterior_draws[,c("ACME_control",
                             "ACME_treat")] |> 
        tidyr::pivot_longer(cols = everything(),
                            names_to = "Treatment",
                            names_prefix = "ACME_",
                            values_to = "acme") |> 
        dplyr::mutate(Treatment = 
                        ifelse(.data$Treatment == "treat",
                               x$treat_value,
                               x$control_value) |> 
                        as.character()) |> 
        ggplot() +
        geom_histogram(aes(x = .data$acme,
                           fill = .data$Treatment),
                       alpha = 0.25,
                       position = "identity") + 
        scale_fill_viridis_d() +
        theme_classic() +
        ggtitle("Avgerage Causal Mediation Effect")
      if( ("ade" %in% type) & (!return_as_list)){
        plot_list$acme = 
          plot_list$acme +
          theme(legend.position = "none")
      }else{
        plot_list$acme = 
          plot_list$acme +
          guides(fill = guide_legend(title = "Treatment held\nconstant at..."))
      }
      
    }#End: Complex case
    
  }#End: ACME plots
  
  # Start ADE plots
  if("ade" %in% type){
    ## Simple case
    if(nrow(x$summary) == 4){
      
      plot_list$ade = 
        x$posterior_draws |> 
        ggplot(aes(x = .data$ADE)) +
        geom_histogram(alpha = 0.5) + 
        theme_classic() +
        ggtitle("Average Direct Effect")
      
    }else{#End: simple case
      ## Complex case
      
      plot_list$ade = 
        x$posterior_draws[,c("ADE_control",
                             "ADE_treat")] |> 
        tidyr::pivot_longer(cols = everything(),
                            names_to = "Treatment",
                            names_prefix = "ADE_",
                            values_to = "ade") |> 
        dplyr::mutate(Treatment = 
                        ifelse(.data$Treatment == "treat",
                               x$treat_value,
                               x$control_value) |> 
                        as.character()) |> 
        ggplot() +
        geom_histogram(aes(x = .data$ade,
                           fill = .data$Treatment),
                       alpha = 0.25,
                       position = "identity") + 
        scale_fill_viridis_d() +
        theme_classic() +
        ggtitle("Avgerage Direct Effect") +
        guides(fill = guide_legend(title = "Treatment held\nconstant at..."))
      
    }#End: Complex case
  }#End: ADE plots
  
  
  if(return_as_list){
    return(plot_list)
  }else{
    return(
      wrap_plots(plot_list)
    )
  }
  
}



#' @rdname plot
#' @method plot survfit_b
#' @export
plot.survfit_b = function(x,
                          n_draws = 1e4,
                          seed = 1,
                          CI_level = 0.95,
                          ...){
  
  alpha_ci = 1.0 - CI_level
  
  times = model.response(x$data)[,1]
  
  if(x$single_group_analysis){
    
    set.seed(seed)
    lambda_draws = 
      sapply(1:nrow(x$intervals),
             function(j){
               rgamma(n_draws,
                      shape = x$posterior_parameters[j,1],
                      rate = x$posterior_parameters[j,2])
             })
    
    intwidths = 
      x$intervals[,2] -
      x$intervals[,1]
    
    if(length(unique(times)) > 250){
      t_seq = 
        seq(.Machine$double.eps,max(times),
            l = 200)
    }else{
      t_seq = 
        c(.Machine$double.eps,unique(times))
    }
    
    j_of_t = 
      sapply(t_seq,function(s){
        max(which(c(-.Machine$double.eps,
                    x$intervals[,2]) < s))
      })
    
    lambda_intwidth = 
      lambda_draws[,-ncol(lambda_draws),drop=FALSE]
    for(j in 1:(ncol(lambda_draws)-1)){
      lambda_intwidth[,j] = 
        lambda_intwidth[,j] * intwidths[j]
    }
    if(ncol(lambda_intwidth) == 1){
      lambda_intwidth_cumsums = 
        cbind(0.0,
              apply(lambda_intwidth,
                    1,
                    cumsum)
        )
    }else{
      lambda_intwidth_cumsums = 
        cbind(0.0,
          apply(lambda_intwidth,
                1,
                cumsum) |> 
          t()
        )
    }
      
    plotting_df = 
      tibble::tibble(Time = t_seq,
                     `S(t)` = 0.0,
                     Lower = 0.0,
                     Upper = 0.0)
    for(tt in 1:length(t_seq)){
      S_t_draws = 
        exp(-lambda_intwidth_cumsums[,j_of_t[tt]] -
              lambda_draws[,j_of_t[tt]] * (t_seq[tt] - x$intervals[j_of_t[tt],1])
        )
      plotting_df$`S(t)`[tt] = 
        mean(S_t_draws)
      plotting_df$Lower[tt] = 
        quantile(S_t_draws,
                 0.5 * alpha_ci)
      plotting_df$Upper[tt] = 
        quantile(S_t_draws,
                 1.0 - 0.5 * alpha_ci)
    }
    
    
    survplot = 
      plotting_df |> 
      ggplot(aes(x = .data$Time)) + 
      geom_ribbon(aes(ymin = .data$Lower,
                      ymax = .data$Upper),
                  fill = "lightsteelblue3",
                  alpha = 0.5) +
      geom_line(aes(y = .data$`S(t)`)) + 
      theme_classic()
    
    print(survplot)
    
    invisible(list(plot = survplot,
                   data = plotting_df))
    
  }else{#End: single group analysis
    
    set.seed(seed)
    G = length(x$group_names)
    plotting_df = list()
    for(g in 1:G){
      lambda_draws =
        sapply(1:nrow(x[[g]]$intervals),
               function(j){
                 rgamma(n_draws,
                        shape = x[[g]]$posterior_parameters[j,1],
                        rate = x[[g]]$posterior_parameters[j,2])
               })
      
      intwidths = 
        x[[g]]$intervals[,2] -
        x[[g]]$intervals[,1]
      
      if(length(unique(times)) > 250){
        t_seq = 
          seq(.Machine$double.eps,max(times),
              l = 200)
      }else{
        t_seq = 
          c(.Machine$double.eps,unique(times))
      }
      
      j_of_t = 
        sapply(t_seq,function(s){
          max(which(c(-.Machine$double.eps,
                      x[[g]]$intervals[,2]) < s))
        })
      
      lambda_intwidth = 
        lambda_draws[,-ncol(lambda_draws),drop=FALSE]
      for(j in 1:(ncol(lambda_draws)-1)){
        lambda_intwidth[,j] = 
          lambda_intwidth[,j] * intwidths[j]
      }
      if(ncol(lambda_intwidth) == 1){
        lambda_intwidth_cumsums = 
          cbind(0.0,
                apply(lambda_intwidth,
                      1,
                      cumsum)
          )
      }else{
        lambda_intwidth_cumsums = 
          cbind(0.0,
                apply(lambda_intwidth,
                      1,
                      cumsum) |> 
                  t()
          )
      }
      
      plotting_df[[g]] = 
        tibble::tibble(Time = t_seq,
                       `S(t)` = 0.0,
                       Lower = 0.0,
                       Upper = 0.0,
                       Group = x$group_names[g])
      for(tt in 1:length(t_seq)){
        S_t_draws = 
          exp(-lambda_intwidth_cumsums[,j_of_t[tt]] -
                lambda_draws[,j_of_t[tt]] * (t_seq[tt] - x[[g]]$intervals[j_of_t[tt],1])
          )
        plotting_df[[g]]$`S(t)`[tt] = 
          mean(S_t_draws)
        plotting_df[[g]]$Lower[tt] = 
          quantile(S_t_draws,
                   0.5 * alpha_ci)
        plotting_df[[g]]$Upper[tt] = 
          quantile(S_t_draws,
                   1.0 - 0.5 * alpha_ci)
      }
      
    }
    
    plotting_df = 
      do.call(dplyr::bind_rows,
              plotting_df)
    
    survplot =
      plotting_df |> 
      ggplot(aes(x = .data$Time)) + 
      geom_ribbon(aes(ymin = .data$Lower,
                      ymax = .data$Upper,
                      fill = .data$Group),
                  alpha = 0.25,
                  color = NA) +
      geom_line(aes(y = .data$`S(t)`,
                    color = .data$Group)) + 
      scale_fill_viridis_d() + 
      scale_color_viridis_d() + 
      theme_classic()
    
    print(survplot)
    
    invisible(list(plot = survplot,
                   data = plotting_df))
    
    
  }#End: multiple group analysis
}


#' @rdname plot
#' @method plot b_procedure
#' @export
plot.b_procedure = function(x,...){
  if(!("plot" %in% names(x))){
    warning("Either no plot is available or the plot argument was set to FALSE.")
  }else{
    print(x$plot)
  }
}
                          