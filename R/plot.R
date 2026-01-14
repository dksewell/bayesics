#' @name plot
#' 
#' @title Plots linreg objects.
#' 
#' @param x A linreg object
#' @param type character. Select any of "diagnostics" ("dx" is also allowed),
#'  "pdp" (partial dependence plot), "cred band", and/or "pred band".  
#'  NOTE: the credible and prediction bands only work for numeric 
#'  variables.  If plotting a \code{mediate_b} object, the valid 
#'  values for \code{type} are "diagnostics" (or "dx"), "acme", 
#'  or "ade".
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
#' @param ... optional arguments.
#' 
#' @rdname plot
#' @method plot lm_b
#' @export
plot.lm_b = function(x,
                     type,
                     variable,
                     exemplar_covariates,
                     combine_pred_cred = TRUE,
                     variable_seq_length = 30,
                     return_as_list = FALSE,
                     CI_level = 0.95,
                     PI_level = 0.95,
                     backtransformation = function(x){x},
                     ...){
  
  alpha_ci = 1.0 - CI_level
  alpha_pi = 1.0 - PI_level
  
  
  if(missing(type)){
    type = 
      c("diagnostics",
        #"pdp",
        "cred band",
        "pred band")
  }
  
  type = c("diagnostics",
           "diagnostics",
           "pdp",
           "cred band",
           "pred band")[pmatch(tolower(type),
                             c("diagnostics",
                               "dx",
                               "pdp",
                               "cred band",
                               "pred band"))]
  
  if(missing(variable)){
    variable = 
      terms(x) |> 
      delete.response() |> 
      all.vars() |> 
      unique()
  }
  
  N = nrow(x$data)
  
  plot_list = list()
  
  
  # Diagnostic plots
  if("diagnostics" %in% type){
    
    dx_data = 
      tibble::tibble(yhat = x$fitted,
                     epsilon = x$residuals)
    
    plot_list[["fitted_vs_residuals"]] =
      dx_data |>
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
      dx_data |>
      ggplot(aes(sample = .data$epsilon)) +
      geom_qq(alpha = 0.3) + 
      geom_qq_line() +
      xlab("Theoretical quantiles") +
      ylab("Empirical quantiles") +
      theme_classic() +
      ggtitle("QQ norm plot")
    
  }# End: diagnostics
  
  
  # Get unique values and x sequences for plots
  if( length(intersect(c("pdp","cred band","pred band"),
                       type)) > 0){
    
    x_unique = 
      lapply(variable,
             function(v) unique(x$data[[v]]))
    x_seq = 
      lapply(x_unique,
             function(xvals){
               if(length(xvals) > variable_seq_length){
                 return( 
                   seq(min(xvals),
                       max(xvals),
                       l = variable_seq_length)
                 )
               }else{
                 if(is.numeric(xvals)){
                   return(sort(xvals))
                 }else{
                   if(is.character(xvals)){
                     return(
                       factor(sort(xvals),
                              levels = sort(xvals))
                     )
                   }else{
                     return(xvals)
                   }
                 }
               }
             })
    
    names(x_unique) = 
      names(x_seq) = variable
  }# End: Get unique values and x_seq
  
  
  # Partial Dependence Plots
  if("pdp" %in% type){
    
    for(v in 1:length(variable)){
      
      newdata = 
        tibble::tibble(var_of_interest = x_seq[[v]],
                       y = 0.0)
      for(i in 1:length(x_seq[[v]])){
        temp_preds = 
          predict(x,
                  newdata = 
                    x$data |>
                    dplyr::mutate(!!variable[v] := newdata$var_of_interest[i]),
                  CI_level = CI_level,
                  PI_level = PI_level)
        newdata$y[i] = backtransformation(mean(temp_preds$`Post Mean`))
      }
      
      plot_list[[paste0("pdp_",variable[v])]] = 
        x$data |>
        ggplot(aes(x = .data[[variable[v]]],
                   y = .data[[all.vars(x$formula)[1]]])) + 
        geom_point(alpha = 0.2)
      if(is.numeric(x_seq[[v]])){
        plot_list[[paste0("pdp_",variable[v])]] = 
          plot_list[[paste0("pdp_",variable[v])]] + 
          geom_line(data = newdata,
                    aes(x = .data$var_of_interest,
                        y = .data$y))
      }else{
        plot_list[[paste0("pdp_",variable[v])]] = 
          plot_list[[paste0("pdp_",variable[v])]] + 
          geom_point(data = newdata,
                     aes(x = .data$var_of_interest,
                         y = .data$y),
                     size = 3)
      }
      
      plot_list[[paste0("pdp_",variable[v])]] = 
        plot_list[[paste0("pdp_",variable[v])]] + 
        xlab(variable[v]) + 
        ylab(all.vars(x$formula)[1]) + 
        theme_classic() +
        ggtitle("Partial dependence plot")
      
    }
  }# End: PDP
  
  
  # If drawing CI/PI bands, get reference covariate values and prediction/CIs
  if( ("pred band" %in% type) | ("cred band" %in% type) ){
    
    # Get other covariate values
    if(missing(exemplar_covariates)){
      message("Missing other covariate values in 'exemplar_covariates.'  Using medoid observation instead.")
      desmat = 
        model.matrix(x$formula,
                     x$data) |> 
        scale()
      exemplar_covariates = 
        x$data[cluster::pam(desmat,k=1)$id.med,]
    }
    
    # Get CI and PI values
    newdata = list()
    for(v in variable){
      newdata[[v]] = 
        tibble::tibble(!!v := x_seq[[v]])
      for(j in setdiff(names(exemplar_covariates),v)){
        if(is.character(exemplar_covariates[[j]])){
          newdata[[v]][[j]] = 
            factor(exemplar_covariates[[j]],
                   levels = unique(x$data[[j]]))
        }else{
          newdata[[v]][[j]] = exemplar_covariates[[j]]
        }
      }
      
      newdata[[v]] = 
        predict(x,
                newdata = newdata[[v]],
                CI_level = CI_level,
                PI_level = PI_level)
      
      newdata[[v]] = 
        newdata[[v]] |> 
        dplyr::mutate(dplyr::across(dplyr::all_of(c("Post Mean",
                                                    "PI_lower",
                                                    "PI_upper",
                                                    "CI_lower",
                                                    "CI_upper")),
                                    backtransformation)) # below causes no visible binding for global variable ‘Post Mean’ note 
        # dplyr::mutate(dplyr::across(`Post Mean`:CI_upper,backtransformation))
    }
    
  }# End: Get exemplar and PI/CI
  
  
  # Prediction Band plots
  if("pred band" %in% type){
    
    # Get starter plots if !combine_pred_cred
    for(v in variable){
      plot_name_v = 
        paste0(ifelse((!combine_pred_cred) | !("cred band" %in% type),
                      "pred_band_","band_"),v)
      
      if(is.numeric(x$data[[v]])){
        plot_list[[plot_name_v]] =
          x$data |> 
          ggplot(aes(x = .data[[v]],
                     y = .data[[all.vars(x$formula)[1]]])) +
          geom_point(alpha = 0.2)
      }else{
        plot_list[[plot_name_v]] =
          x$data |> 
          ggplot(aes(x = .data[[v]],
                     y = .data[[all.vars(x$formula)[1]]])) +
          geom_violin(alpha = 0.2)
      }
    }
    
    for(v in variable){
      plot_name_v = 
        paste0(ifelse((!combine_pred_cred) | !("cred band" %in% type),
                      "pred_band_","band_"),v)
      
      if(is.numeric(x_seq[[v]])){
        plot_list[[plot_name_v]] =
          plot_list[[plot_name_v]] +
          geom_ribbon(data = newdata[[v]],
                      aes(ymin = .data$PI_lower,
                          ymax = .data$PI_upper),
                      fill = "lightsteelblue3",
                      alpha = 0.5) +
          geom_line(data = newdata[[v]],
                    aes(x = .data[[v]],
                        y = .data$`Post Mean`))
      }else{
        plot_list[[plot_name_v]] =
          plot_list[[plot_name_v]] +
          geom_errorbar(data = newdata[[v]],
                        aes(x = .data[[v]],
                            ymin = .data$PI_lower,
                            ymax = .data$PI_upper),
                        color = "lightsteelblue3") +
          geom_point(data = newdata[[v]],
                     aes(x = .data[[v]],
                         y = .data$`Post Mean`),
                     size = 3)
      }
      
      
    }
    
    
    
  }
  
  if("cred band" %in% type){
    
    # Get starter plots if !combine_pred_cred
    if( (!combine_pred_cred) | !("pred band" %in% type)){
      for(v in variable){
        if(is.numeric(x$data[[v]])){
          plot_list[[paste0("cred_band_",v)]] =
            x$data |> 
            ggplot(aes(x = .data[[v]],
                       y = .data[[all.vars(x$formula)[1]]])) +
            geom_point(alpha = 0.2)
        }else{
          plot_list[[paste0("cred_band_",v)]] =
            x$data |> 
            ggplot(aes(x = .data[[v]],
                       y = .data[[all.vars(x$formula)[1]]])) +
            geom_violin(alpha = 0.2)
        }
      }
    }
    
    for(v in variable){
      plot_name_v = 
        paste0(ifelse((!combine_pred_cred) | !("pred band" %in% type),
                      "cred_band_","band_"),v)
      
      
      
      if(is.numeric(x_seq[[v]])){
        plot_list[[plot_name_v]] =
          plot_list[[plot_name_v]] +
          geom_ribbon(data = newdata[[v]],
                      aes(ymin = .data$CI_lower,
                          ymax = .data$CI_upper),
                      fill = "steelblue4",
                      alpha = 0.5) +
          geom_line(data = newdata[[v]],
                    aes(x = .data[[v]],
                        y = .data$`Post Mean`))
      }else{
        plot_list[[plot_name_v]] =
          plot_list[[plot_name_v]] +
          geom_errorbar(data = newdata[[v]],
                        aes(x = .data[[v]],
                            ymin = .data$CI_lower,
                            ymax = .data$CI_upper),
                        color = "steelblue4") +
          geom_point(data = newdata[[v]],
                     aes(x = .data[[v]],
                         y = .data$`Post Mean`),
                     size = 3)
      }
    }
    
    
  }
  
  
  # Polish up plots
  if( ("pred band" %in% type) | ("cred band" %in% type) ){
    for(v in variable){
      
      for(j in names(plot_list)[grepl("band",names(plot_list)) & grepl(v,names(plot_list))]){
        plot_list[[j]] =
          plot_list[[j]] +
          theme_classic() +
          ggtitle(
            paste0(
              ifelse(
                grepl("pred_",j),
                paste0("Prediction band for ",v),
                ifelse(grepl("cred_",j),
                       paste0("Credible band for ",v),
                       paste0("Cred. and Pred. bands for ",v)
                )
              )
            )
          )
      }
      
    }
  }
  
  
  if(return_as_list){
    return(plot_list)
  }else{
    return(
      wrap_plots(plot_list)
    )
  }
  
  
}



#' @rdname plot
#' @method plot aov_b
#' @export
plot.aov_b = function(x,
                      type = c("diagnostics",
                               "cred band",
                               "pred band"),
                      combine_pred_cred = TRUE,
                      return_as_list = FALSE,
                      CI_level = 0.95,
                      PI_level = 0.95,
                      ...){
  
  type = c("diagnostics",
           "diagnostics",
           "cred band",
           "pred band")[pmatch(tolower(type),
                             c("diagnostics",
                               "dx",
                               "cred band",
                               "pred band"))]
  
  
  plot_list = list()
  
  # Diagnostic plots
  if("diagnostics" %in% type){
    
    dx_data =
      tibble::tibble(group = x$data$group,
                     yhat = x$fitted,
                     epsilon = x$residuals)
    
    plot_list[["residuals_by_group"]] =
      dx_data |>
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
      dx_data |>
      ggplot(aes(sample = .data$epsilon)) +
      geom_qq(alpha = 0.3) +
      geom_qq_line() +
      xlab("Theoretical quantiles") +
      ylab("Empirical quantiles") +
      theme_classic() +
      ggtitle("QQ norm plot")
    
  }# End: diagnostics
  
  
  # If drawing CI/PI bands, get newdata for prediction/CIs
  if( ("pred band" %in% type) | ("cred band" %in% type) ){
    
    # Get CI and PI values
    newdata =
      predict(x,
              CI_level = CI_level,
              PI_level = PI_level)
    
  }# End: Get newdata
  
  
  # Prediction Band plots
  if("pred band" %in% type){
    
    # Get starter plots
    plot_name_v =
      ifelse((!combine_pred_cred) | !("cred band" %in% type),
             "pred_intervals","intervals")
    
    plot_list[[plot_name_v]] =
      x$data |>
      ggplot(aes(x = .data$group,
                 y = .data[[all.vars(x$formula)[1]]])) +
      geom_violin(alpha = 0.2) +
      geom_errorbar(data = newdata,
                    aes(x = .data[[all.vars(x$formula)[2]]],
                        y = .data$`Post Mean`,
                        ymin = .data$PI_lower,
                        ymax = .data$PI_upper),
                    color = "lightsteelblue3") +
      geom_point(data = newdata,
                 aes(x = .data[[all.vars(x$formula)[2]]],
                     y = .data$`Post Mean`),
                 size = 3)
    
    
    
    
    
  }
  
  if("cred band" %in% type){
    
    # Get starter plots if !combine_pred_cred
    if( (!combine_pred_cred) | !("pred band" %in% type)){
      plot_list[["cred_intervals"]] =
        x$data |>
        ggplot(aes(x = .data$group,
                   y = .data[[all.vars(x$formula)[1]]])) +
        geom_violin(alpha = 0.2)
    }
    
    plot_name_v =
      ifelse((!combine_pred_cred) | !("pred band" %in% type),
             "cred_intervals","intervals")
    
    plot_list[[plot_name_v]] =
      plot_list[[plot_name_v]] +
      geom_errorbar(data = newdata,
                    aes(x = .data[[all.vars(x$formula)[2]]],
                        y = .data$`Post Mean`,
                        ymin = .data$CI_lower,
                        ymax = .data$CI_upper),
                    color = "steelblue4") +
      geom_point(data = newdata,
                 aes(x = .data[[all.vars(x$formula)[2]]],
                     y = .data$`Post Mean`),
                 size = 3)
    
  }
  
  
  # Polish up plots
  if( ("pred band" %in% type) | ("cred band" %in% type) ){
    
    for(j in names(plot_list)[grepl("intervals",names(plot_list))]){
      plot_list[[j]] =
        plot_list[[j]] +
        theme_classic() +
        xlab(all.vars(x$formula)[2]) +
        ggtitle(
          paste0(
            ifelse(
              grepl("pred_",j),
              "Prediction intervals",
              ifelse(grepl("cred_",j),
                     "Credible intervals",
                     "Cred. and Pred. intervals"
              )
            )
          )
        )
    }
    
  }
  
  
  if(return_as_list){
    return(plot_list)
  }else{
    return(
      wrap_plots(plot_list)
    )
  }
  
}





#' @param bayes_pvalues_quantiles ADD description!
#' @param seed ADD description!
#' @rdname plot
#' @method plot lm_b_bma
#' @export
plot.lm_b_bma = function(x,
                         type = c("diagnostics",
                                  "cred band",
                                  "pred band"),
                         variable,
                         exemplar_covariates,
                         combine_pred_cred = TRUE,
                         bayes_pvalues_quantiles = c(0.01,1:19/20,0.99),
                         variable_seq_length = 30,
                         return_as_list = FALSE,
                         CI_level = 0.95,
                         PI_level = 0.95,
                         seed = 1,
                         backtransformation = function(x){x},
                         ...){
  
  alpha_ci = 1.0 - CI_level
  alpha_pi = 1.0 - PI_level
  
  type = c("diagnostics",
           "diagnostics",
           "pdp",
           "cred band",
           "pred band")[pmatch(tolower(type),
                             c("diagnostics",
                               "dx",
                               "pdp",
                               "cred band",
                               "pred band"))]
  
  if(missing(variable)){
    variable = 
      terms(x) |> 
      delete.response() |> 
      all.vars() |> 
      unique()
  }
  N = nrow(x$data)
  
  plot_list = list()
  
  
  # Diagnostic plots
  if("diagnostics" %in% type){
    set.seed(seed)
    
    message("Bayesian p-values measure GOF via \nPr(T(y_obs) - T(y_pred) > 0 | y_obs).\nThus values close to 0.5 are ideal.  Be concerned if values are near 0 or 1.\nThese Bayesian p-values correspond to quantiles of the distribution of y.")
    
    bayes_pvalues_quantiles = sort(bayes_pvalues_quantiles)
    
    preds = predict(x)
    
    T_pred = 
      preds$posterior_draws$ynew |> 
      apply(1,quantile,probs = bayes_pvalues_quantiles)
    
    T_obs = quantile(x$data[[ all.vars(x$formula)[1] ]],
                     bayes_pvalues_quantiles)
    
    bpvals = 
      rowMeans(T_obs - T_pred > 0)
    
    plot_list$bpvals = 
      tibble::tibble(quants = bayes_pvalues_quantiles,
                     bpvals  = bpvals) |> 
      ggplot(aes(x = .data$quants,
                 y = .data$bpvals)) + 
      geom_line() + 
      geom_point() + 
      geom_polygon(data = tibble::tibble(x = c(0,1,1,0,0),
                                         y = c(0,0,0.05,0.05,0)),
                   aes(x=.data$x,y=.data$y),
                   color = NA,
                   fill = "firebrick3",
                   alpha = 0.25) +
      geom_polygon(data = tibble::tibble(x = c(0,1,1,0,0),
                                         y = c(1,1,0.95,0.95,1)),
                   aes(x=.data$x,y=.data$y),
                   color = NA,
                   fill = "firebrick3",
                   alpha = 0.25) + 
      geom_polygon(data = tibble::tibble(x = c(0,1,1,0,0),
                                         y = c(0,0,0.025,0.025,0)),
                   aes(x=.data$x,y=.data$y),
                   color = NA,
                   fill = "firebrick3",
                   alpha = 0.5) +
      geom_polygon(data = tibble::tibble(x = c(0,1,1,0,0),
                                         y = c(1,1,0.975,0.975,1)),
                   aes(x=.data$x,y=.data$y),
                   color = NA,
                   fill = "firebrick3",
                   alpha = 0.5) + 
      xlab("Quantiles of outcome") +
      ylab("Bayesian p-values") + 
      theme_minimal()
    
    rm(preds,T_pred)
    
  }# End: diagnostics
  
  
  # Get unique values and x sequences for plots
  if( length(intersect(c("pdp","cred band","pred band"),
                       type)) > 0){
    
    x_unique = 
      lapply(variable,
             function(v) unique(x$data[[v]]))
    x_seq = 
      lapply(x_unique,
             function(xvals){
               if(length(xvals) > variable_seq_length){
                 return( 
                   seq(min(xvals),
                       max(xvals),
                       l = variable_seq_length)
                 )
               }else{
                 if(is.numeric(xvals)){
                   return(sort(xvals))
                 }else{
                   if(is.character(xvals)){
                     return(
                       factor(sort(xvals),
                              levels = sort(xvals))
                     )
                   }else{
                     return(xvals)
                   }
                 }
               }
             })
    
    names(x_unique) = 
      names(x_seq) = variable
  }# End: Get unique values and x_seq
  
  
  # Partial Dependence Plots
  if("pdp" %in% type){
    
    message("Partial dependence plots typically require long run times.  Plan accordingly.")
    
    for(v in 1:length(variable)){
      
      newdata = 
        tibble::tibble(var_of_interest = x_seq[[v]],
                       y = 0.0)
      for(i in 1:length(x_seq[[v]])){
        temp_preds = 
          predict(x,
                  newdata = 
                    x$data |>
                    dplyr::mutate(!!variable[v] := newdata$var_of_interest[i]))
        newdata$y[i] = backtransformation(mean(temp_preds$newdata$`Post Mean`))
      }
      
      plot_list[[paste0("pdp_",variable[v])]] = 
        x$data |>
        ggplot(aes(x = .data[[variable[v]]],
                   y = .data[[all.vars(x$formula)[1]]])) + 
        geom_point(alpha = 0.2)
      if(is.numeric(x_seq[[v]])){
        plot_list[[paste0("pdp_",variable[v])]] = 
          plot_list[[paste0("pdp_",variable[v])]] + 
          geom_line(data = newdata,
                    aes(x = .data$var_of_interest,
                        y = .data$y))
      }else{
        plot_list[[paste0("pdp_",variable[v])]] = 
          plot_list[[paste0("pdp_",variable[v])]] + 
          geom_point(data = newdata,
                     aes(x = .data$var_of_interest,
                         y = .data$y),
                     size = 3)
      }
      
      plot_list[[paste0("pdp_",variable[v])]] = 
        plot_list[[paste0("pdp_",variable[v])]] + 
        xlab(variable[v]) + 
        ylab(all.vars(x$formula)[1]) + 
        theme_classic() +
        ggtitle("Partial dependence plot")
      
    }
  }# End: PDP
  
  
  # If drawing CI/PI bands, get reference covariate values and prediction/CIs
  if( ("pred band" %in% type) | ("cred band" %in% type) ){
    
    # Get other covariate values
    if(missing(exemplar_covariates)){
      message("Missing other covariate values in 'exemplar_covariates.'  Using medoid observation instead.")
      desmat = 
        model.matrix(x$formula,
                     x$data) |> 
        scale()
      exemplar_covariates = 
        x$data[cluster::pam(desmat,k=1)$id.med,]
    }
    
    # Get CI and PI values
    newdata = list()
    for(v in variable){
      newdata[[v]] = 
        tibble::tibble(!!v := x_seq[[v]])
      for(j in setdiff(names(exemplar_covariates),v)){
        if(is.character(exemplar_covariates[[j]])){
          newdata[[v]][[j]] = 
            factor(exemplar_covariates[[j]],
                   levels = unique(x$data[[j]]))
        }else{
          newdata[[v]][[j]] = exemplar_covariates[[j]]
        }
      }
      
      newdata[[v]] = 
        predict(x,
                newdata = newdata[[v]],
                CI_level = CI_level,
                PI_level = PI_level)
      newdata[[v]]$newdata = 
        newdata[[v]]$newdata |> 
        dplyr::mutate(dplyr::across(dplyr::all_of(c("Post Mean",
                                                    "PI_lower",
                                                    "PI_upper",
                                                    "CI_lower",
                                                    "CI_upper")),
                                    backtransformation))
    }
    
  }# End: Get exemplar and PI/CI
  
  
  # Prediction Band plots
  if("pred band" %in% type){
    
    # Get starter plots if !combine_pred_cred
    for(v in variable){
      plot_name_v = 
        paste0(ifelse((!combine_pred_cred) | !("cred band" %in% type),
                      "pred_band_","band_"),v)
      
      if(is.numeric(x$data[[v]])){
        plot_list[[plot_name_v]] =
          x$data |> 
          ggplot(aes(x = .data[[v]],
                     y = .data[[all.vars(x$formula)[1]]])) +
          geom_point(alpha = 0.2)
      }else{
        plot_list[[plot_name_v]] =
          x$data |> 
          ggplot(aes(x = .data[[v]],
                     y = .data[[all.vars(x$formula)[1]]])) +
          geom_violin(alpha = 0.2)
      }
    }
    
    for(v in variable){
      plot_name_v = 
        paste0(ifelse((!combine_pred_cred) | !("cred band" %in% type),
                      "pred_band_","band_"),v)
      
      if(is.numeric(x_seq[[v]])){
        plot_list[[plot_name_v]] =
          plot_list[[plot_name_v]] +
          geom_ribbon(data = newdata[[v]]$newdata,
                      aes(ymin = .data$PI_lower,
                          ymax = .data$PI_upper),
                      fill = "lightsteelblue3",
                      alpha = 0.5) +
          geom_line(data = newdata[[v]]$newdata,
                    aes(x = .data[[v]],
                        y = .data$`Post Mean`))
      }else{
        plot_list[[plot_name_v]] =
          plot_list[[plot_name_v]] +
          geom_errorbar(data = newdata[[v]]$newdata,
                        aes(x = .data[[v]],
                            ymin = .data$PI_lower,
                            ymax = .data$PI_upper),
                        color = "lightsteelblue3") +
          geom_point(data = newdata[[v]]$newdata,
                     aes(x = .data[[v]],
                         y = .data$`Post Mean`),
                     size = 3)
      }
      
      
    }
    
    
    
  }
  
  if("cred band" %in% type){
    
    # Get starter plots if !combine_pred_cred
    if( (!combine_pred_cred) | !("pred band" %in% type)){
      for(v in variable){
        if(is.numeric(x$data[[v]])){
          plot_list[[paste0("cred_band_",v)]] =
            x$data |> 
            ggplot(aes(x = .data[[v]],
                       y = .data[[all.vars(x$formula)[1]]])) +
            geom_point(alpha = 0.2)
        }else{
          plot_list[[paste0("cred_band_",v)]] =
            x$data |> 
            ggplot(aes(x = .data[[v]],
                       y = .data[[all.vars(x$formula)[1]]])) +
            geom_violin(alpha = 0.2)
        }
      }
    }
    
    for(v in variable){
      plot_name_v = 
        paste0(ifelse((!combine_pred_cred) | !("pred band" %in% type),
                      "cred_band_","band_"),v)
      
      
      
      if(is.numeric(x_seq[[v]])){
        plot_list[[plot_name_v]] =
          plot_list[[plot_name_v]] +
          geom_ribbon(data = newdata[[v]]$newdata,
                      aes(ymin = .data$CI_lower,
                          ymax = .data$CI_upper),
                      fill = "steelblue4",
                      alpha = 0.5) +
          geom_line(data = newdata[[v]]$newdata,
                    aes(x = .data[[v]],
                        y = .data$`Post Mean`))
      }else{
        plot_list[[plot_name_v]] =
          plot_list[[plot_name_v]] +
          geom_errorbar(data = newdata[[v]]$newdata,
                        aes(x = .data[[v]],
                            ymin = .data$CI_lower,
                            ymax = .data$CI_upper),
                        color = "steelblue4") +
          geom_point(data = newdata[[v]]$newdata,
                     aes(x = .data[[v]],
                         y = .data$`Post Mean`),
                     size = 3)
      }
    }
    
    
  }
  
  
  # Polish up plots
  if( ("pred band" %in% type) | ("cred band" %in% type) ){
    for(v in variable){
      
      for(j in names(plot_list)[grepl("band",names(plot_list)) & grepl(v,names(plot_list))]){
        plot_list[[j]] =
          plot_list[[j]] +
          theme_classic() +
          ggtitle(
            paste0(
              ifelse(
                grepl("pred_",j),
                paste0("Prediction band for ",v),
                ifelse(grepl("cred_",j),
                       paste0("Credible band for ",v),
                       paste0("Cred. and Pred. bands for ",v)
                )
              )
            )
          )
      }
      
    }
  }
  
  
  if(return_as_list){
    return(plot_list)
  }else{
    return(
      wrap_plots(plot_list)
    )
  }
  
  
}

#' @rdname plot
#' @method plot glm_b
#' @export
plot.glm_b = function(x,
                      type,
                      variable,
                      exemplar_covariates,
                      combine_pred_cred = TRUE,
                      variable_seq_length = 30,
                      return_as_list = FALSE,
                      CI_level = 0.95,
                      PI_level = 0.95,
                      seed = 1,
                      ...){
  
  alpha_ci = 1.0 - CI_level
  alpha_pi = 1.0 - PI_level
  
  if(missing(type)){
    type = 
      c("diagnostics",
        #"pdp",
        "cred band",
        "pred band")
    if(x$family$family != "binomial") type = c(type,"pred band")
  }
  
  type = c("diagnostics",
           "diagnostics",
           "pdp",
           "cred band",
           "pred band")[pmatch(tolower(type),
                             c("diagnostics",
                               "dx",
                               "pdp",
                               "cred band",
                               "pred band"))]
  
  if( (x$family$family == "binomial") & 
      ("pred band" %in% type) ){
    type = setdiff(type,"pred band")
    if(length(type) == 0){
      warning("Prediction band cannot be supplied for a binomial outcome.\nResults shown will be credible band instead.")
      type = "cred band"
    }
  }
  
  if(missing(variable)){
    variable = 
      terms(x) |> 
      delete.response() |> 
      all.vars() |> 
      unique()
    offset_to_rm = 
      attr(terms(x),"offset")
    if(!is.null(offset_to_rm)){
      variable = variable[-(offset_to_rm - 1)] # - 1 because response was removed
    }
  }
  
  N = nrow(x$data)
  
  plot_list = list()
  
  
  # Diagnostic plots
  if("diagnostics" %in% type){
    
    # Extract 
    mframe = model.frame(x$formula, x$data)
    y = model.response(mframe)
    X = model.matrix(x$formula,x$data)
    os = model.offset(mframe)
    N = nrow(X)
    p = ncol(X)
    if(is.null(os)) os = numeric(N)
    
    
    message("Bayesian p-values measure GOF via \nPr(T(y_obs) - T(y_pred) > 0 | y_obs).\nThus values close to 0.5 are ideal.  Be concerned if values are near 0 or 1.\nThis Bayesian p-value corresponds to the deviance.")
    
    if("posterior_covariance" %in% names(x)){
      
      # Get posterior draws of E(y)
      theta_draws = 
        mvtnorm::rmvnorm(5e3,
                         x$summary$`Post Mean`,
                         x$posterior_covariance)
      yhat_draws = 
        x$trials * 
        x$family$linkinv(os + tcrossprod(X, theta_draws[,1:ncol(X)]))
      
      # Get posterior draws of y
      if(x$family$family == "binomial"){
        y_draws = 
          future.apply::future_sapply(1:nrow(yhat_draws),
                                      function(i){
                                        rbinom(ncol(yhat_draws),
                                               x$trials[i],
                                               yhat_draws[i,])
                                      },
                                      future.seed = seed)
        
        deviances_pred = 
          future.apply::future_sapply(1:nrow(y_draws),
                                      function(draw){
                                        -2.0 * 
                                          sum(dbinom(y_draws[draw,],
                                                     x$trials,
                                                     yhat_draws[,draw],
                                                     log = TRUE))
                                      })
        deviances_obs = 
          future.apply::future_sapply(1:nrow(y_draws),
                                      function(draw){
                                        -2.0 * 
                                          sum(dbinom(y,
                                                     x$trials,
                                                     yhat_draws[,draw],
                                                     log = TRUE))
                                      })
        
        
      }
      if(x$family$family == "poisson"){
        y_draws = 
          future.apply::future_sapply(1:nrow(yhat_draws),
                                      function(i){
                                        rpois(ncol(yhat_draws),yhat_draws[i,])
                                      },
                                      future.seed = seed)
        deviances_pred = 
          future.apply::future_sapply(1:nrow(y_draws),
                                      function(draw){
                                        -2.0 * 
                                          sum(dpois(y_draws[draw,],
                                                    yhat_draws[,draw],
                                                    log = TRUE))
                                      })
        deviances_obs = 
          future.apply::future_sapply(1:nrow(y_draws),
                                      function(draw){
                                        -2.0 * 
                                          sum(dpois(y,
                                                    yhat_draws[,draw],
                                                    log = TRUE))
                                      })
      }
      
      if(x$family$family == "negbinom"){
        y_draws = 
          future.apply::future_sapply(1:nrow(yhat_draws),
                                      function(i){
                                        rnbinom(ncol(yhat_draws),
                                                mu = yhat_draws[i,],
                                                size = exp(theta_draws[,ncol(X) + 1]))
                                      },
                                      future.seed = seed)
        deviances_pred = 
          future.apply::future_sapply(1:nrow(y_draws),
                                      function(draw){
                                        -2.0 * 
                                          sum(dnbinom(y_draws[draw,],
                                                      mu = yhat_draws[,draw],
                                                      size = exp(theta_draws[draw,ncol(X)+1]),
                                                      log = TRUE))
                                      })
        deviances_obs = 
          future.apply::future_sapply(1:nrow(y_draws),
                                      function(draw){
                                        -2.0 * 
                                          sum(dnbinom(y,
                                                      mu = yhat_draws[,draw],
                                                      size = exp(theta_draws[draw,ncol(X)+1]),
                                                      log = TRUE))
                                      })
      }
      
    }else{#End: Getting pvals for large sample approx
      
      # Get posterior draws of E(y)
      yhat_draws = 
        x$trials * 
        x$family$linkinv(os + tcrossprod(X, x$proposal_draws[,1:ncol(X)]))
      
      # Get posterior draws of y
      if(x$family$family == "binomial"){
        y_draws = 
          future.apply::future_sapply(1:nrow(yhat_draws),
                                      function(i){
                                        rbinom(ncol(yhat_draws),
                                               x$trials[i],
                                               yhat_draws[i,])
                                      },
                                      future.seed = seed)
        
        deviances_pred = 
          future.apply::future_sapply(1:nrow(y_draws),
                                      function(draw){
                                        -2.0 * 
                                          sum(dbinom(y_draws[draw,],
                                                     x$trials,
                                                     yhat_draws[,draw],
                                                     log = TRUE))
                                      })
        deviances_obs = 
          future.apply::future_sapply(1:nrow(y_draws),
                                      function(draw){
                                        -2.0 * 
                                          sum(dbinom(y,
                                                     x$trials,
                                                     yhat_draws[,draw],
                                                     log = TRUE))
                                      })
      }
      if(x$family$family == "poisson"){
        y_draws = 
          future.apply::future_sapply(1:nrow(yhat_draws),
                                      function(i){
                                        rpois(ncol(yhat_draws),yhat_draws[i,])
                                      },
                                      future.seed = seed)
        deviances_pred = 
          future.apply::future_sapply(1:nrow(y_draws),
                                      function(draw){
                                        -2.0 * 
                                          sum(dpois(y_draws[draw,],
                                                    yhat_draws[,draw],
                                                    log = TRUE))
                                      })
        deviances_obs = 
          future.apply::future_sapply(1:nrow(y_draws),
                                      function(draw){
                                        -2.0 * 
                                          sum(dpois(y,
                                                    yhat_draws[,draw],
                                                    log = TRUE))
                                      })
        
      }
      if(x$family$family == "negbinom"){
        y_draws = 
          future.apply::future_sapply(1:nrow(yhat_draws),
                                      function(i){
                                        rnbinom(ncol(yhat_draws),
                                                mu = yhat_draws[i,],
                                                size = exp(x$proposal_draws[,ncol(X) + 1]))
                                      },
                                      future.seed = seed)
        deviances_pred = 
          future.apply::future_sapply(1:nrow(y_draws),
                                      function(draw){
                                        -2.0 * 
                                          sum(dnbinom(y_draws[draw,],
                                                      mu = yhat_draws[,draw],
                                                      size = exp(x$proposal_draws[draw,ncol(X)+1]),
                                                      log = TRUE))
                                      })
        deviances_obs = 
          future.apply::future_sapply(1:nrow(y_draws),
                                      function(draw){
                                        -2.0 * 
                                          sum(dnbinom(y,
                                                      mu = yhat_draws[,draw],
                                                      size = exp(x$proposal_draws[draw,ncol(X)+1]),
                                                      log = TRUE))
                                      })
        
      }
      
      resample_index = 
        sample(1:length(deviances_obs),length(deviances_obs),TRUE,x$importance_sampling_weights)
      deviances_obs = deviances_obs[resample_index]
      deviances_pred = deviances_pred[resample_index]
      
    }#End: Getting pvals for IS
    
    
    dx_data = 
      tibble::tibble(T_obs = deviances_obs,
                     T_pred = deviances_pred) |> 
      dplyr::mutate(obs_gr_pred = .data$T_obs > .data$T_pred)
    
    plot_list$bpvals = 
      dx_data |> 
      ggplot(aes(x = .data$T_pred,
                 y = .data$T_obs,
                 color = .data$obs_gr_pred)) + 
      geom_point() + 
      geom_abline(intercept = 0,
                  slope = 1) + 
      xlab(bquote(T(y[pred] * "," * beta))) +
      ylab(bquote(T(y[obs] * "," * beta))) +
      theme_classic() +
      scale_color_viridis_d() +
      ggtitle(paste0("Bayesian p-value based on deviance = ",
                     round(mean(deviances_obs > deviances_pred),3))) + 
      theme(legend.position = "none")
    
  }# End: diagnostics
  
  # Get unique values and x sequences for plots
  if( length(intersect(c("pdp","cred band","pred band"),
                       type)) > 0){
    
    x_unique = 
      lapply(variable,
             function(v) unique(x$data[[v]]))
    x_seq = 
      lapply(x_unique,
             function(xvals){
               if(length(xvals) > variable_seq_length){
                 return( 
                   seq(min(xvals),
                       max(xvals),
                       l = variable_seq_length)
                 )
               }else{
                 if(is.numeric(xvals)){
                   return(sort(xvals))
                 }else{
                   if(is.character(xvals)){
                     return(
                       factor(sort(xvals),
                              levels = sort(xvals))
                     )
                   }else{
                     return(xvals)
                   }
                 }
               }
             })
    
    names(x_unique) = 
      names(x_seq) = variable
  }# End: Get unique values and x_seq
  
  
  # Partial Dependence Plots
  if("pdp" %in% type){
    
    for(v in 1:length(variable)){
      
      newdata = 
        tibble::tibble(var_of_interest = x_seq[[v]],
                       y = 0.0)
      suppressMessages({
        for(i in 1:length(x_seq[[v]])){
          temp_preds = 
            predict(x,
                    newdata = 
                      x$data |>
                      dplyr::mutate(!!variable[v] := newdata$var_of_interest[i]),
                    CI_level = CI_level,
                    PI_level = PI_level)
          newdata$y[i] = mean(temp_preds$`Post Mean`)
        }
      })
      
      
      if(is.numeric(x_seq[[v]])){
        plot_list[[paste0("pdp_",variable[v])]] = 
          x$data |>
          ggplot(aes(x = .data[[variable[v]]],
                     y = .data[[all.vars(x$formula)[1]]])) + 
          geom_point(alpha = 0.2) +
          geom_line(data = newdata,
                    aes(x = .data$var_of_interest,
                        y = .data$y))
      }else{
        if(x$family$family %in% c("poisson","negbinom")){
          plot_list[[paste0("pdp_",variable[v])]] = 
            x$data |>
            ggplot(aes(x = .data[[variable[v]]],
                       y = .data[[all.vars(x$formula)[1]]])) + 
            geom_violin(alpha = 0.2)
        }
        if(x$family$family == "binomial"){
          plot_list[[paste0("pdp_",variable[v])]] = 
            x$data |>
            dplyr::group_by(get(variable[v])) |> 
            dplyr::summarize(prop1 = mean(dplyr::near(.data[[all.vars(x$formula)[1]]], 1))) |> 
            ggplot(aes(x = .data$`get(variable[v])`,
                       y = .data$prop1)) + 
            geom_col(fill="gray70") + 
            ylab(all.vars(x$formula)[1])
        }
        
        plot_list[[paste0("pdp_",variable[v])]] = 
          plot_list[[paste0("pdp_",variable[v])]] + 
          geom_point(data = newdata,
                     aes(x = .data$var_of_interest,
                         y = .data$y),
                     size = 3)
      }
      
      plot_list[[paste0("pdp_",variable[v])]] = 
        plot_list[[paste0("pdp_",variable[v])]] + 
        xlab(variable[v]) + 
        ylab(all.vars(x$formula)[1]) + 
        theme_classic() +
        ggtitle("Partial dependence plot")
      
    }
  }# End: PDP
  
  # If drawing CI/PI bands, get reference covariate values and prediction/CIs
  if( ("pred band" %in% type) | ("cred band" %in% type) ){
    
    # Get other covariate values
    if(missing(exemplar_covariates)){
      message("Missing other covariate values in 'exemplar_covariates.'  Using medoid observation instead.")
      desmat = 
        model.matrix(x$formula,
                     x$data) |> 
        scale()
      exemplar_covariates = 
        x$data[cluster::pam(desmat,k=1)$id.med,]
    }
    
    # Get CI and PI values
    newdata = list()
    for(v in variable){
      newdata[[v]] = 
        tibble::tibble(!!v := x_seq[[v]])
      for(j in setdiff(names(exemplar_covariates),v)){
        if(is.character(exemplar_covariates[[j]])){
          newdata[[v]][[j]] = 
            factor(exemplar_covariates[[j]],
                   levels = unique(x$data[[j]]))
        }else{
          newdata[[v]][[j]] = exemplar_covariates[[j]]
        }
      }
      
      suppressMessages({
        newdata[[v]] = 
          predict(x,
                  newdata = newdata[[v]],
                  CI_level = CI_level,
                  PI_level = PI_level)
      })
      
      
      # Get starter plots
      two_plots = 
        (!combine_pred_cred) &
        ( ("cred band" %in% type) & ("pred band" %in% type) )
      if(two_plots){
        plot_name_v1 = 
          paste0("pred_band_",v)
        plot_name_v2 = 
          paste0("cred_band_",v)
        
        if(is.numeric(x$data[[v]])){
          plot_list[[plot_name_v1]] =
            plot_list[[plot_name_v2]] =
            x$data |>
            ggplot(aes(x = .data[[v]],
                       y = .data[[all.vars(x$formula)[1]]])) + 
            geom_point(alpha = 0.2)
        }else{
          if(x$family$family %in% c("poisson","negbinom")){
            plot_list[[plot_name_v1]] =
              plot_list[[plot_name_v2]] =
              x$data |>
              ggplot(aes(x = .data[[v]],
                         y = .data[[all.vars(x$formula)[1]]])) + 
              geom_violin(alpha = 0.2)
          }
          if(x$family$family == "binomial"){
            plot_list[[plot_name_v1]] =
              plot_list[[plot_name_v2]] =
              x$data |>
              # dplyr::group_by(get(v)) |> 
              dplyr::group_by(.data[[v]]) |> 
              dplyr::summarize(prop1 = mean(dplyr::near(.data[[all.vars(x$formula)[1]]], 1))) |> 
              # dplyr::rename(!!v := .data$`get(v)`) |> 
              ggplot(aes(x = .data[[v]],
                         y = .data$prop1)) + 
              geom_col(fill="gray70") +
              ylab(all.vars(x$formula)[1])
          }
        }
        
      }else{#End: starting two plots for bands/intervals
        
        band_to_plot = 
          paste0(gsub("\ ","_",type[grep("band",type)]),
                 "_",v)
        if(length(band_to_plot) == 2) band_to_plot = paste0("band_",v)
        
        
        for( plot_name_v in band_to_plot){
          if(is.numeric(x$data[[v]])){
            plot_list[[plot_name_v]] =
              x$data |>
              ggplot(aes(x = .data[[v]],
                         y = .data[[all.vars(x$formula)[1]]])) + 
              geom_point(alpha = 0.2)
          }else{
            if(x$family$family%in% c("poisson","negbinom")){
              plot_list[[plot_name_v]] =
                x$data |>
                ggplot(aes(x = .data[[v]],
                           y = .data[[all.vars(x$formula)[1]]])) + 
                geom_violin(alpha = 0.2)
            }
            if(x$family$family == "binomial"){
              plot_list[[plot_name_v]] =
                x$data |>
                # dplyr::group_by(get(v)) |> 
                dplyr::group_by(.data[[v]]) |> 
                dplyr::summarize(prop1 = mean(dplyr::near(.data[[all.vars(x$formula)[1]]], 1))) |> 
                # dplyr::rename(!!v := .data$`get(v)`) |> 
                ggplot(aes(x = .data[[v]],
                           y = .data$prop1)) + 
                geom_col(fill="gray70") + 
                ylab(all.vars(x$formula)[1])
            }
          }
        }
        
      }
      
      
    }#End: for loop through variables
    
  }# End: Get exemplar and PI/CI
  
  
  # Prediction Band plots
  if("pred band" %in% type){
    
    for(v in variable){
      plot_name_v = 
        paste0(ifelse((!combine_pred_cred) | !("cred band" %in% type),
                      "pred_band_","band_"),v)
      
      if(is.numeric(x_seq[[v]])){
        plot_list[[plot_name_v]] =
          plot_list[[plot_name_v]] +
          geom_ribbon(data = newdata[[v]],
                      aes(ymin = .data$PI_lower,
                          ymax = .data$PI_upper),
                      fill = "lightsteelblue3",
                      alpha = 0.5) +
          geom_line(data = newdata[[v]],
                    aes(x = .data[[v]],
                        y = .data$`Post Mean`))
      }else{
        plot_list[[plot_name_v]] =
          plot_list[[plot_name_v]] +
          geom_errorbar(data = newdata[[v]],
                        aes(x = .data[[v]],
                            ymin = .data$PI_lower,
                            ymax = .data$PI_upper),
                        color = "lightsteelblue3") +
          geom_point(data = newdata[[v]],
                     aes(x = .data[[v]],
                         y = .data$`Post Mean`),
                     size = 3)
      }
      
      
    }
    
    
    
  }
  
  if("cred band" %in% type){
    
    
    for(v in variable){
      plot_name_v = 
        paste0(ifelse((!combine_pred_cred) | !("pred band" %in% type),
                      "cred_band_","band_"),v)
      
      if(is.numeric(x_seq[[v]])){
        plot_list[[plot_name_v]] =
          plot_list[[plot_name_v]] +
          geom_ribbon(data = newdata[[v]],
                      aes(ymin = .data$CI_lower,
                          ymax = .data$CI_upper),
                      fill = "steelblue4",
                      alpha = 0.5) +
          geom_line(data = newdata[[v]],
                    aes(x = .data[[v]],
                        y = .data$`Post Mean`))
      }else{
        plot_list[[plot_name_v]] =
          plot_list[[plot_name_v]] +
          geom_errorbar(data = 
                          newdata[[v]] |> 
                          dplyr::mutate(prop1 = 0.0), # Stupid hack to make ggplot work right.
                        aes(x = .data[[v]],
                            ymin = .data$CI_lower,
                            ymax = .data$CI_upper),
                        color = "steelblue4") +
          geom_point(data = newdata[[v]],
                     aes(x = .data[[v]],
                         y = .data$`Post Mean`),
                     size = 3) +
          ylab(all.vars(x$formula)[1])
      }
    }
    
    
  }
  
  
  # Polish up plots
  if( ("pred band" %in% type) | ("cred band" %in% type) ){
    for(v in variable){
      
      for(j in names(plot_list)[grepl("band",names(plot_list)) & grepl(v,names(plot_list))]){
        plot_list[[j]] =
          plot_list[[j]] +
          theme_classic() +
          ggtitle(
            paste0(
              ifelse(
                grepl("pred_",j),
                paste0("Prediction band for ",v),
                ifelse(grepl("cred_",j),
                       paste0("Credible band for ",v),
                       paste0("Cred. and Pred. bands for ",v)
                )
              )
            )
          )
      }
      
    }
  }
  
  
  if(return_as_list){
    return(plot_list)
  }else{
    return(
      wrap_plots(plot_list)
    )
  }
  
  
}

#' @rdname plot
#' @method plot np_glm_b 
#' @export
plot.np_glm_b = function(x,
                         type,
                         variable,
                         exemplar_covariates,
                         variable_seq_length = 30,
                         return_as_list = FALSE,
                         CI_level = 0.95,
                         seed = 1,
                         backtransformation = function(x){x},
                         ...){
  
  alpha_ci = 1.0 - CI_level
  
  if(missing(type)){
    type = 
      c(#"pdp",
        "cred band")
  }
  
  type = c("pdp",
           "cred band")[pmatch(tolower(type),
                             c("pdp",
                               "cred band"))]
  
  
  if(missing(variable)){
    variable = 
      terms(x) |> 
      delete.response() |> 
      all.vars() |> 
      unique()
  }
  
  N = nrow(x$data)
  
  plot_list = list()
  
  
  
  # Get unique values and x sequences for plots
  x_unique = 
    lapply(variable,
           function(v) unique(x$data[[v]]))
  x_seq = 
    lapply(x_unique,
           function(xvals){
             if(length(xvals) > variable_seq_length){
               return( 
                 seq(min(xvals),
                     max(xvals),
                     l = variable_seq_length)
               )
             }else{
               if(is.numeric(xvals)){
                 return(sort(xvals))
               }else{
                 if(is.character(xvals)){
                   return(
                     factor(sort(xvals),
                            levels = sort(xvals))
                   )
                 }else{
                   return(xvals)
                 }
               }
             }
           })
  
  names(x_unique) = 
    names(x_seq) = variable
  
  
  # Partial Dependence Plots
  if("pdp" %in% type){
    
    for(v in 1:length(variable)){
      
      newdata = 
        tibble::tibble(var_of_interest = x_seq[[v]],
                       y = 0.0)
      suppressMessages({
        for(i in 1:length(x_seq[[v]])){
          temp_preds = 
            predict(x,
                    newdata = 
                      x$data |>
                      dplyr::mutate(!!variable[v] := newdata$var_of_interest[i]),
                    CI_level = CI_level)
          newdata$y[i] = backtransformation(mean(temp_preds$`Post Mean`))
        }
      })
      
      
      if(is.numeric(x_seq[[v]])){
        plot_list[[paste0("pdp_",variable[v])]] = 
          x$data |>
          ggplot(aes(x = .data[[variable[v]]],
                     y = .data[[all.vars(x$formula)[1]]])) + 
          geom_point(alpha = 0.2) +
          geom_line(data = newdata,
                    aes(x = .data$var_of_interest,
                        y = .data$y))
      }else{
        if(x$family$family == "binomial"){
          plot_list[[paste0("pdp_",variable[v])]] = 
            x$data |>
            dplyr::group_by(get(variable[v])) |> 
            dplyr::summarize(prop1 = mean(dplyr::near(.data[[all.vars(x$formula)[1]]], 1))) |> 
            ggplot(aes(x = .data$`get(variable[v])`,
                       y = .data$prop1)) + 
            geom_col(fill="gray70") +
            ylab(all.vars(x$formula)[1])
        }else{
          plot_list[[paste0("pdp_",variable[v])]] = 
            x$data |>
            ggplot(aes(x = .data[[variable[v]]],
                       y = .data[[all.vars(x$formula)[1]]])) + 
            geom_violin(alpha = 0.2)
        }
        
        plot_list[[paste0("pdp_",variable[v])]] = 
          plot_list[[paste0("pdp_",variable[v])]] + 
          geom_point(data = newdata,
                     aes(x = .data$var_of_interest,
                         y = .data$y),
                     size = 3)
      }
      
      plot_list[[paste0("pdp_",variable[v])]] = 
        plot_list[[paste0("pdp_",variable[v])]] + 
        xlab(variable[v]) + 
        ylab(all.vars(x$formula)[1]) + 
        theme_classic() +
        ggtitle("Partial dependence plot")
      
    }
  }# End: PDP
  
  # If drawing CI/PI bands, get reference covariate values and prediction/CIs
  if("cred band" %in% type){
    
    # Get other covariate values
    if(missing(exemplar_covariates)){
      message("Missing other covariate values in 'exemplar_covariates.'  Using medoid observation instead.")
      desmat = 
        model.matrix(x$formula,
                     x$data) |> 
        scale()
      exemplar_covariates = 
        x$data[cluster::pam(desmat,k=1)$id.med,]
    }
    
    # Get CI and PI values
    newdata = list()
    for(v in variable){
      newdata[[v]] = 
        tibble::tibble(!!v := x_seq[[v]])
      for(j in setdiff(names(exemplar_covariates),v)){
        if(is.character(exemplar_covariates[[j]])){
          newdata[[v]][[j]] = 
            factor(exemplar_covariates[[j]],
                   levels = unique(x$data[[j]]))
        }else{
          newdata[[v]][[j]] = exemplar_covariates[[j]]
        }
      }
      
      suppressMessages({
        newdata[[v]] = 
          predict(x,
                  newdata = newdata[[v]],
                  CI_level = CI_level)
      })
      newdata[[v]] = 
        newdata[[v]] |> 
        dplyr::mutate(dplyr::across(dplyr::all_of(c("Post Mean",
                                                    "CI_lower",
                                                    "CI_upper")),
                                    backtransformation))
      
      # Get starter plot
      plot_name_v = paste0("cred_band_",v)
      if(is.numeric(x$data[[v]])){
        plot_list[[plot_name_v]] =
          x$data |>
          ggplot(aes(x = .data[[v]],
                     y = .data[[all.vars(x$formula)[1]]])) + 
          geom_point(alpha = 0.2) +
          geom_ribbon(data = newdata[[v]],
                      aes(ymin = .data$CI_lower,
                          ymax = .data$CI_upper),
                      fill = "steelblue4",
                      alpha = 0.5) +
          geom_line(data = newdata[[v]],
                    aes(x = .data[[v]],
                        y = .data$`Post Mean`))
      }else{
        if(x$family$family == "binomial"){
          plot_list[[plot_name_v]] =
            x$data |>
            # dplyr::group_by(get(v)) |> 
            dplyr::group_by(.data[[v]]) |> 
            dplyr::summarize(prop1 = mean(dplyr::near(.data[[all.vars(x$formula)[1]]], 1))) |> 
            # dplyr::rename(!!v := .data$`get(v)`) |>
            ggplot(aes(x = .data[[v]],
                       y = .data$prop1)) + 
            geom_col(fill="gray70") + 
            ylab(all.vars(x$formula)[1])
        }else{
          plot_list[[plot_name_v]] =
            x$data |>
            ggplot(aes(x = .data[[v]],
                       y = .data[[all.vars(x$formula)[1]]])) + 
            geom_violin(alpha = 0.2)
        }
        
        
        plot_list[[plot_name_v]] =
          plot_list[[plot_name_v]] +
          geom_errorbar(data = 
                          newdata[[v]] |> 
                          dplyr::mutate(prop1 = 0.0), # Stupid hack to make ggplot work right.
                        aes(x = .data[[v]],
                            ymin = .data$CI_lower,
                            ymax = .data$CI_upper),
                        color = "steelblue4") +
          geom_point(data = newdata[[v]],
                     aes(x = .data[[v]],
                         y = .data$`Post Mean`),
                     size = 3) + 
          ylab(all.vars(x$formula)[1])
        
      }
      
      plot_list[[plot_name_v]] =
        plot_list[[plot_name_v]] +
        theme_classic() +
        ggtitle(paste0("Credible band for ",v))
      
    }#End: loop through variables
    
  }#End: CI band code
  
  
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
                          type,
                          return_as_list = FALSE,
                          ...){
  
  if(missing(type)){
    type = c("dx","acme","ade")
  }
  
  type = c("diagnostics",
           "diagnostics",
           "acme",
           "ade")[pmatch(tolower(type),
                         c("diagnostics",
                           "dx",
                           "acme",
                           "ade"))]
  
  plot_list = list()
  
  
  # Start diagnostic plots
  if("diagnostics" %in% type){
    # Mediator model
    ## Get dx plots
    temp = 
      plot(x$model_m,
           type = "dx",
           return_as_list = TRUE)
    ## Make titles more specific
    for(j in names(temp)){
      plot_list[[paste0(j,"_m")]] = 
        temp[[j]] + 
        ggtitle(paste0(temp[[j]]$labels$title,
                       " (Mediator model)"))
    }
    
    # Outcome model
    ## Get dx plots
    temp = 
      plot(x$model_y,
           type = "dx",
           return_as_list = TRUE)
    ## Make titles more specific
    for(j in names(temp)){
      plot_list[[paste0(j,"_y")]] = 
        temp[[j]] + 
        ggtitle(paste0(temp[[j]]$labels$title,
                       " (Outcome model)"))
    }
    
  }#End: diagnostic plots
  
  
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
