#' @rdname plot_bands
#' 
#' @title Plot Credible and Prediction Bands
#' 
#' @param x object of class \code{aov_b}, \code{lm_b}, or \code{glm_b}
#' @param type character. Select "cred band", and/or "pred band".  
#'  NOTE: the credible and prediction bands only work for numeric 
#'  variables.
#' @param combine_pred_cred logical. If type includes both "cred band" and "pred band", 
#' should the credible band be superimposed on the prediction band or 
#' plotted separately?
#' @param CI_level Posterior probability covered by credible interval
#' @param PI_level Posterior probability covered by prediction interval
#' @param backtransformation function.  If a transformation of 
#' the response variable was used, \code{backtransformation} 
#' should be the inverse of this transformation function.  E.g., 
#' if you fit lm_b(log(y) ~ x), then set \code{backtransformation=exp}. 
#' @param return_as_list logical.  If TRUE, a list of ggplots will be returned, 
#' rather than a single plot produced by the patchwork package.
#' @param variable character. If type = "pdp" , which variable should be plotted?
#' @param variable_seq_length integer. Number of points used to draw pdp.
#' @param exemplar_covariates data.frame or tibble with exactly one row.  
#' Used to fix other covariates while varying the variable of interest for the plot.
#' @param ... arguments passed on to plot_bands
#' 
#' @export
plot_bands = function(x,
                      type,
                      combine_pred_cred,
                      CI_level,
                      PI_level,
                      backtransformation,
                      return_as_list,
                      ...){
  UseMethod("plot_bands")
}


#' @rdname plot_bands
#' @exportS3Method plot_bands lm_b
plot_bands.lm_b = function(x,
                           type = c("cred band",
                                    "pred band"),
                           combine_pred_cred = TRUE,
                           CI_level = 0.95,
                           PI_level = 0.95,
                           backtransformation = function(x){x},
                           return_as_list = TRUE,
                           variable,
                           variable_seq_length = 30,
                           exemplar_covariates,
                           ...){
  
  type = match.arg(type,
                   several.ok = TRUE)
  
  plot_list = list()
  
  if( (x$family$family == "binomial") & 
      ("pred band" %in% type) ){
    warning("Prediction band cannot be supplied for a binomial outcome.\nResults shown will be credible band instead.")
    type = "cred band"
  }
  if( (x$model_type == "nonparametric") & 
      ("pred band" %in% type) ){
    warning("Prediction band cannot be supplied for a non-parametric model.\nResults shown will be credible band instead.")
    type = "cred band"
  }
  
  
  
  if(is.character(x$data[[all.vars(x$formula)[1]]])){
    x$data[[all.vars(x$formula)[1]]] = 
      factor(x$data[[all.vars(x$formula)[1]]])
  }
  if(is.factor(x$data[[all.vars(x$formula)[1]]])){
    x$data[[all.vars(x$formula)[1]]] = as.integer(x$data[[all.vars(x$formula)[1]]])
    if(length(unique(x$data[[all.vars(x$formula)[1]]])) == 2) x$data[[all.vars(x$formula)[1]]] = x$data[[all.vars(x$formula)[1]]] - 1
  }
  
  # to01 = function(x) {
  #   if(is.factor(x)){
  #     as.numeric(x) - 1.0  # maps level 1 -> 0, level 2 -> 1
  #   }else{
  #     x                   # leave numeric (or logical) as-is
  #   }
  # }
  # x$data[[all.vars(x$formula)[1]]] = 
  #   to01(x$data[[all.vars(x$formula)[1]]])
  
  
  if(missing(variable)){
    variable = 
      terms(x) |> 
      delete.response() |> 
      all.vars() |> 
      unique()
    if(!is.null(attributes(terms(x))$offset)){
      variable = 
        variable[-(attributes(terms(x))$offset - 1)]
    }
  }
  
  # Get unique values and x sequences for plots
  x_unique = 
    lapply(variable,
           function(v) na.omit(unique(x$data[[v]])))
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
                                                  "CI_lower",
                                                  "CI_upper")),
                                  backtransformation)) # below causes no visible binding for global variable ‘Post Mean’ note 
    # dplyr::mutate(dplyr::across(`Post Mean`:CI_upper,backtransformation))
    if("PI_lower" %in% names(newdata[[v]])){
      newdata[[v]] = 
        newdata[[v]] |> 
        dplyr::mutate(dplyr::across(dplyr::all_of(c("PI_lower",
                                                    "PI_upper")),
                                    backtransformation))
    }
    
    
    
    if(is.character(newdata[[v]][[all.vars(x$formula)[1]]])){
      newdata[[v]][[all.vars(x$formula)[1]]] = 
        factor(newdata[[v]][[all.vars(x$formula)[1]]])
    }
    if(is.factor(newdata[[v]][[all.vars(x$formula)[1]]])){
      newdata[[v]][[all.vars(x$formula)[1]]] = 
        as.integer(newdata[[v]][[all.vars(x$formula)[1]]])
      if(length(unique(newdata[[v]][[all.vars(x$formula)[1]]])) == 2) newdata[[v]][[all.vars(x$formula)[1]]] = newdata[[v]][[all.vars(x$formula)[1]]] - 1
    }
    
  }
  
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


#' @rdname plot_bands
#' @exportS3Method plot_bands aov_b
plot_bands.aov_b = function(x,
                            type = c("cred band",
                                     "pred band"),
                            combine_pred_cred = TRUE,
                            CI_level = 0.95,
                            PI_level = 0.95,
                            backtransformation = function(x){x},
                            return_as_list = TRUE,
                            ...){
  
  type = match.arg(type,
                   several.ok = TRUE)
  
  plot_list = list()
  
  # Get CI and PI values
  newdata =
    predict(x,
            CI_level = CI_level,
            PI_level = PI_level)
  newdata = 
    newdata |> 
    dplyr::mutate(dplyr::across(dplyr::all_of(c("Post Mean",
                                                "CI_lower",
                                                "CI_upper",
                                                "PI_lower",
                                                "PI_upper")),
                                backtransformation)) # below causes no visible binding for global variable ‘Post Mean’ note 
  
  
  
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
                    aes(x = .data$group,
                        y = .data$`Post Mean`,
                        ymin = .data$PI_lower,
                        ymax = .data$PI_upper),
                    color = "lightsteelblue3") +
      geom_point(data = newdata,
                 aes(x = .data$group,
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
                    aes(x = .data$group,
                        y = .data$`Post Mean`,
                        ymin = .data$CI_lower,
                        ymax = .data$CI_upper),
                    color = "steelblue4") +
      geom_point(data = newdata,
                 aes(x = .data$group,
                     y = .data$`Post Mean`),
                 size = 3)
    
  }
  
  # Polish up plots
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
  
  
  if(return_as_list){
    return(plot_list)
  }else{
    return(
      wrap_plots(plot_list)
    )
  }
}


