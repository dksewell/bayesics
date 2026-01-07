#' Mediation using Bayesian methods
#' 
#' 
#' Mediation analysis done in the framework of Imai et al. (2010).  Currently 
#' only applicable to linear models.
#' 
#' @details
#' The model is the same as that of Imai et al. (2010):
#' \deqn{
#'  M_i(X) =  w_i'\alpha_m + X\beta_m + \epsilon_{m,i}, \\
#'  y_i(X, M(\tilde X)) = w_i'\alpha_y + X\beta_y + M(\tilde X)\gamma + \epsilon_{y,i}, \\
#'  \epsilon_{m,i} \overset{iid}{\sim} N(0,\sigma^2_m), \\
#'  \epsilon_{y,i} \overset{iid}{\sim} N(0,\sigma^2_y), \\
#' }
#' where \eqn{M_i(X)} is the mediator as a function of the treatment variable
#' \eqn{X}, and \eqn{w_i} are confounder covariates.
#' 
#' 
#' Unlike the \code{mediation} R package, the estimation in \code{mediate_b} 
#' is fully Bayesian (as opposed to "quasi-Bayesian").
#' 
#' 
#' @references 
#' 
#' Imai, Kosuke, et al. 
#' “A General Approach to Causal Mediation Analysis.” Psychological Methods, 
#' vol. 15, no. 4, 2010, pp. 309–34, https://doi.org/10.1037/a0020761.
#' 
#' @param model_m a fitted model object of class lm_b for mediator.
#' @param model_y a fitted model object of class lm_b for outcome.
#' @param treat a character string indicating the name of the 
#' treatment variable used in the models.  NOTE: Treatment variable must be
#' numeric (even if it's 1's and 0's).
#' @param control_value value of the treatment variable used as the 
#' control condition. Default is the 1st quintile of the treat variable.
#' @param treat_value value of the treatment variable used as the treatment condition. 
#' Default is the 4th quintile of the treat variable.
#' @param n_draws Number of preliminary posterior draws to assess final 
#' number of posterior draws required for accurate interval estimation
#' @param ask_before_full_sampling logical.  If FALSE, the user will not 
#' be asked if they want to complete the full sampling.  Defaults to 
#' TRUE, as this can be a computationally intensive procedure.
#' @param CI_level numeric. Credible interval level.
#' @param mc_error positive scalar.  The number of posterior samples will, 
#' with high probability, estimate the CI bounds up to 
#' \eqn{\pm}\code{mc_error}\eqn{\times}\code{sd(y)}.
#' @param batch_size positive integer.  Number of posterior draws to be 
#' taken at once.  Higher values are more computationally intensive, but 
#' values which are too high might take up significant memory (allocates 
#' on the order of \code{batch_size}\eqn{\times}\code{nrow(model_y$data)}).
#' 
#' @return A list with the following elements:
#' \itemize{
#'  \item \code{summary} - tibble giving results for causal mediation quantities
#'  \item \code{posterior_draws} (of counterfactual expectations)
#'  \item \code{mc_error} absolute error used, including any rescaling 
#'  to match the scale of the outcome
#'  \item other inputs to \code{mediate_b}
#' }
#' 
#' @export

mediate_b = function(model_m,
                     model_y,
                     treat,
                     control_value,
                     treat_value,
                     n_draws = 500,
                     ask_before_full_sampling = TRUE,
                     CI_level = 0.95,
                     seed = 1,
                     mc_error = ifelse("glm_b" %in% model_y,
                                       0.01,0.002),
                     batch_size = 500){
  set.seed(seed)
  alpha_ci = 1 - CI_level
  
  # Get mc_error
  if("lm_b" %in% class(model_y)){
    y =
      model.response(model.frame(terms(model_y),
                                 model_y$data))
    mc_error = mc_error * 4 * sd(y)
  }
  if( ("glm_b" %in% class(model_y)) &&
      (model_y$family$family != "binomial") ){
    y =
      model.response(model.frame(terms(model_y),
                                 model_y$data))
    mc_error = mc_error * 4 * sd(log(y + 1))
  }
  
  
  
  
  if(!all.equal(model_m$data,model_y$data)){
    stop("Data in model_m and model_y must match.")
  }
  
  mediator = as.character(model_m$formula)[[2]]
  
  if(missing(control_value)){
    warning(paste0("control_value missing; set to be the 1st quintile of ",
                   treat))
    control_value = quantile(model_m$data[[treat]],probs = 0.2)
  }
  if(missing(treat_value)){
    warning(paste0("treat_value missing; set to be the 4th quintile of ",
                   treat))
    treat_value = quantile(model_m$data[[treat]],probs = 0.8)
  }
  
  tl = attr(model_y$terms,"term.labels")
  simple = 
    ("lm_b" %in% class(model_m)) &
    ("lm_b" %in% class(model_y)) & 
    !any(grepl(paste0(":",mediator),tl) | 
           grepl(paste0(mediator,":"),tl))
  
  results = list()
  if(simple){
    
    # Get posterior draws for ACME and ADE
    ## Get preliminary draws
    mediator_draws = 
      get_posterior_draws(model_m,
                          n_draws = n_draws)
    outcome_draws = 
      get_posterior_draws(model_y,
                          n_draws = n_draws)
    results$posterior_draws = 
      tibble::tibble(
        ACME = 
          (treat_value - control_value) * 
          mediator_draws[,treat] *
          outcome_draws[,all.vars(model_m$formula)[1]],
        ADE =
          (treat_value - control_value) * 
          outcome_draws[,treat])
    results$posterior_draws$`Total Effect` = 
      results$posterior_draws$ACME + results$posterior_draws$ADE
    
    ## Evaluate number of draws required for accurate CI bounds
    fhats = 
      future.apply::future_lapply(1:2,
                                  function(i){
                                    density(unlist(results$posterior_draws[,i]),adjust = 2)
                                  })
    
    n_more_draws = 
      future.apply::future_sapply(1:2,
                                  function(i){
                                    0.5 * alpha_ci * (1.0 - 0.5 * alpha_ci) *
                                      (
                                        qnorm(0.5 * (1.0 - 0.99)) / 
                                          mc_error /
                                          fhats[[i]]$y[which.min(abs(fhats[[i]]$x - 
                                                                       quantile(unlist(results$posterior_draws[,i]), 0.5 * alpha_ci)))]
                                      )^2
                                  }) |> 
      max() |> 
      round() - n_draws
    
    # Finish sampling
    user_response = TRUE
    give_warning = TRUE
    
    if(n_more_draws <= 0 ){
      ask_before_full_sampling = FALSE
      user_response = FALSE
      give_warning = FALSE
    }
    
    if(ask_before_full_sampling){
      user_response = 
        utils::askYesNo(paste0(n_more_draws,
                               " more draws are required for accurate CI bounds.\nShould sampling proceed? (yes/no)"))
    }
    
    if(user_response){
      cat("Continuing on with ",
          n_more_draws,
          " more posterior samples.\n")
      
      mediator_draws = 
        get_posterior_draws(model_m,
                            n_draws = n_draws)
      outcome_draws = 
        get_posterior_draws(model_y,
                            n_draws = n_draws)
      next_draws = 
        tibble::tibble(
          ACME = 
            (treat_value - control_value) * 
            mediator_draws[,treat] *
            outcome_draws[,all.vars(model_m$formula)[1]],
          ADE =
            (treat_value - control_value) * 
            outcome_draws[,treat])
      next_draws$`Total Effect` = 
        next_draws$ACME + next_draws$ADE
      
      results$posterior_draws = 
        bind_rows(results$posterior_draws,
                  next_draws)
      
    }else{
      results$message = 
        paste0(n_draws + n_more_draws,
               " total draws are required for accurate CI bounds.")
      if(give_warning) warning(results$message)
    }
    
    results$summary = 
      tibble::tibble(Estimand = c("ACME",
                                  "ADE",
                                  "Total Effect",
                                  "Prop. Mediated"),
                     Estimate = 
                       c(mean(results$posterior_draws$ACME),
                         mean(results$posterior_draws$ADE),
                         mean(results$posterior_draws$`Total Effect`),
                         mean(results$posterior_draws$ACME / 
                                results$posterior_draws$`Total Effect`)),
                     Lower = 
                       c(quantile(results$posterior_draws$ACME,
                                  probs = 0.5 * alpha_ci),
                         quantile(results$posterior_draws$ADE,
                                  probs = 0.5 * alpha_ci),
                         quantile(results$posterior_draws$`Total Effect`,
                                  probs = 0.5 * alpha_ci),
                         quantile(results$posterior_draws$ACME / 
                                    results$posterior_draws$`Total Effect`,
                                  probs = 0.5 * alpha_ci)),
                     Upper =
                       c(quantile(results$posterior_draws$ACME,
                                  probs = 1.0 - 0.5 * alpha_ci),
                         quantile(results$posterior_draws$ADE,
                                  probs = 1.0 - 0.5 * alpha_ci),
                         quantile(results$posterior_draws$`Total Effect`,
                                  probs = 1.0 - 0.5 * alpha_ci),
                         quantile(results$posterior_draws$ACME / 
                                    results$posterior_draws$`Total Effect`,
                                  probs = 1.0 - 0.5 * alpha_ci)),
                     `Prob Dir` = 
                       c(mean(results$posterior_draws$ACME > 0),
                         mean(results$posterior_draws$ADE > 0),
                         mean(results$posterior_draws$`Total Effect` > 0),
                         NA)
      )
    
    
  }else{
    
    if(("glm_b" %in% class(model_m)) && (model_m$algorithm == "IS")){
      suppressMessages({
        model_m <- 
          glm_b(formula = model_m$formula,
                data = model_m$data,
                family = model_m$family,
                trials = model_m$trials,
                prior_beta_mean = model_m$hyperparameters$prior_beta_mean,
                prior_beta_precision = model_m$hyperparameters$prior_beta_precision,
                algorithm = "VB")
      })
    }
    if(("glm_b" %in% class(model_y)) && (model_y$algorithm == "IS")){
      suppressMessages({
        model_m <- 
          glm_b(formula = model_y$formula,
                data = model_y$data,
                family = model_y$family,
                trials = model_y$trials,
                prior_beta_mean = model_y$hyperparameters$prior_beta_mean,
                prior_beta_precision = model_y$hyperparameters$prior_beta_precision,
                algorithm = "VB")
      })
    }
    
    # Setup counterfactual data for posterior draws
    counterfactual_data0 = 
      counterfactual_data1 = 
      model_m$data
    counterfactual_data0[[treat]] = control_value
    counterfactual_data1[[treat]] = treat_value
    
    
    # Create sampling function
    draw_y_x_mx = function(n_iter){
      
      ## Draw new mediators
      M_0 =
        predict(model_m,
                newdata = counterfactual_data0,
                n_draws = n_iter)
      M_0 = 
        M_0[,setdiff(colnames(M_0),
                     c(mediator,
                       "Post Mean",
                       "PI_lower",
                       "PI_upper",
                       "CI_lower",
                       "CI_upper"))]
      M_0 = 
        M_0 |> 
        tidyr::pivot_longer(cols = contains("y_new"),
                            names_to = "posterior_draw",
                            values_to = mediator,
                            names_prefix = "y_new")
      M_1 =
        predict(model_m,
                newdata = counterfactual_data1,
                n_draws = n_iter)
      M_1 = 
        M_1[,setdiff(colnames(M_1),
                     c(mediator,
                       "Post Mean",
                       "PI_lower",
                       "PI_upper",
                       "CI_lower",
                       "CI_upper"))]
      M_1 = 
        M_1 |> 
        tidyr::pivot_longer(cols = contains("y_new"),
                            names_to = "posterior_draw",
                            values_to = mediator,
                            names_prefix = "y_new")
      gc_output = 
        utils::capture.output({gc()})
      
      
      y_00 =
        predict(model_y,
                newdata = M_0,
                n_draws = 1)
      gc_output = 
        utils::capture.output({gc()})
      y_11 =
        predict(model_y,
                newdata = M_1,
                n_draws = 1)
      gc_output = 
        utils::capture.output({gc()})
      M_1[[treat]] = control_value
      y_01 =
        predict(model_y,
                newdata = M_1,
                n_draws = 1)
      gc_output = 
        utils::capture.output({gc()})
      M_0[[treat]] = treat_value
      y_10 =
        predict(model_y,
                newdata = M_0,
                n_draws = 1)
      gc_output = 
        utils::capture.output({gc()})
      
      
      E_00 = 
        y_00 |> 
        dplyr::group_by(.data$posterior_draw) |> 
        dplyr::summarize(mean = mean(.data[["y_new1"]],)) |> 
        dplyr::pull(mean)
      E_11 = 
        y_11 |> 
        dplyr::group_by(.data$posterior_draw) |> 
        dplyr::summarize(mean = mean(.data[["y_new1"]])) |> 
        dplyr::pull(mean)
      E_01 = 
        y_01 |> 
        dplyr::group_by(.data$posterior_draw) |> 
        dplyr::summarize(mean = mean(.data[["y_new1"]])) |> 
        dplyr::pull(mean)
      E_10 = 
        y_10 |> 
        dplyr::group_by(.data$posterior_draw) |> 
        dplyr::summarize(mean = mean(.data[["y_new1"]])) |> 
        dplyr::pull(mean)
      
      ret = 
        tibble::tibble(Tot_Eff = E_11 - E_00,
                       ACME_control = E_01 - E_00,
                       ACME_treat = E_11 - E_10,
                       ADE_control = E_10 - E_00,
                       ADE_treat = E_11 - E_01)
      
      rm(M_0,M_1,y_00,y_11,y_10,y_01,E_00,E_11,E_10,E_01)
      gc_output = 
        utils::capture.output({gc()})
      
      return( ret )
    }
    
    # Get preliminary posterior draws
    prelim_draws = 
      draw_y_x_mx(n_draws) |> 
      na.omit()
    cat(paste0("\nFinished with ",
               n_draws,
               " preliminary posterior draws.\n"))
    
    
    ## Evaluate number of draws required for accurate CI bounds
    fhats = 
      future.apply::future_lapply(2:NCOL(prelim_draws),
                                  function(i){
                                    stats::density(unlist(prelim_draws[,i]),adjust = 2)
                                  })
    
    n_more_draws = 
      future.apply::future_sapply(2:NCOL(prelim_draws),
                                  function(i){
                                    0.5 * alpha_ci * (1.0 - 0.5 * alpha_ci) *
                                      (
                                        qnorm(0.5 * (1.0 - 0.99)) / 
                                          mc_error /
                                          fhats[[i - 1]]$y[which.min(abs(fhats[[i - 1]]$x - 
                                                                           quantile(unlist(prelim_draws[,i]), 0.5 * alpha_ci)))]
                                      )^2
                                  }) |> 
      max() |> 
      round() - n_draws
    
    
    user_response = TRUE
    give_warning = TRUE
    
    if(n_more_draws <= 0 ){
      ask_before_full_sampling = FALSE
      user_response = FALSE
      give_warning = FALSE
    }
    
    if(ask_before_full_sampling){
      user_response = 
        utils::askYesNo(paste0(n_more_draws,
                               " more draws are required for accurate CI bounds.\nShould sampling proceed? (yes/no)"))
    }
    
    if(user_response){
      cat("Continuing on with ",
          n_more_draws,
          " more posterior samples.\n")
      
      # Do it in batches to save memory
      batch_size_vector = 
        c(seq(1,n_more_draws,by = batch_size),n_more_draws + 1) |> 
        diff() |> 
        pmax(2)
      results$posterior_draws = 
        do.call(dplyr::bind_rows,
                future.apply::future_lapply(1:length(batch_size_vector),
                                            function(b){
                                              suppressWarnings(suppressPackageStartupMessages(library(bayesics)))
                                              draw_y_x_mx(batch_size_vector[b])
                                            },
                                            future.seed = seed + 1)
        ) |> 
        na.omit()
    }else{
      results$message = 
        paste0(n_draws + n_more_draws,
               " total draws are required for accurate CI bounds.")
      results$posterior_draws = 
        prelim_draws
      if(give_warning) warning(results$message)
    }
    
    
    
    # Put it together to return
    results$summary = 
      tibble::tibble(Estimand = c("ACME (Control)",
                                  "ACME (Treatment)",
                                  "ADE (Control)",
                                  "ADE (Treatment)",
                                  "Total Effect",
                                  "ACME (Average)",
                                  "ADE (Average)",
                                  "Prop. Mediated (Average)"),
                     Estimate = 
                       c(mean(results$posterior_draws$ACME_control),
                         mean(results$posterior_draws$ACME_treat),
                         mean(results$posterior_draws$ADE_control),
                         mean(results$posterior_draws$ADE_treat),
                         mean(results$posterior_draws$Tot_Eff),
                         0.5 * mean(results$posterior_draws$ACME_control + 
                                      results$posterior_draws$ACME_treat),
                         0.5 * mean(results$posterior_draws$ADE_control + 
                                      results$posterior_draws$ADE_treat),
                         mean( (results$posterior_draws$ACME_control + 
                                  results$posterior_draws$ACME_treat) / 
                                 (results$posterior_draws$ACME_control + 
                                    results$posterior_draws$ACME_treat + 
                                    results$posterior_draws$ADE_control + 
                                    results$posterior_draws$ADE_treat) )
                       ),
                     Lower = 
                       c(quantile(results$posterior_draws$ACME_control,0.5 * alpha_ci),
                         quantile(results$posterior_draws$ACME_treat,0.5 * alpha_ci),
                         quantile(results$posterior_draws$ADE_control,0.5 * alpha_ci),
                         quantile(results$posterior_draws$ADE_treat,0.5 * alpha_ci),
                         quantile(results$posterior_draws$Tot_Eff,0.5 * alpha_ci),
                         0.5 * quantile(results$posterior_draws$ACME_control + 
                                          results$posterior_draws$ACME_treat,0.5 * alpha_ci),
                         0.5 * quantile(results$posterior_draws$ADE_control + 
                                          results$posterior_draws$ADE_treat,0.5 * alpha_ci),
                         quantile( (results$posterior_draws$ACME_control + 
                                      results$posterior_draws$ACME_treat) / 
                                     (results$posterior_draws$ACME_control + 
                                        results$posterior_draws$ACME_treat + 
                                        results$posterior_draws$ADE_control + 
                                        results$posterior_draws$ADE_treat), 0.5 * alpha_ci )
                       ),
                     Upper = 
                       c(quantile(results$posterior_draws$ACME_control,1.0 - 0.5 * alpha_ci),
                         quantile(results$posterior_draws$ACME_treat,1.0 - 0.5 * alpha_ci),
                         quantile(results$posterior_draws$ADE_control,1.0 - 0.5 * alpha_ci),
                         quantile(results$posterior_draws$ADE_treat,1.0 - 0.5 * alpha_ci),
                         quantile(results$posterior_draws$Tot_Eff,1.0 - 0.5 * alpha_ci),
                         0.5 * quantile(results$posterior_draws$ACME_control + 
                                          results$posterior_draws$ACME_treat,1.0 - 0.5 * alpha_ci),
                         0.5 * quantile(results$posterior_draws$ADE_control + 
                                          results$posterior_draws$ADE_treat,1.0 - 0.5 * alpha_ci),
                         quantile( (results$posterior_draws$ACME_control + 
                                      results$posterior_draws$ACME_treat) / 
                                     (results$posterior_draws$ACME_control + 
                                        results$posterior_draws$ACME_treat + 
                                        results$posterior_draws$ADE_control + 
                                        results$posterior_draws$ADE_treat), 1.0 - 0.5 * alpha_ci )
                       ),
                     `Prob Dir` = 
                       c(mean(results$posterior_draws$ACME_control > 0),
                         mean(results$posterior_draws$ACME_treat > 0),
                         mean(results$posterior_draws$ADE_control > 0),
                         mean(results$posterior_draws$ADE_treat > 0),
                         mean(results$posterior_draws$Tot_Eff > 0),
                         mean(results$posterior_draws$ACME_control + 
                                results$posterior_draws$ACME_treat > 0), # No need to multiply by 0.5 for PDir
                         mean(results$posterior_draws$ADE_control + 
                                results$posterior_draws$ADE_treat > 0),
                         NA)
      )
  }
  
  results$summary$`Prob Dir` = 
    sapply(results$summary$`Prob Dir`,
           function(x) pmax(x, 1.0 - x)
    )
  
  # Don't report negative or >1 proportions
  results$summary$Estimate[nrow(results$summary)] = 
    ifelse(results$summary$Estimate[nrow(results$summary)] < 0,
           0,
           ifelse(results$summary$Estimate[nrow(results$summary)] > 1,
                  1,
                  results$summary$Estimate[nrow(results$summary)]))
  results$summary$Lower[nrow(results$summary)] = 
    min(max(results$summary$Lower[nrow(results$summary)],0.0),1.00)
  results$summary$Upper[nrow(results$summary)] = 
    max(min(results$summary$Upper[nrow(results$summary)],1.0),0.0)
  
  results$treat_value = treat_value
  results$control_value = control_value
  results$model_m = model_m
  results$model_y = model_y
  results$CI_level = CI_level
  results$mc_error = mc_error
  
  
  return(structure(results,
                   class = "mediate_b"))
}


