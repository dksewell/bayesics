#' Predict method for lm_b model fits
#' 
#' 
#' @param object Object of class \code{lm_b}, \code{glm_b}, \code{np_glm_b}, or \code{lm_b_bma}
#' @param newdata An optional data.frame in which to look for variables with which 
#' to predict. 
#' @param trials Integer vector giving the number of trials for each 
#' observation if family = binomial().
#' @param CI_level Posterior probability covered by credible interval
#' @param PI_level Posterior probability covered by prediction interval
#' @param seed integer.  Always set your seed!!!
#' @param n_draws integer.  Number of posterior draws used for prediction.  
#' Ignored if estimation method already relies on posterior sampling, in which case 
#' \code{n_draws} will match the number of posterior draws in the fitted object.
#' @param ... optional arguments.
#' 
#' @returns tibble with estimate (posterior mean), prediction intervals, and credible intervals 
#' for the mean.
#' 
#' @examples
#' \donttest{
#' 
#' # lm_b
#' ## Create data
#' N = 500
#' test_data <-
#'   data.frame(x1 = rnorm(N),
#'              x2 = rnorm(N),
#'              x3 = letters[1:5])
#' test_data$outcome <-
#'   rnorm(N,-1 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) )
#' 
#' ## Fit linear model
#' fit <-
#'   lm_b(outcome ~ x1 + x2 + x3,
#'        data = test_data)
#' predict(fit)
#' 
#' 
#' # glm_b
#' ## Generate some negative binomial data
#' set.seed(2025)
#' N = 500
#' test_data =
#'   data.frame(x1 = rnorm(N),
#'              x2 = rnorm(N),
#'              x3 = letters[1:5],
#'              time = rexp(N))
#' test_data$outcome =
#'   rnbinom(N,
#'           mu = exp(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e"))) * test_data$time,
#'           size = 0.7)
#' 
#' ## Fit using variational Bayes (default)
#' fit_vb1 <-
#'   glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
#'         data = test_data,
#'         family = negbinom(),
#'         seed = 2025)
#' # Predict
#' predict(fit_vb1)
#' 
#' ## Fit the GLM via the (non-parametric) loss-likelihood bootstrap.
#' fit_np <-
#'   np_glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
#'            data = test_data,
#'            family = negbinom())
#' predict(fit_np)
#' 
#' 
#' # bma_inference
#' ## Create data
#' set.seed(2025)
#' N = 500
#' test_data = 
#'   data.frame(x1 = rnorm(N),
#'              x2 = rnorm(N),
#'              x3 = letters[1:5],
#'              x4 = rnorm(N),
#'              x5 = rnorm(N),
#'              x6 = rnorm(N),
#'              x7 = rnorm(N),
#'              x8 = rnorm(N),
#'              x9 = rnorm(N),
#'              x10 = rnorm(N))
#' test_data$outcome = 
#'   rnorm(N,-1 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) )
#' 
#' ## Fit linear model using Bayesian model averaging
#' fit <-
#'   bma_inference(outcome ~ .,
#'                 test_data,
#'                 user.int = FALSE)
#' predict(fit)
#' 
#' 
#' }
#' 

#' @rdname predict
#' @exportS3Method predict lm_b
predict.lm_b = function(object,
                        newdata,
                        trials,
                        CI_level = 0.95,
                        PI_level = 0.95,
                        seed = 1,
                        n_draws = 5e3,
                        ...){
  
  alpha_ci = 1.0 - CI_level
  alpha_pi = 1.0 - PI_level
  
  
  if(missing(newdata)){
    newdata = object$data
  }
  
  if(!is.null(object$xlevels)){
    for(j in names(object$xlevels)){
      if(!("factor" %in% class(newdata[[j]]))){
        newdata[[j]] = 
          factor(newdata[[j]],
                 levels = object$xlevels[[j]])
      }
    }
  }
  
  # Extract 
  mframe = model.frame(delete.response(terms(object)),
                       data = newdata)
  X = model.matrix(delete.response(terms(object)),
                   data = newdata)
  os = model.offset(mframe)
  N = nrow(X)
  p = ncol(X)
  alpha = 1 - CI_level
  if(is.null(os)) os = numeric(N)
  
  # Get trials variables sorted
  if(object$family$family == "binomial"){
    if(missing(trials)){
      message("Assuming all observations correspond to Bernoulli, i.e., Binomial with one trial.")
      trials = rep(1.0,N)
    }else{
      if(is.character(trials)) trials = newdata[[trials]]
      trials = as.numeric(trials)
    }
  }else{
    trials = rep(1.0,N)
  }
  
  
  # Assign new values
  ## Get estimates
  yhats = 
    trials * 
    object$family$linkinv(eta = 
                            drop(X %*% 
                                   object$summary$`Post Mean`[1:p]) + 
                            os)
  
  
  ## Get CI bounds
  if("posterior_covariance" %in% names(object)){ # Handles lm, glm\IS, np_glm\bootstrapping
    
    if(object$family$family == "gaussian"){
      grad_ginv_xbeta = 
        X
    }
    if(object$family$family %in% c("poisson","negbinom")){
      grad_ginv_xbeta = 
        drop(exp(X %*% object$summary$`Post Mean`[1:p] + os)) * X
    }
    if(object$family$family == "binomial"){
      probs = 
        1.0 / (1.0 + drop(exp(-X %*% object$summary$`Post Mean` - os)))
      grad_ginv_xbeta = 
        trials * probs * (1.0 - probs) * X
    }
    
    yhats_covar =
      rowSums(
        grad_ginv_xbeta * 
          (grad_ginv_xbeta %*%
             object$posterior_covariance[1:ncol(X),1:ncol(X)]
          )
      )
    # The above equals
    # yhats_covar = 
    #       tcrossprod(grad_ginv_xbeta %*% object$posterior_covariance[1:p,1:p],
    #                  grad_ginv_xbeta)
    yhats_sds = 
      sqrt(yhats_covar)
    
    
    newdata =
      newdata |> 
      dplyr::mutate(`Post Mean` = yhats,
                    CI_lower = 
                      qlst(alpha / 2.0,
                           object$df,
                           yhats,
                           yhats_sds),
                    CI_upper = 
                      qlst(1.0 - alpha / 2.0,
                           object$df,
                           yhats,
                           yhats_sds))
    
    
    if(object$model_type == "parametric"){
      
      if(object$family$family == "gaussian"){
        
        newdata =
          newdata |> 
          dplyr::mutate(PI_lower = 
                          extraDistr::qlst(alpha_pi / 2.0,
                                           df = object$posterior_parameters$a_tilde,
                                           mu = .data$`Post Mean`,
                                           sigma = sqrt(yhats_sds^2 +
                                                          object$posterior_parameters$b_tilde / 
                                                          object$posterior_parameters$a_tilde) ),
                        PI_upper = 
                          extraDistr::qlst(1.0 - alpha_pi / 2.0,
                                           df = object$posterior_parameters$a_tilde,
                                           mu = .data$`Post Mean`,
                                           sigma = sqrt(yhats_sds^2 +
                                                          object$posterior_parameters$b_tilde / 
                                                          object$posterior_parameters$a_tilde) )
          )
        
        if(n_draws > 0){
          y_draws = matrix(0.0,
                             nrow(newdata),
                             n_draws,
                             dimnames = list(NULL,
                                             paste("y_new",1:n_draws,sep="")))
          for(it in 1:n_draws){
            y_draws[,it] = 
              extraDistr::rlst(nrow(newdata),
                               df = object$posterior_parameters$a_tilde,
                               mu = newdata$`Post Mean`,
                               sigma = sqrt(yhats_sds^2 +
                                              object$posterior_parameters$b_tilde / 
                                              object$posterior_parameters$a_tilde) )
          }
          
          colnames(y_draws) = paste("y_new",1:ncol(y_draws),sep="")
          newdata =
            newdata |> 
            dplyr::bind_cols(y_draws |> 
                               as_tibble())
        }
        
      }else{
      
        if(object$family$family == "poisson"){
          y_draws =
            future.apply::future_sapply(1:n_draws, 
                                        function(i){
                                          rpois(nrow(newdata),
                                                pmax(rnorm(nrow(newdata),
                                                           newdata$`Post Mean`,
                                                           yhats_sds),
                                                     .Machine$double.eps)
                                          )
                                        },
                                        future.seed = seed)
          if( ((n_draws > 1) && (NCOL(y_draws) == 1)) || (length(y_draws) == 1) )
            y_draws = matrix(y_draws,nrow = 1)
          
        }
        if(object$family$family == "binomial"){
          y_draws =
            future.apply::future_sapply(1:n_draws, 
                                        function(i){
                                          rbinom(nrow(newdata),
                                                 trials,
                                                 pmin(
                                                   pmax(
                                                     rnorm(nrow(newdata),
                                                           newdata$`Post Mean`,
                                                           yhats_sds),
                                                     .Machine$double.eps),
                                                   1.0 - .Machine$double.eps)
                                          )
                                        },
                                        future.seed = seed)
          if( ((n_draws > 1) && (NCOL(y_draws) == 1)) || (length(y_draws) == 1) )
            y_draws = matrix(y_draws,nrow = 1)
        }
        if(object$family$family == "negbinom"){
          theta_draws =
            mvtnorm::rmvnorm(n_draws,
                             mean = object$summary$`Post Mean`,
                             sigma = object$posterior_covariance)
          y_draws =
            future.apply::future_sapply(1:n_draws, 
                                        function(i){
                                          rnbinom(nrow(newdata),
                                                  mu = pmax(drop(exp(X %*% theta_draws[i,1:ncol(X)] + os)),
                                                            .Machine$double.eps),
                                                  size = exp(theta_draws[i,ncol(X) + 1])
                                          )
                                        },
                                        future.seed = seed)
          if( ((n_draws > 1) && (NCOL(y_draws) == 1)) || (length(y_draws) == 1) )
            y_draws = matrix(y_draws,nrow = 1)
        }
        
        if(n_draws > 1){
          PI_bounds = 
            y_draws |> 
            future.apply::future_apply(1,quantile,probs = c(0.5 * alpha_pi,
                                                            1.0 - 0.5 * alpha_pi))
          newdata$PI_lower = 
            PI_bounds[1,]
          newdata$PI_upper = 
            PI_bounds[2,]
        }
        
        colnames(y_draws) = paste("y_new",1:ncol(y_draws),sep="")
        newdata =
          newdata |> 
          dplyr::bind_cols(y_draws |> 
                             as_tibble())
        
      }#End: Prediction for non-gaussian models
      
    }#End: Prediction (for parametric models)
    
  }else{#End: if asymptotic approx or VB was used.
    
    if("importance_sampling_weights" %in% names(object)){ # Handles glm IS
    
      yhat_draws = 
        trials * 
        object$family$linkinv(os + tcrossprod(X, as.matrix(object$proposal_draws)[,1:ncol(X)]))
      
      CI_from_weighted_sample = function(x,w,level = alpha_ci){
        w = cumsum(w[order(x)])
        x = x[order(x)]
        LB = max(which(w <= 0.5 * level))
        UB = min(which(w >= 1.0 - 0.5 * level))
        return(c(lower = x[LB],
                 upper = x[UB]))
      }
      CI_bounds = 
        apply(yhat_draws,1,
              CI_from_weighted_sample,
              w = object$importance_sampling_weights)
      
      newdata =
        newdata |> 
        tibble::as_tibble() |> 
        dplyr::mutate(`Post Mean` = yhats,
                      CI_lower = 
                        CI_bounds["lower",],
                      CI_upper = 
                        CI_bounds["upper",])
      
      if(object$family$family == "poisson"){
        y_draws = 
          future.apply::future_sapply(1:ncol(yhat_draws), 
                                      function(i){
                                        rpois(nrow(yhat_draws),yhat_draws[,i])
                                      },
                                      future.seed = seed)
        if(NCOL(y_draws) == 1)
          y_draws = matrix(y_draws,nrow = 1)
      }
      
      if(object$family$family == "binomial"){
        y_draws = 
          future.apply::future_sapply(1:ncol(yhat_draws),
                                      function(i){
                                        rbinom(nrow(yhat_draws),
                                               trials,
                                               yhat_draws[,i] / trials)
                                      },
                                      future.seed = seed)
        if(NCOL(y_draws) == 1)
          y_draws = matrix(y_draws,nrow = 1)
      }
      
      if(object$family$family == "negbinom"){
        y_draws =
          future.apply::future_sapply(1:ncol(yhat_draws),
                                      function(i){
                                        rnbinom(nrow(newdata),
                                                mu = pmax(drop(exp(X %*% object$proposal_draws[i,1:ncol(X)] + os)),
                                                          .Machine$double.eps),
                                                size = exp(object$proposal_draws[i,ncol(X) + 1])
                                        )
                                      },
                                      future.seed = seed)
        if(NCOL(y_draws) == 1)
          y_draws = matrix(y_draws,nrow = 1)
        
      }
      
      PI_bounds = NULL
      try({
        PI_bounds = 
          y_draws |> 
          future.apply::future_apply(1,CI_from_weighted_sample,
                                     w = object$importance_sampling_weights,
                                     level = alpha_pi)
      },silent=TRUE)
      if(is.null(PI_bounds)){
        PI_bounds = 
          y_draws |> 
          apply(1,CI_from_weighted_sample,
                w = object$importance_sampling_weights,
                level = alpha_pi)
      }
      newdata$PI_lower = 
        PI_bounds[1,]
      newdata$PI_upper = 
        PI_bounds[2,]
      
      if(n_draws > 0){
        y_sir_draws = NULL
        try({
          y_sir_draws = 
            future.apply::future_apply(y_draws,1,
                                       function(x){
                                         x[sample(ncol(y_draws),
                                                  n_draws,
                                                  replace = TRUE,
                                                  prob = object$importance_sampling_weights)]
                                       },
                                       future.seed = seed) |> 
            t()
        },silent = TRUE)
        if(is.null(y_sir_draws)){
          y_sir_draws = 
            apply(y_draws,1,
                  function(x){
                    x[sample(ncol(y_draws),
                             n_draws,
                             replace = TRUE,
                             prob = object$importance_sampling_weights)]
                  }) |> 
            t()
        }
        
        colnames(y_sir_draws) = paste("y_new",1:n_draws,sep="")
        newdata =
          newdata |> 
          dplyr::bind_cols(y_sir_draws |> 
                             as_tibble())
        
      }
      
    }else{#End: IS
      if("posterior_draws" %in% names(object)){ # Handles np_glm bootstrapping, bma_inference
        # Get means of E(y|X)
        yhats = 
          trials * 
          object$family$linkinv(eta = 
                                  drop(X %*% 
                                         object$summary$`Post Mean`[1:p]) + 
                                  os)
        
        # Get draws of y_new
        yhat_draws = 
          trials * 
          object$family$linkinv(os + tcrossprod(X, as.matrix(object$posterior_draws)[,1:p]))
        
        if(NCOL(yhat_draws) == 1)
          yhat_draws = matrix(yhat_draws,nrow = 1)
        
        newdata =
          newdata |> 
          tibble::as_tibble() |> 
          dplyr::mutate(`Post Mean` = yhats,
                        CI_lower = 
                          yhat_draws |> 
                          apply(1,quantile, probs = alpha / 2.0),
                        CI_upper = 
                          yhat_draws |> 
                          apply(1,quantile, probs = 1.0 - alpha / 2.0))
        
        
        # Get prediction intervals
        if(object$model_type == "parametric"){
          
          y_draws = 
            yhat_draws + 
            sweep(matrix(rnorm(prod(dim(yhat_draws))),
                         nrow(yhat_draws),
                         ncol(yhat_draws)),
                  2,
                  sqrt(object$posterior_draws$s2),
                  "*")
            
          
          newdata =
            newdata |> 
            dplyr::mutate(PI_lower = 
                            apply(y_draws,1,quantile,probs = 0.5 * alpha_pi),
                          PI_upper = 
                            apply(y_draws,1,quantile,probs = 1.0 - 0.5 * alpha_pi))
          
          colnames(y_draws) = paste("y_new",1:ncol(y_draws),sep="")
          newdata =
            newdata |> 
            dplyr::bind_cols(y_draws |> 
                               as_tibble())
          
        }#End: PI intervals for bma_inference objects (gaussian only)
        
      }#End: np_glm bootstrapping, bma_inference
    }#End: IS, np_glm bootstrapping, bma_inference
    
  }#End: all other estimation algos
  
  
  # Correct for delta method bound errors
  if(object$family$family == "binomial"){
    newdata = 
      newdata |>
      dplyr::mutate(across(c(CI_lower,
                             CI_upper),
                           ~ ifelse(.x < 0, 0,
                                    ifelse(.x > 1,
                                           1,
                                           .x))))
    if("PI_lower" %in% names(newdata)){
      newdata = 
        newdata |>
        dplyr::mutate(across(c(PI_lower,
                               PI_upper),
                             ~ ifelse(.x < 0, 0,
                                      ifelse(.x > 1,
                                             1,
                                             .x))))
    }
    
  }
  if(object$family$family %in% c("poisson","negbinom")){
    newdata = 
      newdata |>
      dplyr::mutate(across(c(CI_lower,
                             CI_upper),
                           ~ ifelse(.x < 0, 0,.x)))
    if("PI_lower" %in% names(newdata)){
      newdata = 
        newdata |>
        dplyr::mutate(across(c(PI_lower,
                               PI_upper),
                             ~ ifelse(.x < 0, 0,.x)))
    }
  }
  
  
  return(newdata)
}




#' @rdname predict
#' @exportS3Method predict aov_b
predict.aov_b = function(object,
                         CI_level = 0.95,
                         PI_level = 0.95,
                        ...){
  
  alpha_ci = 1.0 - CI_level
  alpha_pi = 1.0 - PI_level
  
  G = length(object$posterior_parameters$nu_g)
  
  newdata = 
    tibble(group = 
             levels(object$data$group),
           `Post Mean` = object$summary$`Post Mean`[1:G],
           CI_lower = 
             extraDistr::qlst(alpha_ci / 2.0, 
                              df = object$posterior_parameters$a_g,
                              mu = object$posterior_parameters$mu_g,
                              sigma = sqrt(object$posterior_parameters$b_g / 
                                             object$posterior_parameters$a_g / 
                                             object$posterior_parameters$nu_g)),
           CI_upper = 
             extraDistr::qlst(1.0 - alpha_ci / 2.0, 
                              df = object$posterior_parameters$a_g,
                              mu = object$posterior_parameters$mu_g,
                              sigma = sqrt(object$posterior_parameters$b_g / 
                                             object$posterior_parameters$a_g / 
                                             object$posterior_parameters$nu_g)),
           PI_lower = 
             extraDistr::qlst(alpha_ci / 2.0, 
                              df = object$posterior_parameters$a_g,
                              mu = object$posterior_parameters$mu_g,
                              sigma = sqrt(object$posterior_parameters$b_g / 
                                             object$posterior_parameters$a_g * 
                                             (1.0 + 1.0 / object$posterior_parameters$nu_g))),
           PI_upper = 
             extraDistr::qlst(1.0 - alpha_ci / 2.0, 
                              df = object$posterior_parameters$a_g,
                              mu = object$posterior_parameters$mu_g,
                              sigma = sqrt(object$posterior_parameters$b_g / 
                                             object$posterior_parameters$a_g * 
                                             (1.0 + 1.0 / object$posterior_parameters$nu_g)))
    )
  
  
}

