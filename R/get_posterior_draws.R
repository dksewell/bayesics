#' @name get_posterior_draws
#' 
#' @title Get posterior samples from lm_b object
#' 
#' @param object Object of class lm_b
#' @param n_draws integer.  Number of posterior draws to obtain.
#' @param seed integer.
#' @param ... optional arguments.
#' 
#' @returns matrix of posterior draws
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
#' pdraws <-
#'   get_posterior_draws(fit1)
#'   }
#' 
#' 
#' @export
get_posterior_draws = function(object,...){
  UseMethod("get_posterior_draws")
}

#' @rdname credint
#' @exportS3Method get_posterior_draws lm_b 
get_posterior_draws.lm_b = function(object, 
                                    n_draws = 1e4,
                                    seed = 1){
  set.seed(seed)
  
  if("posterior_covariance" %in% names(object)){ # Handles lm, glm\IS, np_glm\bootstrapping
    
    if("sigma_sq" %in% names(object)){ # Handles lm specifically
      
      p = nrow(object$summary)
      
      # Get sqrt of unscaled covariance matrix
      V_tilde_eig = 
        eigen(object$posterior_parameters$V_tilde)
      Vinv_sqrt = 
        tcrossprod(diag(x = 1 / sqrt(V_tilde_eig$values),
                        nrow = length(V_tilde_eig$values),
                        ncol = length(V_tilde_eig$values)),
                   V_tilde_eig$vectors)
      
      post_draws = 
        matrix(0.0,n_draws,nrow(object$summary) + 1,
               dimnames = list(NULL,
                               c(object$summary$Variable,
                                 "s2")))
      post_draws[,"s2"] = 
        extraDistr::rinvgamma(n_draws,
                              0.5 * object$posterior_parameters$a_tilde,
                              0.5 * object$posterior_parameters$b_tilde)
      post_draws[,1:p] = 
        matrix(1.0,n_draws,1) %*% matrix(object$summary$`Post Mean`,nrow=1) +
        matrix(rnorm(n_draws*p,
                     sd = sqrt(rep(post_draws[,"s2"],p))),n_draws,p) %*% Vinv_sqrt
      
    }else{
      
      post_draws = 
        mvtnorm::rmvt(n_draws,
                      type = "shifted",
                      delta = object$summary$`Post Mean`,
                      sigma = object$posterior_covariance,
                      df = object$df)
      
    }#End: posterior_covariance without having to add in s2
    
  }else{
    
    if("importance_sampling_weights" %in% names(object)){ # Handles glm IS
      
      post_draws = 
        object$proposal_draws[sample(nrow(object$proposal_draws),
                                     n_draws,
                                     replace = TRUE,
                                     prob = object$importance_sampling_weights),]
      
    }else{
      if("posterior_draws" %in% names(object)){ # Handles aov, np_glm bootstrapping, bma_inference
        
        post_draws = 
          object$posterior_draws[sample(nrow(object$posterior_draws),
                                        n_draws,
                                        replace = TRUE),]
        
      }#End: posterior_draws if
    }#End: importance_sampling_weights ifelse
  }#End: posterior_covariance ifelse
  
  
  return(post_draws)
}




#' @rdname credint
#' @exportS3Method get_posterior_draws aov_b 
get_posterior_draws.aov_b = function(object,
                                     n_draws = 1e4,
                                     seed = 1){
  set.seed(seed)
  
  G = length(object$posterior_parameters$nu_g)
  
  heteroscedastic = 
    (length(object$posterior_parameters$a_g) > 1)
  
  if(heteroscedastic){
    
    
    s2_g_draws = 
      future.apply::future_sapply(1:G,
                                  function(g){
                                    extraDistr::rinvgamma(n_draws,
                                                          alpha = object$posterior_parameters$a_g[g]/2,
                                                          beta = object$posterior_parameters$b_g[g]/2)
                                  },
                                  future.seed = seed)
    mu_g_draws = 
      future.apply::future_sapply(1:G,
                                  function(g){
                                    rnorm(n_draws,
                                          mean = object$posterior_parameters$mu_g[g],
                                          sd = sqrt(s2_g_draws[,g] / object$posterior_parameters$nu_g[g]))
                                  },
                                  future.seed = seed)
    
    post_draws = 
      cbind(mu_g_draws,
            s2_g_draws)
    colnames(post_draws) =
      object$summary$Variable
    
  }else{
    
    s2_G_draws =
      extraDistr::rinvgamma(n_draws,
                            alpha = object$posterior_parameters$a_g/2,
                            beta = object$posterior_parameters$b_g/2)
    mu_g_draws = 
      future.apply::future_sapply(1:G,
                                  function(g){
                                    rnorm(n_draws,
                                          mean = object$posterior_parameters$mu_g[g],
                                          sd = sqrt(s2_G_draws / object$posterior_parameters$nu_g[g]))
                                  },
                                  future.seed = seed)
    
    post_draws = 
      cbind(mu_g_draws,
            s2_G_draws)
    colnames(post_draws) =
      object$summary$Variable
    
  }
  
  return(post_draws)
}