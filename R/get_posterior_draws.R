#' Get posterior samples from lm_b object
#' 
#' @param object Object of class lm_b
#' @param n_draws integer.  Number of posterior draws to obtain.
#' @param seed integer.  Always set your seed!!!
#' 
#' @returns matrix of posterior draws
#' 
#' @examples
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
#' 
#' 
#' @export

get_posterior_draws = function(object, 
                               n_draws = 1e4,
                               seed = 1){
  set.seed(seed)
  
  y = model.frame(object$formula,
                  object$data)[,1]
  X = model.matrix(object$formula,
                   object$data)
  
  p = nrow(object$summary)
  
  V_tilde_eig = eigen(object$posterior_parameters$V_tilde)
  Vinv_sqrt = tcrossprod(diag(x = 1 / sqrt(V_tilde_eig$values),
                              nrow = length(V_tilde_eig$values),
                              ncol = length(V_tilde_eig$values)),
                         V_tilde_eig$vectors)
  
  post_draws = 
    matrix(0.0,n_draws,p + 1,
           dimnames = list(NULL,
                           c(object$summary$Variable,"s2")))
  post_draws[,"s2"] = 
    extraDistr::rinvgamma(n_draws,
                          0.5 * object$posterior_parameters$a_tilde,
                          0.5 * object$posterior_parameters$b_tilde)
  post_draws[,1:p] = 
    matrix(1.0,n_draws,1) %*% matrix(object$summary$`Post Mean`,nrow=1) +
    matrix(rnorm(n_draws*p,
                 sd = sqrt(rep(post_draws[,"s2"],p))),n_draws,p) %*% Vinv_sqrt
  
  
  return(post_draws)
}