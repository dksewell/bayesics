#' Negative-binomial family
#' 
#' The \code{negbinom()} is an additional family to be considered 
#' alongside others documented under \code{stats::family}.
#' 
#' @export
negbinom = function(){
  list(family = "negbinom",
       link = "log",
       linkfun = 
         function(mu){
           log(mu)
         },
       linkinv = function(eta){
         pmax(exp(eta), .Machine$double.eps)
       },
       variance = function(mu,phi){
         mu + mu^2 / phi
       },
       aic = function(y,n,mu,wt,dev){
         -2.0 * sum(dnbinom(y,mu = mu,size=dev,log=TRUE) * wt)
       }
  ) |> 
    structure(class = "family")
}

