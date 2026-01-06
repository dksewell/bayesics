#' Find parameters for Beta prior based on prior mean and one quantile
#' 
#' @param mean numeric between 0 and 1 giving the prior mean
#' @param quantile numeric between 0 and 1 giving the quantile lying at left_tail_prob
#' @param left_tail_prob numeric between 0 and 1 giving the prior probability of theta 
#' being less than or equal to quantile
#' @param search_bounds bounds with which to search.  Sometimes you need to adjust this to get a good solution.
#' @param plot_results logical. Should the resulting inverse gamma distribution be plotted?
#' @examples
#' find_beta_parms(2/5,0.68,0.9)
#' 2/ (2 + 3)
#' qbeta(0.9,2,3)
#'
#'
#' @export 


find_beta_parms = function(mean,quantile,left_tail_prob, plot_results = TRUE, search_bounds = c(1e-3,100)){
  if(missing(mean)) stop("Must provide mean, quantile, and the CDF at that quantile")
  helper = function(a){
    b = a * (1.0 - mean) / mean
    (stats::pbeta(quantile,a,b) - left_tail_prob)^2
  }
  opt = stats::optimize(f = helper,
                 interval = search_bounds)
  a = opt$min
  b = a * (1.0 - mean) / mean
  
  if(plot_results){
    graphics::curve(dbeta(x,a,b),
          from= 0, 
          to = 1,
          lwd = 3,
          yaxt = "n",
          xlab = expression(theta))
  }
  
  
  return(c(shape1 = a, shape2 = b))
}
