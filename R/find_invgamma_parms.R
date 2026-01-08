#' Find parameters for Inverse gamma prior based on prior mean and one quantile
#' 
#' @details
#' Either provide the lower and upper quantiles that contain \code{probability} 
#' of the inverse gamma distribution, or if this is for linear regression, you 
#' can specify that you are a priori \code{probability} sure that the coefficient 
#' of determination (\eqn{R^2}) falls within the two bounds provided, assuming 
#' that the residual variance is \eqn{1-R^2} times the total variance.
#' 
#' @param lower_quantile lower quantile desired
#' @param upper_quantile upper quantile desired
#' @param response_variance variance of the response variable of the regression model
#' @param lower_R2,upper_R2 We are a priori \code{probability} sure that the 
#' coefficient of determination (\eqn{R^2}) falls within these lower and upper bounds.
#' @param probability prior probability to be contained within the lower and upper quantiles
#' @param plot_results logical.  Should the resulting inverse gamma distribution be plotted?
#' 
#' @returns twice the shape and rate of the inverse gamma distribution.
#' 
#' @examples
#' # When aimed at linear regression via coefficient of determination...
#' hypothetical_s2_y = 2.0
#' lower_R2 = 0.05
#' upper_R2 = 0.85
#' find_invgamma_parms(response_variance = hypothetical_s2_y,
#'                     lower_R2 = lower_R2,
#'                     upper_R2 = upper_R2,
#'                     probability = 0.8)
#' 
#' # More arbitrary task...
#' find_invgamma_parms(0.3, # hypothetical_s2_y * (1.0 - upper_R2)
#'                     1.9, #hypothetical_s2_y * (1.0 - lower_R2)
#'                     probability = 0.8)
#' 
#' 
#' 
#' 
#' @export 


find_invgamma_parms = function(lower_quantile,
                               upper_quantile,
                               response_variance,
                               lower_R2,
                               upper_R2,
                               probability,
                               plot_results = TRUE){
  tail_prob = 0.5 * (1.0 - probability)
  
  if(missing(lower_quantile) | 
     missing(upper_quantile)){
    
    if(missing(response_variance) |
       missing(lower_R2) |
       missing(upper_R2))
      stop("Must provide EITHER\nboth lower_quantile and upper_quantile, \nOR\nprovide response_variance, lower_R2, and upper_R2.")
    
    if( (lower_R2 < 0) | (lower_R2 >= 1) )
      stop("lower_R2 must be between 0 and 1")
    
    if( (upper_R2 <= 0) | (upper_R2 > 1) )
      stop("upper_R2 must be between 0 and 1")
    
    if(lower_R2 >= upper_R2)
      stop("lower_R2 must be less than upper_R2")
    
    
    lower_quantile = (1.0 - upper_R2) * response_variance
    upper_quantile = (1.0 - lower_R2) * response_variance
  }
  
  if(missing(probability))
    stop("Must provide the probability to be contained within the two quantiles")
  
  helper = function(x){
    a = exp(x[1])
    b = exp(x[2])
    (extraDistr::pinvgamma(upper_quantile,a,b) - (1.0 - tail_prob))^2 + 
      (extraDistr::pinvgamma(lower_quantile,a,b) - tail_prob)^2
  }
  opt = stats::optim(numeric(2),
                     fn = helper)
  a = exp(opt$par[1])
  b = exp(opt$par[2])
  
  if(plot_results){
    temp_dinvgamma <- function(x) extraDistr::dinvgamma(x,a,b)
    graphics::curve(temp_dinvgamma,
                    from= extraDistr::qinvgamma(0.01,a,b),
                    to = extraDistr::qinvgamma(0.99,a,b),
                    lwd = 3,
                    yaxt = "n",
                    ylab = "",
                    xlab = expression(sigma^2))
    if(!missing(response_variance))
      graphics::abline(v = response_variance)
  }
  
  
  return(c(shape = a, rate = b) * 2)
}

