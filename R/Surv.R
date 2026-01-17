#' Create a Survival Object
#' 
#' Create a survival object, usually used as a response variable in a 
#' model formula. Argument matching is special for this function, see 
#' Details under \code{\link[survival]{Surv}}. This is a restricted wrapper around
#' \code{\link[survival]{Surv}} and currently supports only right-censored data.
#' 
#' @param ... arguments to be passed into \code{survival::Surv}.  Currently, 
#' the input must be of the form \code{Surv(time,event)} for right censored 
#' data.
#' @return An object of class "\code{Surv}".
#' @export
Surv = function(...){
  s <- survival::Surv(...)
  # check if right censored
  if (attr(s, "type") != "right" || ncol(s)>2) { 
    stop("Only right-censored Surv objects are supported.") # stop if not right censored or if time 2 is also supplied
  }
  return(s)
}
