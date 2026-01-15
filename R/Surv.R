#' Create a Survival Object
#' 
#' Create a survival object, usually used as a response variable in a 
#' model formula. Argument matching is special for this function, see 
#' Details under \code{survival::Surv}, as this is purely a wrapper for it.
#' 
#' @param ... arguments to be passed into \code{survival::Surv}.  Currently, 
#' the input must be of the form \code{Surv(time,event)} for right censored 
#' data.
#' 
#' @export
Surv = function(...){
  survival::Surv(...)
}
