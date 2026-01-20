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
#' 
#' @returns An object of class "\code{Surv}".
#' 
#' @examples
#' \donttest{
#' set.seed(2025)
#' N = 300
#' test_data = 
#'   data.frame(outcome = 
#'                rweibull(N,2,5))
#' test_data$observed = 
#'   ifelse(test_data$outcome >= 7, 0, 1)
#' test_data$outcome =
#'   ifelse(dplyr::near(test_data$observed,1), test_data$outcome, 7)
#' Surv(test_data$outcome,
#'      test_data$observed)
#' }
#' 
#' @export
Surv = function(...){
  s <- survival::Surv(...)
  # check if right censored
  if (attr(s, "type") != "right" || ncol(s)>2) { 
    stop("Only right-censored Surv objects are supported.") # stop if not right censored or if time 2 is also supplied
  }
  return(s)
}
