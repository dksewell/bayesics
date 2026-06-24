#' \code{lm_b} Objects
#'
#' Objects of Class \code{lm_b} 
#'
#' @details
#' An object of class \code{lm_b} contains at least the following:
#'
#' \describe{
#'   \item{summary}{Tibble giving the summary of the model parameters of 
#'   having a minimum of:
#'      \describe{
#'        \item \code{Variable}: character
#'        \item \code{Post Mean}: Numeric
#'        \item \code{Lower}: Numeric
#'        \item \code{Upper}: Numeric
#'        \item \code{Prob Dir}: Numeric
#'      }
#'      
#'   }
#'   \item{formula}{}
#'   \item{data}{A tibble containing the data used in the analysis.}
#'   \item{CI_level}{Numeric scalar giving the credible interval level as 
#'   provided by the user.}
#'   \item{fitted}{Vector of fitted values}
#'   \item{residuals}{Vector of Pearson residuals}
#'   \item{family}{}
#'   \item{xlevels}{Named list, giving the levels for each factor covariate.}
#'   \item{model_type}{ Character, either "parametric" or "nonparametric".}
#' }
#' 
#' Objects from \code{\link{lm_b}} have the following additional elements:
#' \describe{
#'    \item{}{}
#'    \item{}{}
#'    \item{}{}
#'    \item{}{}
#'    \item{}{}
#'    \item{}{}
#'    \item{}{}
#'    \item{}{}
#'    \item{}{}
#'    \item{}{}
#'    \item{}{}
#' }
#' 
#' 
#' @section S3 methods:
#' Methods are available for \code{print()} and \code{plot()},
#' depending on which components are present.
#'
#' @seealso
#' \code{\link{print.b_procedure}},
#' \code{\link{plot.b_procedure}}
#'
#' @examples
#' \dontrun{
#' cc_fit <- case_control_b(matrix(c(8,47,1,26),2,2))
#' cc_fit
#' plot(cc_fit)
#' }
#'
#' @name lm_b-class
#' @docType class
NULL