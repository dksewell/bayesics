#' b_procedure objects
#'
#' Objects of class \code{b_procedure} represent the result of a Bayesian
#' procedure, including the data, prior specification, posterior summaries, 
#' and optional plotting output.
#'
#' @details
#' A \code{b_procedure} object is a named list with the following components:
#'
#' \describe{
#'   \item{name}{Character string giving the name of the procedure.}
#'   \item{data}{A tibble containing the data used in the analysis.}
#'   \item{print_data}{Logical; whether the data should be printed by 
#'   \code{print.b_procedure}.}
#'   \item{CI_level}{Numeric scalar giving the credible interval level as 
#'   provided by the user.}
#'   \item{prior}{Character string describing the prior used.}
#'   \item{posterior_summaries}{A tibble containing posterior summaries with
#'     columns:
#'     \itemize{
#'       \item \code{Quantity}: character
#'       \item \code{Post Mean}: numeric
#'       \item \code{Lower}: numeric
#'       \item \code{Upper}: numeric
#'       \item \code{ROPE}: optional numeric
#'       \item \code{ROPE_lower_bound},\code{ROPE_upper_bound}: optional numeric
#'     }
#'   }
#'   \item{PDir}{If applicable, list containing:
#'     \itemize{
#'       \item \code{description}: character
#'       \item \code{pdir}: numeric scalar giving the probability of 
#'   direction.
#'     }
#'   }
#'   \item{BF}{If applicable, list containing:
#'     \itemize{
#'       \item \code{description}: character
#'       \item \code{BF}: numeric scalar giving the Bayes factor
#'       \item \code{interpretation}: character
#'     }
#'   }
#'   \item{PDir_description}{If applicable, character string describing the 
#'   probability of direction.}
#'   \item{plot}{A \code{ggplot} object associated with the procedure 
#'   (optional).}
#'   \item{plot_description}{Character string describing the plot (optional).}
#'   \item{object_fit}{If applicable, the underlying fitted model object 
#'   (e.g., \code{aov_b}).}
#' }
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
#' @name b_procedure
#' @docType class
NULL