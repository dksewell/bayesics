#' Fractional Bayes factors
#' 
#' Compute fractional Bayes factors for lm_b objects
#' 
#' @details
#' Fractional Bayes factors, devised by O'Hagan, are a way to use flat, 
#' even improper, priors to obtain valid Bayes factors.  The idea is built 
#' on the notion of partial Bayes factors, where a part of the data is used 
#' to determine the prior, and the remaining is used to compare the models.
#' 
#' @param object1 object of class \code{lm_b}
#' @param object2 object of class \code{lm_b}
#' @param fractional_proportion The fraction of the data used to create the
#' prior in turn used to compute the marginal likelihood.  By default, 
#' O'Hagan's recommendation of \eqn{max(log(n),ncol(X) + 1) / n)} is used.
#' 
#' 
#' @references 
#' O’Hagan, Anthony. “Fractional Bayes Factors for Model Comparison.” Journal of the Royal Statistical Society. Series B (Methodological), vol. 57, no. 1, 1995, pp. 99–138. https://doi.org/10.1111/j.2517-6161.1995.tb02017.x
#' 
#' 
#' @examples
#' \donttest{
#' set.seed(2026)
#' N = 500
#' test_data <-
#'   data.frame(x1 = rnorm(N),
#'              x2 = rnorm(N),
#'              x3 = letters[1:5])
#' test_data$outcome <-
#'   rnorm(N,-1 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) )
#' fit_full <-
#'   lm_b(outcome ~ x1 + x2 + x3,
#'        data = test_data)
#' fit_no_x1 <-
#'   lm_b(outcome ~ x2 + x3,
#'        data = test_data)
#' fit_no_x2 <-
#'   lm_b(outcome ~ x1 + x3,
#'        data = test_data)
#' 
#' frac_bayes_factors(fit_full,
#'                    fit_no_x1)
#' frac_bayes_factors(fit_full,
#'                    fit_no_x2)
#' 
#' }
#' 
#' 
#' @importFrom stats lm
#' @import Matrix
#' @export



frac_bayes_factors = function(object1,
                              object2,
                              fractional_proportion){
  
  if(!all.equal(object1$data,
                object2$data))
    stop("Data must be the same in both objects.")
  
  if(missing(fractional_proportion)){
    message('By default, the fraction of data "used" is max(ncol(X) + 1,log(n)) / n.')
    fractional_proportion = 
      max(nrow(object1$posterior_parameters$V_tilde) + 1, # Minimal dataset for proper posterior
          nrow(object2$posterior_parameters$V_tilde) + 1,
          log(nrow(object1$data))) / 
      nrow(object1$data)
  }
  
  # Get key quantities from data.
  m1 = model.frame(terms(object1),
                   data = object1$data)
  y = model.response(m1)
  X1 = model.matrix(delete.response(terms(object1)),
                    data = object1$data)
  
  m2 = model.frame(terms(object2),
                   data = object2$data)
  X2 = model.matrix(delete.response(terms(object2)),
                    data = object2$data)
  check_y = model.response(m2)
  if(!isTRUE(all.equal(y,check_y)))
    stop("Response variable must be the same in both models.")
  
  N = nrow(X1)
  r1 = ncol(X1)
  r2 = ncol(X2)
  lm1 = 
    lm(object1$formula,
       object1$data)
  SSE1 = sum(lm1$residuals^2)
  lm2 = 
    lm(object2$formula,
       object2$data)
  SSE2 = sum(lm2$residuals^2)
  
  
  # Compute Bayes factor
  BF = 
    exp( 
      lgamma(0.5 * (N - r1)) +
        lgamma(0.5 * (N * fractional_proportion - r2)) - 
        lgamma(0.5 * (N - r2)) -
        lgamma(0.5 * (N * fractional_proportion - r1)) -
        0.5 * N *(1.0 - fractional_proportion) * (log(SSE1) - log(SSE2))
    )
  
  bf_max = 
    pmax(BF, 1.0 / BF)
  
  Interpretation = 
    ifelse(bf_max <= 3.2,
           "Not worth more than a bare mention",
           ifelse(bf_max <= 10,
                  "Substantial",
                  ifelse(bf_max <= 100,
                         "Strong",
                         "Decisive"))) |> 
    paste(ifelse(BF > 1,
                 "(in favor of the first model)",
                 "(against the first model)"))
  
  message("\n----------\n")
  message(paste0("The fractional Bayes factor equaled ",
                 format(signif(BF, 3)),
                 ".\nInterpretation: ",
                 Interpretation))
  message("\n----------\n\n")
  
  invisible(tibble::tibble(`fractional BF` = BF,
                           Interpretation = Interpretation))
}
