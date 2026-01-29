#' @rdname heteroscedasticity_test
#' @aliases var_test_b
#' 
#' @title 
#' Test for heteroscedasticity in AOV models
#' 
#' @description
#' Use Chib's method to compute the Bayes factor to test for 
#' heteroscedasticity in analysis of variance models.
#' 
#' @param hetero_model aov_b object where the heteroscedastic argument 
#' has been set to TRUE (the default)
#' @param homo_model aov_b object where the heteroscedastic argument 
#' has been set to FALSE
#' 
#' @returns (returned invisible) A tibble with Bayes factors and interpretations.
#' 
#' @examples
#' \donttest{
#' # Test homoscedastic case
#' ## Generate some data
#' set.seed(2025)
#' N = 200
#' test_data = 
#'   data.frame(x1 = rep(letters[1:5],N/5))
#' test_data$outcome = 
#'   rnorm(N,-1 + 2 * (test_data$x1 %in% c("d","e")) )
#' 
#' ## Fit the anova models
#' hetero_model = 
#'   aov_b(outcome ~ x1,
#'         test_data)
#' homo_model = 
#'   aov_b(outcome ~ x1,
#'         test_data,
#'         heteroscedastic = FALSE)
#' 
#' ## Perform test for heteroscedasticity using Bayes factors
#' heteroscedasticity_test(hetero_model,
#'                         homo_model)
#' 
#' # Test heteroscedastic case
#' ## Generate some data
#' set.seed(2025)
#' N = 200
#' test_data = 
#'   data.frame(x1 = rep(letters[1:5],N/5))
#' test_data$outcome = 
#'   rnorm(N,
#'         -1 + 2 * (test_data$x1 %in% c("d","e")),
#'         sd = 3 - 2 * (test_data$x1 %in% c("d","e")))
#' 
#' ## Fit the anova models
#' hetero_model = 
#'   aov_b(outcome ~ x1,
#'         test_data)
#' homo_model = 
#'   aov_b(outcome ~ x1,
#'         test_data,
#'         heteroscedastic = FALSE)
#' 
#' ## Perform test for heteroscedasticity using Bayes factors
#' heteroscedasticity_test(hetero_model,
#'                         homo_model)
#' }
#' 
#' 
#' @references 
#' 
#' Kass, R. E., & Raftery, A. E. (1995). Bayes Factors. Journal of the American Statistical Association, 90(430), 773–795.
#' 
#' @export


heteroscedasticity_test = function(hetero_model,
                                   homo_model){
  
  # Checks
  if( !inherits(hetero_model,"aov_b") |
      !inherits(homo_model,"aov_b") ) stop("Models must be of class aov_b.")
  if(homo_model$summary$Variable[nrow(homo_model$summary)] != "Var") 
    stop("homo_model must be homoscedastic.")
  if(sum(grepl("Var : ",
               hetero_model$summary$Variable)) <= 1) 
    stop("homo_model must be homoscedastic.")
  if(!isTRUE(
    all.equal(homo_model$data[[all.vars(homo_model$formula)[1]]],
              hetero_model$data[[all.vars(hetero_model$formula)[1]]])))
    stop("Outcomes must be the same in both models")
  if( isTRUE(is.na(homo_model$hyperparameters)) | 
      isTRUE(is.na(hetero_model$hyperparameters)) )
    stop("Both models must have a proper prior.")
  
  # Extract
  y = homo_model$data[[all.vars(homo_model$formula)[1]]]
  
  marg_lik = list()
  
  # Homoscedastic model
  s2_hat = homo_model$summary$`Post Mean`[nrow(homo_model$summary)]
  marg_lik$homo = 
    ## Likelihood
    sum(
      dnorm(y,
            homo_model$fitted,
            sqrt(s2_hat),
            log = TRUE)
    ) +
    ## Prior
    sum(
      dnorm(homo_model$posterior_parameters$mu_g,
            homo_model$hyperparameters$mu,
            sqrt(s2_hat / homo_model$hyperparameters$nu),
            log = TRUE)
    ) + 
    extraDistr::dinvgamma(s2_hat,
                          0.5 * homo_model$hyperparameters$a,
                          0.5 * homo_model$hyperparameters$b,
                          log = TRUE) - 
    ## Posterior
    sum(
      dnorm(homo_model$posterior_parameters$mu_g,
            homo_model$posterior_parameters$mu_g,
            sqrt(s2_hat / homo_model$posterior_parameters$nu_g),
            log = TRUE)
    ) - 
    extraDistr::dinvgamma(s2_hat,
                          0.5 * homo_model$posterior_parameters$a_g,
                          0.5 * homo_model$posterior_parameters$b_g,
                          log = TRUE)
  
  # Heteroscedastic model
  s2g_hat = 
    tibble::tibble(s2 = 
                     hetero_model$summary$`Post Mean`[grep("Var : ",hetero_model$summary$Variable)],
                   group = 
                     hetero_model$summary$Variable[grep("Var : ",hetero_model$summary$Variable)] |>
                     strsplit(":") |> 
                     sapply(function(x) x[length(x)]) |> 
                     trimws()
    )
  s2i = 
    dplyr::left_join(
      hetero_model$data,
      s2g_hat,
      by = "group")
  marg_lik$hetero = 
    ## Likelihood
    sum(
      dnorm(y,
            hetero_model$fitted,
            sqrt(s2i$s2),
            log = TRUE)
    ) +
    ## Prior
    sum(
      dnorm(hetero_model$posterior_parameters$mu_g,
            hetero_model$hyperparameters$mu,
            sqrt(s2g_hat$s2 / hetero_model$hyperparameters$nu),
            log = TRUE)
    ) + 
    sum(
      extraDistr::dinvgamma(s2g_hat$s2,
                            0.5 * hetero_model$hyperparameters$a,
                            0.5 * hetero_model$hyperparameters$b,
                            log = TRUE)
    ) - 
    ## Posterior
    sum(
      dnorm(hetero_model$posterior_parameters$mu_g,
            hetero_model$posterior_parameters$mu_g,
            sqrt(s2g_hat$s2 / hetero_model$posterior_parameters$nu_g),
            log = TRUE)
    ) - 
    sum(
      extraDistr::dinvgamma(s2g_hat$s2,
                            0.5 * hetero_model$posterior_parameters$a_g,
                            0.5 * hetero_model$posterior_parameters$b_g,
                            log = TRUE)
    )
  
  
  log_BF = marg_lik$homo - marg_lik$hetero
  BF = exp(log_BF)
  bf_max = max(BF, 1.0 / BF)
  Interpretation = 
    ifelse(bf_max <= 3.2,
           "Not worth more than a bare mention",
           ifelse(bf_max <= 10,
                  "Substantial",
                  ifelse(bf_max <= 100,
                         "Strong",
                         "Decisive")))
  Interpretation = 
    paste0(Interpretation,
           ifelse(BF > 1,
                  " in favor of homoscedasticity",
                  " in favor of heteroscedasticity")
    )
  
  
  message("\n----------\n\nTest for heteroscedasticity in 1-way ANOVA models.\n")
  message("\n\n")
  message(paste0("Bayes factor in favor of homoscedasticity = ",
             BF,
             "\nLevel of evidence: ",
             Interpretation))
  message("\n\n----------\n\n")
  
  
  invisible(tibble::tibble(BF = BF,
                           Interpretation = Interpretation))
}
