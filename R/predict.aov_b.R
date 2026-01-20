#' Predict method for aov_b model fits
#' 
#' 
#' @param object Object of class aov_b
#' @param CI_level Posterior probability covered by credible interval
#' @param PI_level Posterior probability covered by prediction interval
#' @param ... optional arguments.
#'  
#' @returns tibble with estimate (posterior mean), prediction intervals, and credible intervals 
#' for the mean.
#' 
#' @examples
#' \donttest{
#' set.seed(2025)
#' N = 500
#' test_data = 
#'   data.frame(x1 = rep(letters[1:5],N/5))
#' test_data$outcome = 
#'   rnorm(N,-1 + 2 * (test_data$x1 %in% c("d","e")) )
#' 
#' # Fit 1-way ANOVA model
#' fit1 <-
#'   aov_b(outcome ~ x1,
#'         test_data,
#'         prior_mean_mu = 2,
#'         prior_mean_nu = 0.5,
#'         prior_var_shape = 0.01,
#'         prior_var_rate = 0.01)
#' predict(fit1)
#' }
#' 
#' 
#' 
#' @exportS3Method predict aov_b


predict.aov_b = function(object, 
                         CI_level = 0.95, 
                         PI_level = 0.95,
                         ...){
  
  alpha_ci = 1.0 - CI_level
  alpha_pi = 1.0 - PI_level
  G = length(levels(object$data$group))
  
  return(
    object$summary |>
      dplyr::filter(dplyr::row_number() <= G) |> 
      dplyr::select(dplyr::all_of(c("Variable", "Post Mean"))) |> 
      dplyr::mutate(Variable = gsub("Mean : ",
                                    "",
                                    gsub(paste0(all.vars(object$formula)[2],
                                                " : "),
                                         "",
                                         .data$Variable))) |> 
      dplyr::rename(!!all.vars(object$formula)[2] := dplyr::all_of("Variable")) |> 
      dplyr::mutate(PI_lower = 
                      extraDistr::qlst(alpha_pi/2.0,
                                       df = object$posterior_parameters$a_g,
                                       mu = .data$`Post Mean`,
                                       sigma = sqrt(object$posterior_parameters$b_g / 
                                                      object$posterior_parameters$a_g * 
                                                      (1.0/object$posterior_parameters$nu_g + 1.0) ) ),
                    PI_upper = 
                      extraDistr::qlst(1.0 - alpha_pi/2.0,
                                       df = object$posterior_parameters$a_g,
                                       mu = .data$`Post Mean`,
                                       sigma = sqrt(object$posterior_parameters$b_g / 
                                                      object$posterior_parameters$a_g * 
                                                      (1.0/object$posterior_parameters$nu_g + 1.0) ) ),
                    CI_lower = 
                      extraDistr::qlst(alpha_ci/2.0,
                                       df = object$posterior_parameters$a_g,
                                       mu = .data$`Post Mean`,
                                       sigma = sqrt(object$posterior_parameters$b_g / 
                                                      object$posterior_parameters$a_g /
                                                      object$posterior_parameters$nu_g ) ),
                    CI_upper = 
                      extraDistr::qlst(1.0 - alpha_ci/2.0,
                                       df = object$posterior_parameters$a_g,
                                       mu = .data$`Post Mean`,
                                       sigma = sqrt(object$posterior_parameters$b_g / 
                                                      object$posterior_parameters$a_g /
                                                      object$posterior_parameters$nu_g ) )
      )
  )
}




