#' @rawNamespace import(Matrix, except = image)
#' @import ggplot2
#' @importFrom dplyr rename group_by summarize mutate left_join n relocate near pull bind_rows bind_cols across filter select row_number all_of contains everything where
#' @importFrom extraDistr rinvgamma pinvgamma dinvgamma qinvgamma qlst plst dlst rlst rdirichlet pbbinom
#' @importFrom future.apply future_sapply future_lapply
#' @importFrom future plan multisession sequential
#' @importFrom tibble tibble as_tibble
#' @importFrom mvtnorm dmvt rmvt dmvnorm rmvnorm
#' @importFrom janitor clean_names
#' @importFrom patchwork wrap_plots
#' @importFrom cluster pam
#' @importFrom DFBA dfba_bivariate_concordance dfba_wilcoxon dfba_mann_whitney
#' @importFrom graphics curve
#' @importFrom utils capture.output askYesNo combn
#' @importFrom tidyr pivot_longer
#' @importFrom rlang .data :=
#' @importFrom survival Surv
#' @importFrom stats AIC BIC as.formula binomial coef cov dbeta dbinom delete.response density dgamma dlnorm dnbinom dnorm dpois gaussian glm lm model.frame model.matrix model.offset model.response na.omit optim optimize pbeta pgamma pnorm poisson predict qbeta qgamma qlnorm qnorm quantile rbeta rbinom resid rgamma rnbinom rnorm rpois sd sigma terms var vcov weighted.mean
# #' @import coda
NULL

