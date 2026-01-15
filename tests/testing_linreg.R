# Change mc_error to scale s_y.  Be consistent on defaults across functions
# Changes need to be made to lm_b/glm_b/np_glm_b/aov_b to handle splines, 
#   as well as their plot and predict generics
# Augment SDRatio so that sets of columns of X (e.g., factor variables, 
#    splines,...) can be evaluated at once.
# Make sure np_glm_b works with vector inputs, and with intercept only models.
# homogeneity_b
# independence_b
# get_posterior_draws.lm_b
# get_posterior_draws.glm_b
# mediate_b for glm_b objects?
# Add example code to documentation
# Add SUBSET
# Maybe add AR(p) 


run <- FALSE

if (run) {
  
  library(future)
  
  
  # Helper functions --------------------------------------------------------
  
  beta_parms = find_beta_parms(2/5,0.68,0.9)
  beta_parms[1] / sum(beta_parms)
  pbeta(0.68,beta_parms[1],beta_parms[2])
  
  hypothetical_s2_y = 2.0
  lower_R2 = 0.05
  upper_R2 = 0.85
  igamma_parms = find_invgamma_parms(response_variance = hypothetical_s2_y,
                                     lower_R2 = lower_R2,
                                     upper_R2 = upper_R2,
                                     probability = 0.8)
  extraDistr::pinvgamma(hypothetical_s2_y * (1.0 - upper_R2),
                        igamma_parms[1] * 0.5,
                        igamma_parms[2] * 0.5)
  extraDistr::pinvgamma(hypothetical_s2_y * (1.0 - lower_R2),
                        igamma_parms[1] * 0.5,
                        igamma_parms[2] * 0.5,
                        lower.tail = FALSE)
  
  igamma_parms
  find_invgamma_parms(hypothetical_s2_y * (1.0 - upper_R2),
                      hypothetical_s2_y * (1.0 - lower_R2),
                      probability = 0.8)
  
  
  rm(list=ls())
  
  
  
  # Linear regression -------------------------------------------------------
  
  # Create data
  
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5])
  test_data$outcome = 
    rnorm(N,-1 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) )
  
  
  # Conjugate prior
  ## Make sure CI_level works (and print.lm_b works)
  fita = 
    lm_b(outcome ~ x1 + x2 + x3,
         data = test_data,
         prior = "conj")
  fitb = 
    lm_b(outcome ~ x1 + x2 + x3,
         data = test_data,
         prior = "conj",
         CI_level = 0.9)
  fita
  fitb
  
  ## Test input format
  fitc = 
    lm_b(test_data$outcome ~ test_data$x1 + test_data$x2 + test_data$x3,
         prior = "conj")
  fitc
  
  ## Test number of inputs
  fitd = 
    lm_b(test_data$outcome ~ test_data$x1,
         prior = "conj")
  fitd
  fite = 
    lm_b(test_data$outcome ~ 1,
         prior = "conj")
  fite
  fitf = 
    lm_b(outcome ~ x1,
         data = test_data,
         prior = "conj")
  fitf
  fitg = 
    lm_b(outcome ~ 1,
         data = test_data,
         prior = "conj")
  fitg
  
  
  
  ## Make sure summary.lm_b works
  summary(fita)
  
  ## Make sure coef.lm_b works
  coef(fita)
  
  ## Make sure credint works
  credint(fita)
  
  ## Make sure vcov works
  vcov(fita)
  
  ## Make sure prior hyperparameters work
  fith = 
    lm_b(outcome ~ x1 + x2 + x3,
         data = test_data,
         prior = "conj",
         prior_beta_mean = rep(2,7),
         prior_beta_precision = diag(7) * 0.1,
         prior_var_shape = 2,
         prior_var_rate = 2,
         CI_leve = 0.9)
  fita
  fith
  
  ## Make sure prediction function works
  head(predict(fita))
  
  ## Make sure information criteria work
  null_model = 
    lm_b(outcome ~ 1,
         data = test_data)
  AIC(fita)
  AIC(null_model)
  BIC(fita)
  BIC(null_model)
  DIC(fita)
  DIC(null_model)
  WAIC(fita)
  WAIC(null_model)
  
  ## Make sure Savage-Dickey ratio works
  SDratio(fita)
  
  ## Make sure plotting function works
  plot(fita,
       type = "diagnostics")
  plot(fita,
       type = "pdp")
  test = 
    plot(fita,
         type = "pdp",
         return_as_list = TRUE)
  patchwork::wrap_plots(test)
  plot(fita,
       type = "pdp",
       variable = "x1")
  plot(fita,
       type = c("ci","pi"),
       variable = "x1",
       combine_pi_ci = TRUE,
       exemplar_covariates = fita$data[1,])
  plot(fita,
       type = c("ci","pi"),
       combine_pi_ci = TRUE,
       exemplar_covariates = fita$data[1,])
  plot(fita,
       type = c("ci","pi"),
       combine_pi_ci = TRUE)
  plot(fita,
       type = "pi",
       variable = "x1")
  plot(fita,
       type = "pi",
       variable = "x1",
       PI_level = 0.5)
  plot(fita,
       type = "pi")
  plot(fita,
       type = "pi",
       exemplar_covariates = fita$data[1,])
  plot(fita,
       type = "ci",
       variable = "x1")
  plot(fita,
       type = "ci")
  plot(fita,
       type = "ci",
       CI_level = 0.999)
  plot(fita)
  
  # Not all generics will work if data was not supplied.  Not sure if this is 
  # worth fixing.
  plot(fitc,
       type = "diagnostics")
  plot(fitc,
       type = "pdp")
  plot(fitc,
       type = "pdp",
       variable = "x1")
  plot(fitc,
       type = c("ci","pi"),
       variable = "x1",
       combine_pi_ci = TRUE,
       exemplar_covariates = fitc$data[1,])
  plot(fitc,
       type = c("ci","pi"),
       combine_pi_ci = TRUE,
       exemplar_covariates = fitc$data[1,])
  plot(fitc,
       type = c("ci","pi"),
       combine_pi_ci = TRUE)
  plot(fitc,
       type = "pi",
       variable = "x1")
  plot(fitc,
       type = "pi",
       variable = "x1",
       PI_level = 0.5)
  plot(fitc,
       type = "pi")
  plot(fitc,
       type = "pi",
       exemplar_covariates = fitc$data[1,])
  plot(fitc,
       type = "ci",
       variable = "x1")
  plot(fitc,
       type = "ci")
  plot(fitc,
       type = "ci",
       CI_level = 0.999)
  plot(fitc)
  
  
  
  # Zellner-g prior
  ## Make sure CI_level works (and print.lm_b works)
  fita = 
    lm_b(outcome ~ x1 + x2 + x3,
         data = test_data,
         prior = "zellner")
  fitb = 
    lm_b(outcome ~ x1 + x2 + x3,
         data = test_data,
         prior = "zellner",
         CI_level = 0.9)
  fita
  fitb
  
  ## Make sure summary.lm_b works
  summary(fita)
  
  ## Make sure coef.lm_b works
  coef(fita)
  
  ## Make sure credint works
  credint(fita)
  
  ## Make sure vcov works
  vcov(fita)
  
  ## Make sure prior hyperparameters work
  fitc = 
    lm_b(outcome ~ x1 + x2 + x3,
         data = test_data,
         prior = "zellner",
         zellner_g = 10)
  fita
  fitc
  
  
  ## Test input format
  fitc = 
    lm_b(test_data$outcome ~ test_data$x1 + test_data$x2 + test_data$x3,
         prior = "conj")
  fitc
  
  ## Test number of inputs
  fitd = 
    lm_b(test_data$outcome ~ test_data$x1,
         prior = "zellner")
  fitd
  fite = 
    lm_b(test_data$outcome ~ 1,
         prior = "zellner")
  fite
  fitf = 
    lm_b(outcome ~ x1,
         data = test_data,
         prior = "zellner")
  fitf
  fitg = 
    lm_b(outcome ~ 1,
         data = test_data,
         prior = "zellner")
  fitg
  
  
  ## Make sure prediction function works
  head(predict(fita))
  
  
  ## Make sure information criteria work
  AIC(fita)
  BIC(fita)
  DIC(fita)
  WAIC(fita)
  
  ## Make sure Savage-Dickey ratio works
  SDratio(fita)
  
  ## Make sure plotting function works
  plot(fita,
       type = "diagnostics")
  plot(fita,
       type = "pdp")
  test = 
    plot(fita,
         type = "pdp",
         return_as_list = TRUE)
  patchwork::wrap_plots(test)
  plot(fita,
       type = "pdp",
       variable = "x1")
  plot(fita,
       type = c("ci","pi"),
       variable = "x1",
       combine_pi_ci = TRUE,
       exemplar_covariates = fita$data[1,])
  plot(fita,
       type = c("ci","pi"),
       combine_pi_ci = TRUE,
       exemplar_covariates = fita$data[1,])
  plot(fita,
       type = c("ci","pi"),
       combine_pi_ci = TRUE)
  plot(fita,
       type = "pi",
       variable = "x1")
  plot(fita,
       type = "pi")
  plot(fita,
       type = "pi",
       exemplar_covariates = fita$data[1,])
  plot(fita,
       type = "ci",
       variable = "x1")
  plot(fita,
       type = "ci")
  
  
  rm(list = setdiff(ls(),"test_data"))
  
  
  # Improper prior
  ## Make sure CI_level works (and print.lm_b works)
  fita = 
    lm_b(outcome ~ x1 + x2 + x3,
         data = test_data,
         prior = "impr")
  fitb = 
    lm_b(outcome ~ x1 + x2 + x3,
         data = test_data,
         prior = "impr",
         CI_level = 0.9)
  fita
  fitb
  
  ## Make sure summary.lm_b works
  summary(fita)
  
  ## Make sure coef.lm_b works
  coef(fita)
  
  ## Make sure credint works
  credint(fita)
  
  ## Make sure vcov works
  vcov(fita)
  
  ## Test number of inputs
  fitd = 
    lm_b(test_data$outcome ~ test_data$x1,
         prior = "impr")
  fitd
  fite = 
    lm_b(test_data$outcome ~ 1,
         prior = "impr")
  fite
  fitf = 
    lm_b(outcome ~ x1,
         data = test_data,
         prior = "impr")
  fitf
  fitg = 
    lm_b(outcome ~ 1,
         data = test_data,
         prior = "impr")
  fitg
  
  
  ## Make sure prior hyperparameters DO NOT work
  fitc = 
    lm_b(outcome ~ x1 + x2 + x3,
         data = test_data,
         prior = "impr",
         zellner_g = 5,
         prior_beta_mean = rep(2,7),
         prior_beta_precision = diag(7) * 0.1,
         prior_var_shape = 2,
         prior_var_rate = 2,
         CI_leve = 0.95)
  all.equal(summary(fita),
            summary(fitc))
  
  ## Make sure prediction function works
  head(predict(fita))
  predict(fita,
          newdata = fita$data[1,])
  
  ## Make sure information criteria work
  AIC(fita)
  BIC(fita)
  DIC(fita)
  WAIC(fita)
  
  ## Make sure Savage-Dickey ratio DOES NOT works
  SDratio(fita)
  
  ## Make sure plotting function works
  plot(fita,
       type = "diagnostics")
  plot(fita,
       type = "pdp")
  test = 
    plot(fita,
         type = "pdp",
         return_as_list = TRUE)
  patchwork::wrap_plots(test)
  plot(fita,
       type = "pdp",
       variable = "x1")
  plot(fita,
       type = c("ci","pi"),
       variable = "x1",
       combine_pi_ci = TRUE,
       exemplar_covariates = fita$data[1,])
  plot(fita,
       type = c("ci","pi"),
       combine_pi_ci = TRUE,
       exemplar_covariates = fita$data[1,])
  plot(fita,
       type = c("ci","pi"),
       combine_pi_ci = TRUE)
  plot(fita,
       type = "pi",
       variable = "x1")
  plot(fita,
       type = "pi")
  plot(fita,
       type = "pi",
       exemplar_covariates = fita$data[1,])
  plot(fita,
       type = "ci",
       variable = "x1")
  plot(fita,
       type = "ci")
  
  
  # Make sure WLS works
  
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5],
               w = rbeta(N,2,2))
  test_data$outcome = 
    rnorm(N,
          -1 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")),
          1.0 / sqrt(test_data$w))
  
  
  # Conjugate prior
  ## Make sure CI_level works (and print.lm_b works)
  fita = 
    lm_b(outcome ~ x1 + x2 + x3,
         data = test_data,
         prior = "conj")
  fitb = 
    lm_b(outcome ~ x1 + x2 + x3,
         data = test_data,
         weights = test_data$w,
         prior = "conj",
         CI_level = 0.9)
  summary(fita)
  summary(fitb)
  
  
  
  
  # Check if quadratic terms work
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5])
  test_data$outcome = 
    rnorm(N,-1 + test_data$x1 + test_data$x1^2 + 2 * (test_data$x3 %in% c("d","e")) )
  fita = 
    lm_b(outcome ~ x1 + I(x1^2) + x2 + x3,
         data = test_data)
  fita
  summary(fita)
  plot(fita)
  
  
  # Check if interaction terms work
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5])
  test_data$outcome = 
    rnorm(N,-1 + test_data$x1 + test_data$x1 - 2 * test_data$x1 * test_data$x2 + 2 * (test_data$x3 %in% c("d","e")) )
  fita = 
    lm_b(outcome ~ x1 * x2 + x3,
         data = test_data)
  fita
  summary(fita)
  plot(fita)
  
  
  
  # Check if splines terms work
  library(splines)
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5])
  test_data$outcome = 
    rnorm(N,-1 + test_data$x1 + test_data$x1^2 + 2 * (test_data$x3 %in% c("d","e")) )
  fita = 
    lm_b(outcome ~ ns(x1,df = 5) + x2,
         data = test_data)
  fita
  summary(fita)
  predict(fita,
          newdata = 
            test_data[1,])
  plot(fita)
  
  # Check if factors work
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5])
  test_data$outcome = 
    rnorm(N,-1 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) )
  fitb = 
    lm_b(outcome ~ x1 + x2 + x3,
         data = test_data)
  fitb
  summary(fitb)
  predict(fitb,
          newdata = 
            test_data[1,])
  plot(fitb)
  
  
  rm(list=ls())
  
  
  # Check SDratio for complicated terms
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5])
  test_data$outcome = 
    rnorm(N,-1 + 
            test_data$x1 + 
            0.25 * I(test_data$x1^2) - 
            #0.5 * test_data$x1 * test_data$x2 + 
            2 * (test_data$x3 %in% c("d","e")) )
  SDratio(lm_b(outcome ~ x1 * x2 + x3 + I(x1^2),
               data = test_data),
          by = "coef")
  SDratio(lm_b(outcome ~ x1 * x2 + x3 + I(x1^2),
               data = test_data),
          by = "vari")
  library(splines)
  SDratio(lm_b(outcome ~ ns(x1,df = 5) + x2,
               data = test_data),
          by = "coef")
  SDratio(lm_b(outcome ~ ns(x1,df = 5) + x2,
               data = test_data),
          by = "var")
  
  rm(list=ls())
  
  # Check BMA ---------------------------------------------------------------
  
  
  # Create data
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5],
               x4 = rnorm(N),
               x5 = rnorm(N),
               x6 = rnorm(N),
               x7 = rnorm(N),
               x8 = rnorm(N),
               x9 = rnorm(N),
               x10 = rnorm(N))
  test_data$outcome = 
    rnorm(N,-1 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) )
  
  # Check to make sure print.lm_b_bma and coef.lm_b_bma works
  fita = 
    bma_inference(outcome ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10,
                  test_data,
                  user.int = FALSE)
  fita
  coef(fita)
  summary(fita)
  credint(fita)
  
  # Make sure parallelization works. 
  plan(multisession, workers = 10)
  fitb = 
    bma_inference(outcome ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10,
                  test_data,
                  mc_draws = 5e3,
                  user.int = FALSE)
  plan(sequential)
  fitb
  
  ## Make sure summary.lm_b_bma works
  summary(fita)
  summary(fita,
          CI_level = 0.8)
  
  
  # Test predictions
  preds = 
    predict(fita)
  colnames(preds$newdata)
  plot(outcome ~ `Post Mean`,
       data = preds$newdata)
  dim(preds$posterior_draws$mean_of_ynew)
  dim(preds$posterior_draws$ynew)
  
  
  # Test plotting
  plot(fita,
       "dx")
  
  plot(fita,
       "pdp",
       variable = c("x1","x2","x3"))
  
  plot(fita,
       "ci",
       variable = c("x1","x2","x3"))
  plot(fita,
       "ci",
       variable = c("x1","x2","x3"),
       CI_level = 0.99)
  plot(fita,
       "pi",
       variable = c("x1","x2","x3"))
  plot(fita,
       "pi",
       variable = c("x1","x2","x3"),
       PI_level = 0.5)
  plot(fita,
       c("ci","pi"),
       variable = c("x1","x2","x3"))
  plot(fita,
       c("ci","pi"),
       combine_pi_ci = FALSE,
       variable = c("x1","x2","x3"))
  
  # Check if splines and factors work
  library(splines)
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5],
               x4 = rnorm(N),
               x5 = rnorm(N),
               x6 = rnorm(N),
               x7 = rnorm(N),
               x8 = rnorm(N),
               x9 = rnorm(N),
               x10 = rnorm(N))
  test_data$outcome = 
    rnorm(N,-1 + test_data$x1 + test_data$x1^2 + 2 * (test_data$x3 %in% c("d","e")) )
  fita = 
    bma_inference(outcome ~ ns(x1,df = 5) + x2 + x3,
                  data = test_data,
                  mc_draws = 5e3,
                  user.int = FALSE)
  fita
  summary(fita)
  predict(fita,
          newdata = 
            test_data[1,])$newdata
  plot(fita)
  
  rm(list = ls())
  
  
  # ANOVA -------------------------------------------------------------------
  
  # Create data
  
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rep(letters[1:5],N/5))
  test_data$outcome = 
    rnorm(N,-1 + 2 * (test_data$x1 %in% c("d","e")) )
  
  
  # Proper, heteroscedastic
  
  ## Make sure CI_level works (and print.lm_b works)
  fita = 
    aov_b(outcome ~ x1,
          test_data,
          prior_mean_mu = 2,
          prior_mean_nu = 0.5,
          prior_var_shape = 0.01,
          prior_var_rate = 0.01)
  
  fitb = 
    aov_b(outcome ~ x1,
          test_data,
          prior_mean_nu = 0.5,
          prior_var_shape = 0.01,
          prior_var_rate = 0.01,
          CI_level = 0.9)
  fita
  fitb
  
  ## Make sure summary.aov_b works
  summary(fita)
  
  ## Make sure coef.aov_b works
  coef(fita)
  
  ## Make sure credint works
  credint(fita)
  credint(fita,which = "mean")
  credint(fita, which = "pairw")
  
  ## Make sure vcov works
  vcov(fita)
  
  ## Make sure prior hyperparameters work
  fitc = 
    aov_b(outcome ~ x1,
          test_data,
          prior_mean_mu = 200,
          prior_mean_nu = 0.5,
          prior_var_shape = 0.01,
          prior_var_rate = 0.01,
          CI_level = 0.9)
  fita
  fitc
  
  ## Make sure prediction function works
  predict(fita)
  predict(fita,
          newdata = fita$data[1,])
  predict(fita,
          CI_level = 0.8)
  predict(fita,
          PI_level = 0.8)
  
  ## Make sure information criteria work
  AIC(fita)
  BIC(fita)
  DIC(fita)
  WAIC(fita)
  
  ## Make sure plotting function works
  plot(fita,
       type = "diagnostics")
  plot(fita,
       type = c("ci","pi"),
       combine_pi_ci = TRUE)
  plot(fita,
       type = "pi")
  plot(fita,
       type = "pi",
       PI_level = 0.8)
  plot(fita,
       type = "ci")
  plot(fita,
       type = "ci",
       CI_level = 0.999)
  plot(fita)
  
  
  ## Make sure parallelization works
  plan(multisession,workers = 5)
  fitd = 
    aov_b(outcome ~ x1,
          test_data,
          prior_mean_mu = 2,
          prior_mean_nu = 0.5,
          prior_var_shape = 0.01,
          prior_var_rate = 0.01)
  fitd
  fita
  plan(sequential)
  
  rm(list = setdiff(ls(),"test_data"))
  
  
  # Proper, homoscedastic
  
  ## Make sure CI_level works (and print.lm_b works)
  fita = 
    aov_b(outcome ~ x1,
          test_data,
          heteroscedastic = FALSE,
          prior_mean_mu = 2,
          prior_mean_nu = 0.5,
          prior_var_shape = 0.01,
          prior_var_rate = 0.01)
  
  fitb = 
    aov_b(outcome ~ x1,
          test_data,
          heteroscedastic = FALSE,
          prior_mean_mu = 2,
          prior_mean_nu = 0.5,
          prior_var_shape = 0.01,
          prior_var_rate = 0.01,
          CI_level = 0.8)
  fita
  fitb
  
  ## Make sure summary.aov_b works
  summary(fita)
  
  ## Make sure coef.aov_b works
  coef(fita)
  
  ## Make sure credint works
  credint(fita)
  
  ## Make sure vcov works
  vcov(fita)
  
  ## Make sure prior hyperparameters work
  fitc = 
    aov_b(outcome ~ x1,
          test_data,
          heteroscedastic = FALSE,
          prior_mean_mu = 200,
          prior_mean_nu = 0.5,
          prior_var_shape = 0.01,
          prior_var_rate = 0.01,
          CI_level = 0.9)
  fita
  fitc
  
  ## Make sure prediction function works
  predict(fita)
  predict(fita,
          CI_level = 0.8)
  predict(fita,
          PI_level = 0.8)
  
  ## Make sure information criteria work
  AIC(fita)
  BIC(fita)
  DIC(fita)
  WAIC(fita)
  
  rm(list = setdiff(ls(),"test_data"))
  
  
  
  # Improper, heteroscedastic
  
  ## Make sure CI_level works (and print.lm_b works)
  fita = 
    aov_b(outcome ~ x1,
          test_data,
          heteroscedastic = TRUE,
          improper = TRUE)
  
  fitb = 
    aov_b(outcome ~ x1,
          test_data,
          heteroscedastic = TRUE,
          improper = TRUE,
          CI_level = 0.8)
  fita
  fitb
  
  ## Make sure summary.aov_b works
  summary(fita)
  
  ## Make sure coef.aov_b works
  coef(fita)
  
  ## Make sure credint works
  credint(fita)
  
  ## Make sure vcov works
  vcov(fita)
  
  ## Make sure prior hyperparameters DO NOT work
  fitc = 
    aov_b(outcome ~ x1,
          test_data,
          heteroscedastic = TRUE,
          prior_mean_mu = 200,
          prior_mean_nu = 0.5,
          prior_var_shape = 0.01,
          prior_var_rate = 0.01,
          improper = TRUE)
  all.equal(summary(fita)$summary,
            summary(fitc)$summary)
  
  ## Make sure prediction function works
  predict(fita)
  predict(fita,
          newdata = fita$data[1,])
  predict(fita,
          CI_level = 0.8)
  predict(fita,
          PI_level = 0.8)
  
  
  ## Make sure information criteria work
  AIC(fita)
  BIC(fita)
  DIC(fita)
  WAIC(fita)
  
  rm(list = setdiff(ls(),"test_data"))
  
  
  
  # Improper, homoscedastic
  
  ## Make sure CI_level works (and print.lm_b works)
  fita = 
    aov_b(outcome ~ x1,
          test_data,
          heteroscedastic = FALSE,
          improper = TRUE)
  
  fitb = 
    aov_b(outcome ~ x1,
          test_data,
          heteroscedastic = FALSE,
          improper = TRUE,
          CI_level = 0.8)
  fita
  fitb
  
  ## Make sure summary.aov_b works
  summary(fita)
  
  ## Make sure coef.aov_b works
  coef(fita)
  
  ## Make sure vcov works
  vcov(fita)
  
  ## Make sure prior hyperparameters DO NOT work
  fitc = 
    aov_b(outcome ~ x1,
          test_data,
          heteroscedastic = FALSE,
          prior_mean_mu = 200,
          prior_mean_nu = 0.5,
          prior_var_shape = 0.01,
          prior_var_rate = 0.01,
          improper = TRUE)
  all.equal(summary(fita)$summary,
            summary(fitc)$summary)
  
  ## Make sure prediction function works
  predict(fita)
  predict(fita,
          CI_level = 0.8)
  predict(fita,
          PI_level = 0.8)
  
  ## Make sure information criteria work
  AIC(fita)
  BIC(fita)
  DIC(fita)
  WAIC(fita)
  
  
  
  # Test Bayes factor
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rep(letters[1:5],N/5))
  test_data$outcome = 
    rnorm(N,-1 + 2 * (test_data$x1 %in% c("d","e")) )
  
  set.seed(2025)
  test_data_null = 
    data.frame(x1 = rep(letters[1:5],N/5))
  test_data_null$outcome = 
    rnorm(N,-1)
  
  set.seed(2025)
  N = 25
  test_data_small = 
    data.frame(x1 = rep(letters[1:5],N/5))
  test_data_small$outcome = 
    rnorm(N,-1 + 2 * (test_data_small$x1 %in% c("d","e")) )
  test_data_null_small = 
    data.frame(x1 = rep(letters[1:5],N/5))
  test_data_null_small$outcome = 
    rnorm(N,-1)
  
  
  
  fita = 
    aov_b(outcome ~ x1,
          test_data,
          heteroscedastic = FALSE,
          prior_mean_mu = 2,
          prior_mean_nu = 0.5,
          prior_var_shape = 0.01,
          prior_var_rate = 0.01)
  fitb = 
    aov_b(outcome ~ x1,
          test_data_null,
          heteroscedastic = FALSE,
          prior_mean_mu = 2,
          prior_mean_nu = 0.5,
          prior_var_shape = 0.01,
          prior_var_rate = 0.01)
  fitc = 
    aov_b(outcome ~ x1,
          test_data_small,
          heteroscedastic = FALSE,
          prior_mean_mu = 2,
          prior_mean_nu = 0.5,
          prior_var_shape = 0.01,
          prior_var_rate = 0.01)
  fitd = 
    aov_b(outcome ~ x1,
          test_data_null_small,
          heteroscedastic = FALSE,
          prior_mean_mu = 2,
          prior_mean_nu = 0.5,
          prior_var_shape = 0.01,
          prior_var_rate = 0.01)
  summary(fita)
  summary(fitb)
  summary(fitc)
  summary(fitd)
  
  rm(list = ls())
  
  
  # Heteroscedasticity test for AOV -----------------------------------------
  
  
  # Test homoscedastic case
  set.seed(2025)
  N = 100
  test_data = 
    data.frame(x1 = rep(letters[1:5],N/5))
  test_data$outcome = 
    rnorm(N,-1 + 2 * (test_data$x1 %in% c("d","e")) )
  
  hetero_model = 
    aov_b(outcome ~ x1,
          test_data,
          prior_mean_mu = 2,
          prior_mean_nu = 0.5,
          prior_var_shape = 0.01,
          prior_var_rate = 0.01)
  homo_model = 
    aov_b(outcome ~ x1,
          test_data,
          heteroscedastic = FALSE,
          prior_mean_mu = 2,
          prior_mean_nu = 0.5,
          prior_var_shape = 0.01,
          prior_var_rate = 0.01)
  heteroscedasticity_test(hetero_model,
                          homo_model)
  
  # Test heteroscedastic case
  set.seed(2025)
  N = 100
  test_data = 
    data.frame(x1 = rep(letters[1:5],N/5))
  test_data$outcome = 
    rnorm(N,
          -1 + 2 * (test_data$x1 %in% c("d","e")),
          sd = 3 - 2 * (test_data$x1 %in% c("d","e")))
  
  hetero_model = 
    aov_b(outcome ~ x1,
          test_data,
          prior_mean_mu = 2,
          prior_mean_nu = 0.5,
          prior_var_shape = 0.01,
          prior_var_rate = 0.01)
  homo_model = 
    aov_b(outcome ~ x1,
          test_data,
          heteroscedastic = FALSE,
          prior_mean_mu = 2,
          prior_mean_nu = 0.5,
          prior_var_shape = 0.01,
          prior_var_rate = 0.01)
  heteroscedasticity_test(hetero_model,
                          homo_model)
  
  
  # Mediation ---------------------------------------------------------------
  
  # Binary treatment
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(tr = rep(0:1,N/2),
               x1 = rnorm(N))
  test_data$m = 
    rnorm(N, 0.4 * test_data$tr - 0.25 * test_data$x1)
  test_data$outcome = 
    rpois(N,exp(-1 + 0.6 * test_data$tr + 1.5 * test_data$m + 0.25 * test_data$x1))
  
  
  m1 = 
    lm_b(m ~ tr + x1,
         data = test_data)
  m2 = 
    glm_b(outcome ~ m + tr + x1,
          data = test_data,
          family = poisson())
  m3 = 
    mediate_b(m1,m2,
              treat = "tr",
              control_value = 0,
              treat_value = 1,
              n_draws = 500)
  m3
  summary(m3)
  plot(m3)
  plot(m3,type = "dx")
  plot(m3,type = "acme")
  plot(m3,type = "ade")
  plot(m3,type = c("acme","ade"))
  
  
  ## check answers 
  test_fit = 
    mediation::mediate(
      lm(m ~ tr, data = test_data),
      glm(outcome ~ m + tr, data = test_data, family = poisson()),
      treat = "tr",
      mediator = "m",
      control.value = 0,
      treat.value = 1
    )
  summary(test_fit)
  
  ## test batching
  mediate_b(m1,m2,
            treat = "tr",
            control_value = 0,
            treat_value = 1,
            n_draws = 1e3,
            batch_size = 1e3)
  
  ## Test interaction terms
  m2 = 
    glm_b(outcome ~ m * tr + x1,
          data = test_data,
          family = poisson())
  m3 = 
    mediate_b(m1,m2,
              treat = "tr",
              control_value = 0,
              treat_value = 1,
              n_draws = 500)
  m3
  summary(m3)
  
  m2 = 
    glm_b(outcome ~ m * x1 + tr,
          data = test_data,
          family = poisson())
  m3 = 
    mediate_b(m1,m2,
              treat = "tr",
              control_value = 0,
              treat_value = 1,
              n_draws = 500)
  m3
  summary(m3)
  
  
  # Test no mediation
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(tr = rep(0:1,N/2),
               x1 = rnorm(N))
  test_data$m = 
    rnorm(N,  -0.25 * test_data$x1)
  test_data$outcome = 
    rpois(N,exp(-1 + 0.6 * test_data$tr + 1.5 * test_data$m + 0.25 * test_data$x1))
  
  test_fit = 
    mediation::mediate(
      lm(m ~ tr, data = test_data),
      glm(outcome ~ m + tr, data = test_data, family = poisson()),
      treat = "tr",
      mediator = "m",
      control.value = 0,
      treat.value = 1
    )
  
  
  m1 = 
    lm_b(m ~ tr + x1,
         data = test_data)
  m2 = 
    glm_b(outcome ~ m + tr + x1,
          data = test_data,
          family = poisson())
  m3 = 
    mediate_b(m1,m2,
              treat = "tr",
              control_value = 0,
              treat_value = 1,
              ask_before_full_sampling = TRUE,
              n_draws = 1e3)
  
  summary(test_fit)
  m3
  
  
  # Continuous treatment
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(tr = rnorm(N),
               x1 = rnorm(N))
  test_data$m = 
    rnorm(N, 0.4 * test_data$tr - 0.25 * test_data$x1)
  test_data$outcome = 
    rnorm(N,-1 + 0.6 * test_data$tr + 1.5 * test_data$m + 0.25 * test_data$x1)
  
  m1 = 
    lm_b(m ~ tr + x1,
         data = test_data)
  m2 = 
    lm_b(outcome ~ m + tr + x1,
         data = test_data)
  m3 = 
    mediate_b(m1,m2,
              treat = "tr",
              control_value = -1.5,
              treat_value = 1.5)
  test_fit = 
    mediation::mediate(
      lm(m ~ tr, data = test_data),
      lm(outcome ~ m + tr + x1, data = test_data),
      treat = "tr",
      mediator = "m",
      control.value = -1.5,
      treat.value = 1.5
    )
  m3
  summary(test_fit)
  summary(m3)
  
  
  plot(m3)
  plot(m3,type = "dx")
  plot(m3,type = "acme")
  plot(m3,type = "ade")
  plot(m3,type = c("acme","ade"))
  
  
  ## Test interaction
  m2 = 
    lm_b(outcome ~ m * tr + x1,
         data = test_data)
  m3 = 
    mediate_b(m1,m2,
              treat = "tr",
              control_value = 0,
              treat_value = 1,
              n_draws = 500)
  m3
  summary(m3)
  
  m2 = 
    lm_b(outcome ~ m * x1 + tr,
         data = test_data)
  m3 = 
    mediate_b(m1,m2,
              treat = "tr",
              control_value = 0,
              treat_value = 1,
              n_draws = 500)
  m3
  summary(m3)
  
  rm(list=ls())
  
  
  # Test glm_b (binary) -----------------------------------------------------
  
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5])
  test_data$outcome = 
    rbinom(N,1,1.0 / (1.0 + exp(-(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) ))))
  
  table(test_data$outcome,test_data$x3) |>
    prop.table(2)
  boxplot(x1 ~ outcome,test_data)
  
  # VB
  ## Zellner
  fita =
    glm_b(outcome ~ x1 + x2 + x3,
          data = test_data,
          family = binomial(),
          seed = 2025)
  fitb =
    glm_b(outcome ~ x1 + x2 + x3,
          data = test_data,
          family = binomial(),
          seed = 2025,
          CI_level = 0.8)
  plan(multisession,workers = 10)
  fitc =
    glm_b(outcome ~ x1 + x2 + x3,
          data = test_data,
          seed = 2025,
          CI_level = 0.8)
  fitc
  plan(sequential)
  fita
  fitb
  coef(fita)
  credint(fita)
  summary(fita)
  summary(fita,
          CI_level = 0.8)
  summary(fita,
          interpretable_scale = FALSE)
  vcov(fita)
  preds = predict(fita)
  predict(fita,newdata = fita$data[1,])
  boxplot(`Post Mean` ~ outcome, data = preds)
  SDratio(fita)
  SDratio(fita, by = 'var')
  
  ## Make sure information criteria work
  null_model = 
    glm_b(outcome ~ 1,
          data = test_data,
          family = binomial(),
          seed = 2025)
  AIC(fita)
  AIC(null_model)
  BIC(fita)
  BIC(null_model)
  DIC(fita)
  DIC(null_model)
  WAIC(fita)
  WAIC(null_model)
  
  
  
  ## Test number of inputs
  fitd = 
    glm_b(test_data$outcome ~ test_data$x1,
          family = binomial(),
          prior = "normal")
  fitd
  fite = 
    glm_b(test_data$outcome ~ 1,
          family = binomial(),
          prior = "normal")
  fite
  fitf = 
    glm_b(outcome ~ x1,
          data = test_data,
          family = binomial(),
          prior = "normal")
  fitf
  fitg = 
    glm_b(outcome ~ 1,
          data = test_data,
          family = binomial(),
          prior = "normal")
  fitg
  
  ## Normal
  fitc =
    glm_b(outcome ~ x1 + x2 + x3,
          data = test_data,
          family = binomial(),
          prior = "normal",
          seed = 2025)
  fita
  fitc
  coef(fitc)
  summary(fitc)
  SDratio(fitc)
  SDratio(fitc,by="var")
  
  
  ## improper
  fitc =
    glm_b(outcome ~ x1 + x2 + x3,
          data = test_data,
          family = binomial(),
          prior = "improper")
  fita
  fitc
  coef(fitc)
  summary(fitc)
  vcov(fitc)
  try({
    SDratio(fitc)
    print("If you see this message, something went wrong.  SDratio SHOULD throw an error if used on an improper prior.")
  }, silent=T)
  
  
  
  # IS
  ## Zellner
  fita =
    glm_b(outcome ~ x1 + x2 + x3,
          data = test_data,
          family = binomial(),
          algorithm = "IS",
          seed = 2025)
  fitb =
    glm_b(outcome ~ x1 + x2 + x3,
          data = test_data,
          family = binomial(),
          algorithm = "IS",
          seed = 2025,
          CI_level = 0.8)
  fita
  fitb
  coef(fita)
  credint(fita)
  summary(fita)
  summary(fita,
          CI_level = 0.8)
  summary(fita,
          interpretable_scale = FALSE)
  vcov(fita)
  preds = predict(fita)
  predict(fita, newdata = fita$data[1,])
  boxplot(`Post Mean` ~ outcome, data = preds)
  SDratio(fita)
  SDratio(fita, by = 'var')
  
  ## Make sure information criteria work
  null_model = 
    glm_b(outcome ~ 1,
          data = test_data,
          family = binomial(),
          seed = 2025,
          algorithm = "IS")
  AIC(fita)
  AIC(null_model)
  BIC(fita)
  BIC(null_model)
  DIC(fita)
  DIC(null_model)
  WAIC(fita)
  WAIC(null_model)
  
  
  
  ## Normal
  fitc =
    glm_b(outcome ~ x1 + x2 + x3,
          data = test_data,
          family = binomial(),
          prior = "normal",
          seed = 2025)
  fita
  fitc
  coef(fitc)
  summary(fitc)
  SDratio(fitc)
  
  
  ## improper
  fitc =
    glm_b(outcome ~ x1 + x2 + x3,
          data = test_data,
          family = binomial(),
          prior = "improper")
  fita
  fitc
  coef(fitc)
  summary(fitc)
  try({
    SDratio(fitc)
    print("If you see this message, something went wrong.  SDratio SHOULD throw an error if used on an improper prior.")
  }, silent=T)
  
  
  
  # Large sample approx
  ## Zellner
  fitb =
    glm_b(outcome ~ x1 + x2 + x3,
          data = test_data,
          family = binomial(),
          algorithm = "LSA",
          seed = 2025)
  fitc =
    glm_b(outcome ~ x1 + x2 + x3,
          data = test_data,
          family = binomial(),
          algorithm = "LSA",
          seed = 2025,
          CI_level = 0.8)
  fita
  fitb
  fitc
  coef(fitb)
  credint(fitb)
  summary(fitb)
  vcov(fitb)
  preds = predict(fitb)
  predict(fitb,newdata=fitb$data[1,])
  boxplot(`Post Mean` ~ outcome, data = preds)
  SDratio(fitb)
  SDratio(fitb, by="variabl")
  
  ## Make sure information criteria work
  null_model = 
    glm_b(outcome ~ 1,
          data = test_data,
          family = binomial(),
          algorithm = "LSA",
          seed = 2025)
  AIC(fitb)
  AIC(null_model)
  BIC(fitb)
  BIC(null_model)
  DIC(fitb)
  DIC(null_model)
  WAIC(fitb)
  WAIC(null_model)
  
  
  
  
  ## Normal
  fitd =
    glm_b(outcome ~ x1 + x2 + x3,
          data = test_data,
          family = binomial(),
          prior = "normal",
          algorithm = "LSA",
          seed = 2025)
  fitd
  coef(fitb)
  summary(fitb)
  
  ## improper
  fite =
    glm_b(outcome ~ x1 + x2 + x3,
          data = test_data,
          family = binomial(),
          algorithm = "LSA",
          prior = "improper")
  fite
  coef(fite)
  summary(fite)
  
  # Make sure other link functions work, including failing to give a default rope
  fitf = 
    glm_b(outcome ~ x1 + x2 + x3,
          data = test_data,
          family = binomial("probit"),
          seed = 2025,
          algorithm = "LSA")
  fitf
  summary(fitf)
  fitg = 
    glm_b(outcome ~ x1 + x2 + x3,
          data = test_data,
          family = binomial("probit"),
          seed = 2025,
          algorithm = "LSA")
  fitg
  summary(fitg)
  
  
  
  # Check plotting functionality
  fita =
    glm_b(outcome ~ x1 + x2 + x3,
          data = test_data,
          family = binomial(),
          seed = 2025)
  fitb =
    glm_b(outcome ~ x1 + x2 + x3,
          data = test_data,
          family = binomial(),
          algorithm = "IS",
          seed = 2025)
  plot(fita,
       type = "diagnostics")
  plot(fitb,
       type = "diagnostics")
  plot(fita,
       type = "pdp")
  plot(fitb,
       type = "pdp")
  plot(fita,
       type = "pdp",
       variable = "x1")
  plot(fitb,
       type = "pdp",
       variable = "x1")
  plot(fita,
       type = c("ci","pi"),
       variable = "x1",
       combine_pi_ci = TRUE,
       exemplar_covariates = fita$data[1,])
  plot(fitb,
       type = c("ci","pi"),
       variable = "x1",
       combine_pi_ci = TRUE,
       exemplar_covariates = fita$data[1,])
  plot(fita,
       type = c("ci","pi"),
       variable = "x1",
       combine_pi_ci = FALSE,
       exemplar_covariates = fita$data[1,])
  plot(fitb,
       type = c("ci","pi"),
       variable = "x1",
       combine_pi_ci = FALSE,
       exemplar_covariates = fita$data[1,])
  plot(fita,
       type = c("ci","pi"),
       combine_pi_ci = TRUE,
       exemplar_covariates = fita$data[1,])
  plot(fitb,
       type = c("ci","pi"),
       combine_pi_ci = TRUE,
       exemplar_covariates = fita$data[1,])
  plot(fita,
       type = c("ci","pi"),
       combine_pi_ci = FALSE,
       exemplar_covariates = fita$data[1,])
  plot(fitb,
       type = c("ci","pi"),
       combine_pi_ci = FALSE,
       exemplar_covariates = fita$data[1,])
  plot(fita,
       type = c("ci","pi"),
       combine_pi_ci = TRUE)
  plot(fitb,
       type = c("ci","pi"),
       combine_pi_ci = TRUE)
  plot(fita,
       type = "pi",
       variable = "x1")
  plot(fitb,
       type = "pi",
       variable = "x1")
  plot(fita,
       type = "ci",
       variable = "x1")
  plot(fitb,
       type = "ci",
       variable = "x1")
  plot(fita,
       type = "ci")
  plot(fitb,
       type = "ci")
  plot(fita,
       type = "ci",
       CI_level = 0.999)
  plot(fitb,
       type = "ci",
       CI_level = 0.999)
  plot(fita)
  plot(fitb)
  
  
  # Check if splines and factors work
  library(splines)
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5])
  test_data$outcome = 
    rbinom(N,1,1.0 / (1.0 + exp(-(-2 + test_data$x1 + test_data$x1^2 + 2 * (test_data$x3 %in% c("d","e")) ))))
  fita = 
    glm_b(outcome ~ ns(x1,df = 5) + x2 + x3,
          data = test_data,
          family = binomial())
  fita
  summary(fita)
  predict(fita,
          newdata = 
            test_data[1,])
  plot(fita)
  
  
  # Check SDratio for complicated terms
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5])
  test_data$outcome = 
    rbinom(N,1,1.0 / (1.0 + exp(-(-2 + 
                                    test_data$x1 + 
                                    test_data$x1^2 + 
                                    2 * (test_data$x3 %in% c("d","e")) ))))
  SDratio(glm_b(outcome ~ x1 * x2 + x3 + I(x1^2),
                data = test_data,
                family = binomial()),
          by = "coef")
  SDratio(glm_b(outcome ~ x1 * x2 + x3 + I(x1^2),
                data = test_data,
                family = binomial()),
          by = "vari")
  library(splines)
  SDratio(glm_b(outcome ~ ns(x1,df = 5) + x2,
                data = test_data,
                family = binomial()),
          by = "coef")
  SDratio(glm_b(outcome ~ ns(x1,df = 5) + x2,
                data = test_data,
                family = binomial()),
          by = "var")
  
  
  
  # Test out trials
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5],
               n_trials = rpois(N,15))
  test_data$outcome = 
    rbinom(N,
           test_data$n_trials,
           1.0 / (1.0 + exp(-(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) ))))
  fita = 
    glm_b(outcome ~ x1 + x2 + x3,
          data = test_data,
          trials = "n_trials",
          family = binomial())
  fita
  fita = 
    glm_b(outcome ~ x1 + x2 + x3,
          data = test_data,
          trials = test_data$n_trials,
          family = binomial())
  fita
  
  
  rm(list=ls())
  
  
  
  # Test glm_b (Poisson) ----------------------------------------------------
  
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5],
               time = rexp(N))
  test_data$outcome = 
    rpois(N,exp(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e"))) * test_data$time)
  
  
  
  # VB
  ## Zellner
  fita =
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = poisson(),
          seed = 2025)
  fitb =
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = poisson(),
          seed = 2025,
          CI_level = 0.8)
  fita
  fitb
  coef(fita)
  credint(fita)
  summary(fita)
  summary(fita,
          CI_level = 0.8)
  summary(fita,
          interpretable_scale = FALSE)
  vcov(fita)
  preds = predict(fita)
  colnames(preds)
  predict(fita,newdata = fita$data[1:10,])
  plot(`Post Mean` ~ outcome, 
       data = preds |> dplyr::arrange(outcome))
  SDratio(fita)
  SDratio(fita, by = "Var")
  
  
  ## Make sure information criteria work
  null_model = 
    glm_b(outcome ~ 1 + offset(log(time)),
          data = test_data,
          family = poisson(),
          seed = 2025)
  AIC(fita)
  AIC(null_model)
  BIC(fita)
  BIC(null_model)
  DIC(fita)
  DIC(null_model)
  WAIC(fita)
  WAIC(null_model)
  
  
  
  ## Test number of inputs
  fitd = 
    glm_b(test_data$outcome ~ test_data$x1 + offset(log(test_data$time)),
          family = poisson(),
          prior = "normal")
  fitd
  fite = 
    glm_b(test_data$outcome ~ 1 + offset(log(test_data$time)),
          family = poisson(),
          prior = "normal")
  fite
  fitf = 
    glm_b(outcome ~ x1 + offset(log(time)),
          data = test_data,
          family = poisson(),
          prior = "normal")
  fitf
  fitg = 
    glm_b(outcome ~ 1 + offset(log(time)),
          data = test_data,
          family = poisson(),
          prior = "normal")
  fitg
  
  ## Normal
  fitc =
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = poisson(),
          prior = "normal",
          seed = 2025)
  fita
  fitc
  coef(fitc)
  credint(fitc)
  summary(fitc)
  
  ## improper
  fitc =
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = poisson(),
          prior = "improper")
  fita
  fitc
  coef(fitc)
  summary(fitc)
  try({
    SDratio(fitc)
    cat("Shouldn't see this message!")
  }, silent = TRUE)
  
  # IS
  ## Zellner
  fita =
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = poisson(),
          algorithm = "IS",
          seed = 2025)
  fitb =
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = poisson(),
          algorithm = "IS",
          seed = 2025,
          CI_level = 0.8)
  fita
  fitb
  coef(fita)
  credint(fita)
  summary(fita)
  summary(fita,
          CI_level = 0.8)
  summary(fita,
          interpretable_scale = FALSE)
  vcov(fita)
  preds = predict(fita)
  colnames(preds)
  predict(fita,newdata = fita$data[1,])
  plot(`Post Mean` ~ outcome, 
       data = preds |> dplyr::arrange(outcome))
  SDratio(fita)
  SDratio(fita, by = "Var")
  
  ## Test number of inputs
  fitd = 
    glm_b(test_data$outcome ~ test_data$x1 + offset(log(test_data$time)),
          family = poisson(),
          algorithm = "IS",
          prior = "normal")
  fitd
  fite = 
    glm_b(test_data$outcome ~ 1 + offset(log(test_data$time)),
          family = poisson(),
          algorithm = "IS",
          prior = "normal")
  fite
  
  ## Make sure information criteria work
  null_model = 
    glm_b(outcome ~ 1 + offset(log(time)),
          data = test_data,
          family = poisson(),
          algorithm = "IS",
          seed = 2025)
  AIC(fita)
  AIC(null_model)
  BIC(fita)
  BIC(null_model)
  DIC(fita)
  DIC(null_model)
  WAIC(fita)
  WAIC(null_model)
  
  ## Normal
  fitc =
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = poisson(),
          prior = "normal",
          algorithm = "IS",
          seed = 2025)
  fita
  fitc
  coef(fitc)
  summary(fitc)
  
  ## improper
  fitc =
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = poisson(),
          algorithm = "IS",
          prior = "improper")
  fita
  fitc
  coef(fitc)
  summary(fitc)
  try({
    SDratio(fitc)
    print("If you see this message, something went wrong.  SDratio SHOULD throw an error if used on an improper prior.")
  }, silent=T)
  
  
  
  # Large sample approx
  ## Zellner
  fitb =
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = poisson(),
          algorithm = "LSA",
          seed = 2025)
  fitc =
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = poisson(),
          algorithm = "LSA",
          seed = 2025,
          CI_level = 0.8)
  fita
  fitb
  fitc
  coef(fitb)
  credint(fitb)
  summary(fitb)
  summary(fitb,
          CI_level = 0.8)
  summary(fitb,
          interpretable_scale = FALSE)
  vcov(fitb)
  preds = predict(fitb)
  str(preds)
  predict(fitb,newdata=fitb$data[1,])
  plot(`Post Mean` ~ outcome, data = preds |> dplyr::arrange(outcome))
  SDratio(fitb)
  SDratio(fitb,by="v")
  
  ## Test number of inputs
  fitd = 
    glm_b(test_data$outcome ~ test_data$x1 + offset(log(test_data$time)),
          family = poisson(),
          algorithm = "LSA",
          prior = "normal")
  fitd
  fite = 
    glm_b(test_data$outcome ~ 1 + offset(log(test_data$time)),
          family = poisson(),
          algorithm = "LSA",
          prior = "normal")
  fite
  
  
  ## Make sure information criteria work
  null_model = 
    glm_b(outcome ~ 1 + offset(log(time)),
          data = test_data,
          family = poisson(),
          algorithm = "LSA",
          seed = 2025)
  AIC(fitb)
  AIC(null_model)
  BIC(fitb)
  BIC(null_model)
  DIC(fitb)
  DIC(null_model)
  WAIC(fitb)
  WAIC(null_model)
  
  
  
  ## Normal
  fitd =
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = poisson(),
          prior = "normal",
          algorithm = "LSA",
          seed = 2025)
  fitd
  coef(fitb)
  summary(fitb)
  
  ## improper
  fite =
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = poisson(),
          algorithm = "LSA",
          prior = "improper")
  fite
  coef(fite)
  summary(fite)
  
  
  
  # Check plotting functionality
  fita =
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = poisson(),
          seed = 2025)
  fitb =
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = poisson(),
          seed = 2025,
          algorithm = "IS")
  fitc =
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = poisson(),
          seed = 2025,
          algorithm = "LSA")
  plot(fita,
       type = "diagnostics")
  plot(fitb,
       type = "diagnostics")
  plot(fitc,
       type = "diagnostics")
  plot(fita,
       type = "pdp")
  plot(fitb,
       type = "pdp")
  plot(fitc,
       type = "pdp")
  plot(fita,
       type = "pdp",
       variable = "x1")
  plot(fitb,
       type = "pdp",
       variable = "x1")
  plot(fita,
       type = c("ci","pi"),
       variable = "x1",
       combine_pi_ci = TRUE,
       exemplar_covariates = fita$data[1,])
  plot(fitb,
       type = c("ci","pi"),
       variable = "x1",
       combine_pi_ci = TRUE,
       exemplar_covariates = fita$data[1,])
  plot(fita,
       type = c("ci","pi"),
       combine_pi_ci = TRUE,
       exemplar_covariates = fita$data[1,])
  plot(fitb,
       type = c("ci","pi"),
       combine_pi_ci = TRUE,
       exemplar_covariates = fita$data[1,])
  plot(fita,
       type = c("ci","pi"),
       variable = "x1",
       combine_pi_ci = FALSE,
       exemplar_covariates = fita$data[1,])
  plot(fitb,
       type = c("ci","pi"),
       variable = "x1",
       combine_pi_ci = FALSE,
       exemplar_covariates = fita$data[1,])
  plot(fita,
       type = c("ci","pi"),
       combine_pi_ci = TRUE)
  plot(fitb,
       type = c("ci","pi"),
       combine_pi_ci = TRUE)
  plot(fita,
       type = c("ci","pi"),
       combine_pi_ci = FALSE,
       exemplar_covariates = fita$data[1,])
  plot(fitb,
       type = c("ci","pi"),
       combine_pi_ci = FALSE,
       exemplar_covariates = fita$data[1,])
  plot(fita,
       type = "pi",
       variable = "x1")
  plot(fitb,
       type = "pi",
       variable = "x1")
  plot(fita,
       type = "pi",
       variable = "x1",
       PI_level = 0.5)
  plot(fitb,
       type = "pi",
       variable = "x1",
       PI_level = 0.5)
  plot(fita,
       type = "pi")
  plot(fitb,
       type = "pi")
  plot(fita,
       type = "pi",
       exemplar_covariates = fita$data[1,])
  plot(fitb,
       type = "pi",
       exemplar_covariates = fita$data[1,])
  plot(fita,
       type = "ci",
       variable = "x1")
  plot(fitb,
       type = "ci",
       variable = "x1")
  plot(fita,
       type = "ci")
  plot(fitb,
       type = "ci")
  plot(fita,
       type = "ci",
       CI_level = 0.999)
  plot(fitb,
       type = "ci",
       CI_level = 0.999)
  plot(fita)
  plot(fitb)
  
  
  
  
  # Check if splines and factors work
  library(splines)
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5],
               time = rexp(N))
  test_data$outcome = 
    rpois(N,exp(-2 + test_data$x1 + test_data$x1^2 + 2 * (test_data$x3 %in% c("d","e"))) * test_data$time)
  fita = 
    glm_b(outcome ~ ns(x1,df = 5) + x2 + x3,
          data = test_data,
          family = poisson())
  fita
  summary(fita)
  predict(fita,
          newdata = 
            test_data[1,])
  plot(fita)
  SDratio(fita)
  SDratio(fita,by="v")
  
  
  rm(list=ls())
  
  
  
  
  # Test glm_b (NegBinom) ---------------------------------------------------
  
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5],
               time = rexp(N))
  test_data$outcome = 
    rnbinom(N,
            mu = exp(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e"))) * test_data$time,
            size = 0.7)
  
  # VB
  ## Zellner
  fita =
    glm_b(formula = outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = negbinom(),
          seed = 2025)
  fitb =
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = negbinom(),
          seed = 2025,
          CI_level = 0.8)
  fita
  fitb
  coef(fita)
  credint(fita)
  summary(fita)
  summary(fita,
          CI_level = 0.8)
  summary(fita,
          interpretable_scale = FALSE)
  vcov(fita)
  predict(fita,
          newdata = fita$data[1,])
  preds = predict(fita)
  colnames(preds)
  plot(`Post Mean` ~ outcome, 
       data = preds |> dplyr::arrange(outcome))
  SDratio(fita)
  SDratio(fita, by = "var")
  
  
  ## Make sure information criteria work
  null_model = 
    glm_b(outcome ~ 1 + offset(log(time)),
          data = test_data,
          family = negbinom(),
          seed = 2025)
  AIC(fita)
  AIC(null_model)
  BIC(fita)
  BIC(null_model)
  DIC(fita)
  DIC(null_model)
  WAIC(fita)
  WAIC(null_model)
  
  
  
  ## Test number of inputs
  fitd = 
    glm_b(test_data$outcome ~ test_data$x1 + offset(log(test_data$time)),
          family = negbinom(),
          prior = "normal")
  fitd
  fite = 
    glm_b(test_data$outcome ~ 1 + offset(log(test_data$time)),
          family = negbinom(),
          prior = "normal")
  fite
  fitf = 
    glm_b(outcome ~ x1 + offset(log(time)),
          data = test_data,
          family = negbinom(),
          prior = "normal")
  fitf
  fitg = 
    glm_b(outcome ~ 1 + offset(log(time)),
          data = test_data,
          family = negbinom(),
          prior = "normal")
  fitg
  
  ## Normal
  fitc =
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = negbinom(),
          prior = "normal",
          seed = 2025)
  fita
  fitc
  coef(fitc)
  credint(fitc)
  summary(fitc)
  
  ## improper
  fitc =
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = negbinom(),
          prior = "improper")
  fita
  fitc
  coef(fitc)
  summary(fitc)
  try({
    SDratio(fitc)
    print("If you see this message, something went wrong.  SDratio SHOULD throw an error if used on an improper prior.")
  }, silent=T)
  
  
  # IS
  ## Zellner
  fita =
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = negbinom(),
          algorithm = "IS",
          seed = 2025)
  fitb =
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = negbinom(),
          algorithm = "IS",
          seed = 2025,
          CI_level = 0.8)
  fita
  fitb
  coef(fita)
  credint(fita)
  summary(fita)
  summary(fita,
          CI_level = 0.8)
  summary(fita,
          interpretable_scale = FALSE)
  vcov(fita)
  predict(fita,
          newdata = fita$data[1,])
  preds = predict(fita)
  colnames(preds)
  plot(`Post Mean` ~ outcome, 
       data = preds |> dplyr::arrange(outcome))
  SDratio(fita)
  SDratio(fita,by='v')
  
  ## Test number of inputs
  fitd = 
    glm_b(test_data$outcome ~ test_data$x1 + offset(log(test_data$time)),
          family = negbinom(),
          algorithm = "IS",
          prior = "normal")
  fitd
  fite = 
    glm_b(test_data$outcome ~ 1 + offset(log(test_data$time)),
          family = negbinom(),
          algorithm = "IS",
          prior = "normal")
  fite
  
  ## Make sure information criteria work
  null_model = 
    glm_b(outcome ~ 1 + offset(log(time)),
          data = test_data,
          family = negbinom(),
          algorithm = "IS",
          seed = 2025)
  AIC(fita)
  AIC(null_model)
  BIC(fita)
  BIC(null_model)
  DIC(fita)
  DIC(null_model)
  WAIC(fita)
  WAIC(null_model)
  
  ## Normal
  fitc =
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = negbinom(),
          prior = "normal",
          algorithm = "IS",
          seed = 2025)
  fita
  fitc
  coef(fitc)
  summary(fitc)
  
  ## improper
  fitc =
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = negbinom(),
          algorithm = "IS",
          prior = "improper")
  fita
  fitc
  coef(fitc)
  summary(fitc)
  try({
    SDratio(fitc)
    print("If you see this message, something went wrong.  SDratio SHOULD throw an error if used on an improper prior.")
  }, silent=T)
  
  
  
  # Large sample approx
  ## Zellner
  fitb =
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = negbinom(),
          algorithm = "LSA",
          seed = 2025)
  fitc =
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = negbinom(),
          algorithm = "LSA",
          seed = 2025,
          CI_level = 0.8)
  glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
        data = test_data,
        family = "negbinom",
        algorithm = "LSA",
        seed = 2025)
  glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
        data = test_data,
        algorithm = "LSA",
        seed = 2025)
  fita
  fitb
  fitc
  coef(fitb)
  credint(fitb)
  summary(fitb)
  summary(fitb,
          CI_level = 0.8)
  summary(fitb,
          interpretable_scale = FALSE)
  vcov(fitb)
  predict(fitb,
          newdata = fitb$data[1,])
  preds = predict(fitb)
  str(preds)
  plot(`Post Mean` ~ outcome, data = preds |> dplyr::arrange(outcome))
  SDratio(fitb)
  SDratio(fitb,by = "V")
  
  ## Test number of inputs
  fitd = 
    glm_b(test_data$outcome ~ test_data$x1 + offset(log(test_data$time)),
          family = negbinom(),
          algorithm = "LSA",
          prior = "normal")
  fitd
  fite = 
    glm_b(test_data$outcome ~ 1 + offset(log(test_data$time)),
          family = negbinom(),
          algorithm = "LSA",
          prior = "normal")
  fite
  
  
  ## Make sure information criteria work
  null_model = 
    glm_b(outcome ~ 1 + offset(log(time)),
          data = test_data,
          family = negbinom(),
          algorithm = "LSA",
          seed = 2025)
  AIC(fitb)
  AIC(null_model)
  BIC(fitb)
  BIC(null_model)
  DIC(fitb)
  DIC(null_model)
  WAIC(fitb)
  WAIC(null_model)
  
  
  
  ## Normal
  fitd =
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = negbinom(),
          prior = "normal",
          algorithm = "LSA",
          seed = 2025)
  fitd
  coef(fitb)
  summary(fitb)
  
  ## improper
  fite =
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = negbinom(),
          algorithm = "LSA",
          prior = "improper")
  fite
  coef(fite)
  summary(fite)
  
  
  
  # Check plotting functionality
  fita =
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = negbinom(),
          seed = 2025)
  fitb =
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = negbinom(),
          seed = 2025,
          algorithm = "IS")
  plot(fita,
       type = "diagnostics")
  plot(fitb,
       type = "diagnostics")
  plot(fita,
       type = "pdp")
  plot(fitb,
       type = "pdp")
  plot(fitc,
       type = "pdp")
  plot(fita,
       type = "pdp",
       variable = "x1")
  plot(fitb,
       type = "pdp",
       variable = "x1")
  plot(fita,
       type = c("ci","pi"),
       variable = "x1",
       combine_pi_ci = TRUE,
       exemplar_covariates = fita$data[1,])
  plot(fitb,
       type = c("ci","pi"),
       variable = "x1",
       combine_pi_ci = TRUE,
       exemplar_covariates = fita$data[1,])
  plot(fita,
       type = c("ci","pi"),
       combine_pi_ci = TRUE,
       exemplar_covariates = fita$data[1,])
  plot(fitb,
       type = c("ci","pi"),
       combine_pi_ci = TRUE,
       exemplar_covariates = fita$data[1,])
  plot(fita,
       type = c("ci","pi"),
       variable = "x1",
       combine_pi_ci = FALSE,
       exemplar_covariates = fita$data[1,])
  plot(fitb,
       type = c("ci","pi"),
       variable = "x1",
       combine_pi_ci = FALSE,
       exemplar_covariates = fita$data[1,])
  plot(fita,
       type = c("ci","pi"),
       combine_pi_ci = TRUE)
  plot(fitb,
       type = c("ci","pi"),
       combine_pi_ci = TRUE)
  plot(fita,
       type = c("ci","pi"),
       combine_pi_ci = FALSE,
       exemplar_covariates = fita$data[1,])
  plot(fitb,
       type = c("ci","pi"),
       combine_pi_ci = FALSE,
       exemplar_covariates = fita$data[1,])
  plot(fita,
       type = "pi",
       variable = "x1")
  plot(fitb,
       type = "pi",
       variable = "x1")
  plot(fita,
       type = "pi",
       variable = "x1",
       PI_level = 0.5)
  plot(fitb,
       type = "pi",
       variable = "x1",
       PI_level = 0.5)
  plot(fita,
       type = "pi")
  plot(fitb,
       type = "pi")
  plot(fita,
       type = "pi",
       exemplar_covariates = fita$data[1,])
  plot(fitb,
       type = "pi",
       exemplar_covariates = fita$data[1,])
  plot(fita,
       type = "ci",
       variable = "x1")
  plot(fitb,
       type = "ci",
       variable = "x1")
  plot(fita,
       type = "ci")
  plot(fitb,
       type = "ci")
  plot(fita,
       type = "ci",
       CI_level = 0.999)
  plot(fitb,
       type = "ci",
       CI_level = 0.999)
  plot(fita)
  plot(fitb)
  
  
  
  
  # Check if splines and factors work
  library(splines)
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5],
               time = rexp(N))
  test_data$outcome = 
    rnbinom(N,
            mu = exp(-2 + test_data$x1 + test_data$x1^2 + 2 * (test_data$x3 %in% c("d","e"))) * test_data$time,
            size = 0.7)
  fita = 
    glm_b(outcome ~ ns(x1,df = 5) + x2 + x3,
          data = test_data,
          family = negbinom())
  fita
  summary(fita)
  predict(fita,
          newdata = 
            test_data[1,])
  plot(fita)
  
  
  
  rm(list=ls())
  
  
  # Test glm_b (Gaussian) ----------------------------------------------------
  
  
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5])
  test_data$outcome = 
    rnorm(N,-1 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) )
  
  fita = 
    glm_b(outcome ~ x1 + x2 + x3,
          data = test_data,
          prior = "normal",
          family = "gaussian")
  fitb = 
    glm_b(outcome ~ x1 + x2 + x3,
          data = test_data,
          prior = "zellner",
          family = "gaussian")
  fitb = 
    glm_b(outcome ~ x1 + x2 + x3,
          data = test_data,
          prior = "zellner",
          zellner_g = 10,
          family = "gaussian")
  
  
  rm(list = ls())
  
  
  # Test loss-likelihood approach (Gaussian) --------------------------------
  
  # Create data
  
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5])
  test_data$outcome = 
    rnorm(N,-1 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) )
  
  
  # Bootstrapping approach - sequential
  plan(sequential)
  ## Make sure CI_level works (and print.lm_b works)
  fita = 
    np_glm_b(outcome ~ x1 + x2 + x3,
             data = test_data,
             family = gaussian(),
             n_draws = 50,
             seed = 2025)
  fitb = 
    np_glm_b(outcome ~ x1 + x2 + x3,
             data = test_data,
             family = gaussian(),
             n_draws = 50,
             seed = 2025,
             CI_level = 0.8)
  fita
  fitb
  
  ## Make sure summary.np_glm_b works
  summary(fita)
  
  ## Make sure coef.lm_b works
  coef(fita)
  
  ## Make sure credint works
  credint(fita)
  
  ## Make sure vcov works
  vcov(fita)
  
  ## Make sure prediction function works
  head(predict(fita))
  predict(fita,
          newdata = 
            fita$data[1,])
  
  ## Test number of inputs
  fitd = 
    np_glm_b(test_data$outcome ~ test_data$x1,
             family = gaussian(),
             n_draws = 100)
  fitd
  fite = 
    np_glm_b(test_data$outcome ~ 1,
             family = gaussian(),
             n_draws = 100)
  fite
  fitf = 
    np_glm_b(outcome ~ x1,
             data = test_data,
             family = gaussian(),
             n_draws = 100)
  fitf
  fitg = 
    np_glm_b(outcome ~ 1,
             data = test_data,
             family = gaussian(),
             n_draws = 100)
  fitg
  
  
  
  # Bootstrapping approach - parallel
  plan(multisession, workers = 10)
  ## Make sure future.seed is reproducible, CI_level works (and print.lm_b works)
  fita = 
    np_glm_b(outcome ~ x1 + x2 + x3,
             data = test_data,
             family = gaussian(),
             n_draws = 100,
             seed = 2025)
  fitb = 
    np_glm_b(outcome ~ x1 + x2 + x3,
             data = test_data,
             family = gaussian(),
             n_draws = 100,
             seed = 2025)
  fitc = 
    np_glm_b(outcome ~ x1 + x2 + x3,
             data = test_data,
             family = gaussian(),
             n_draws = 100,
             seed = 2025,
             CI_level = 0.8)
  all.equal(fita$summary,
            fitb$summary)
  fita
  fitc
  plan(sequential)
  
  ## Make sure summary.np_glm_b works
  summary(fita)
  
  ## Make sure coef.lm_b works
  coef(fita)
  
  ## Make sure prediction function works
  head(predict(fita))
  
  ## Check that parallelization is speeding things up
  
  system.time({
    fitd = 
      np_glm_b(outcome ~ x1 + x2 + x3,
               data = test_data,
               family = gaussian(),
               n_draws = 500,
               seed = 2025)
  })
  # user  system elapsed 
  # 11.76    6.28   23.77 
  
  plan(sequential)
  system.time({
    fite = 
      np_glm_b(outcome ~ x1 + x2 + x3,
               data = test_data,
               family = gaussian(),
               n_draws = 500,
               seed = 2025)
  })
  # user  system elapsed 
  # 6.84    0.07    7.00 # Eh?  Not sure what's happening.
  fitd
  fite
  
  plan(sequential)
  
  
  rm(list = setdiff(ls(),c("fita","test_data")))
  
  
  
  
  # large sample approach
  ## Make sure CI_level works (and print.lm_b works)
  fitb = 
    np_glm_b(outcome ~ x1 + x2 + x3,
             data = test_data,
             family = gaussian(),
             seed = 2025)
  fitc = 
    np_glm_b(outcome ~ x1 + x2 + x3,
             data = test_data,
             family = gaussian(),
             seed = 2025,
             CI_level = 0.8)
  fita
  fitb
  fitc
  
  
  ## Test number of inputs
  fitd = 
    np_glm_b(test_data$outcome ~ test_data$x1,
             family = gaussian())
  fitd
  fite = 
    np_glm_b(test_data$outcome ~ 1,
             family = gaussian())
  fite
  fitf = 
    np_glm_b(outcome ~ x1,
             data = test_data,
             family = gaussian())
  fitf
  fitg = 
    np_glm_b(outcome ~ 1,
             data = test_data,
             family = gaussian())
  fitg
  
  ## Make sure summary.np_glm_b works
  summary(fitb)
  
  ## Make sure coef.lm_b works
  coef(fitb)
  
  ## Make sure credint works
  credint(fitb)
  
  ## Make sure vcov works
  vcov(fitb)
  
  ## Make sure prediction function works
  head(predict(fitb))
  predict(fitb,
          newdata = fitb$data[1,])
  
  # Check plotting functionality
  plot(fita,
       type = "pdp")
  plot(fitb,
       type = "pdp")
  plot(fita,
       type = "pdp",
       variable = "x1")
  plot(fitb,
       type = "pdp",
       variable = "x1")
  plot(fita,
       type = "pdp",
       variable = "x3")
  plot(fitb,
       type = "pdp",
       variable = "x3")
  plot(fita,
       type = "ci",
       variable = "x1")
  plot(fitb,
       type = "ci",
       variable = "x1")
  plot(fita,
       type = "ci")
  plot(fita,
       type = "ci",
       CI_level = 0.999)
  plot(fitb,
       type = "ci")
  plot(fitb,
       type = "ci",
       CI_level = 0.999)
  plot(fita)
  plot(fitb)
  
  
  
  # Test GLS
  ## Bootstrapping
  fitd = 
    np_glm_b(outcome ~ x1 + x2 + x3,
             data = test_data,
             family = gaussian(),
             loss = "gls",
             n_draws = 50,
             seed = 2025)
  fitd
  summary(fitd)
  coef(fitd)
  vcov(fitd)
  preds0d = 
    predict(fitd)
  head(preds0d)
  
  ## large sample approx
  fite = 
    np_glm_b(outcome ~ x1 + x2 + x3,
             data = test_data,
             family = gaussian(),
             loss = "gls",
             seed = 2025)
  fitd
  fite
  summary(fite)
  coef(fite)
  vcov(fite)
  head(predict(fite))
  predict(fite,
          newdata = fite$data[1,])
  
  
  
  # Check if splines and factors work
  library(splines)
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5],
               time = rexp(N))
  test_data$outcome = 
    rnorm(N,-1 + test_data$x1 + test_data$x1^2 + 2 * (test_data$x3 %in% c("d","e")) )
  fita = 
    np_glm_b(outcome ~ ns(x1,df = 5) + x2 + x3,
             data = test_data,
             n_draws = 100,
             family = gaussian())
  fita
  summary(fita)
  predict(fita,
          newdata = 
            test_data[1,])
  plot(fita)
  
  fitb = 
    np_glm_b(outcome ~ ns(x1,df = 5) + x2 + x3,
             data = test_data,
             family = gaussian())
  fitb
  summary(fitb)
  predict(fitb,
          newdata = 
            test_data[1,])
  plot(fitb)
  
  
  rm(list = ls())
  
  
  # Test loss-likelihood approach (Binomial) --------------------------------
  
  # Create data
  
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5])
  test_data$outcome = 
    rbinom(N,1,1.0 / (1.0 + exp(-(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) ))))
  table(test_data$outcome,test_data$x3) |>
    prop.table(2)
  boxplot(x1 ~ outcome,test_data)
  
  # Bootstrapping approach - sequential
  plan(sequential)
  ## Make sure CI_level works (and print.lm_b works)
  fita = 
    np_glm_b(outcome ~ x1 + x2 + x3,
             data = test_data,
             family = binomial(),
             n_draws = 50,
             seed = 2025)
  fitb = 
    np_glm_b(outcome ~ x1 + x2 + x3,
             data = test_data,
             family = binomial(),
             n_draws = 50,
             trials = rep(1,N),
             seed = 2025,
             CI_level = 0.8)
  fita
  fitb
  
  ## Make sure summary.np_glm_b works
  summary(fita)
  
  ## Make sure coef.lm_b works
  coef(fita)
  
  ## Make sure vcov works
  vcov(fita)
  
  ## Make sure credint works
  credint(fita)
  
  ## Make sure prediction function works
  head(predict(fita))
  
  ## Test number of inputs
  fitd = 
    np_glm_b(test_data$outcome ~ test_data$x1,
             family = binomial(),
             n_draws = 50)
  fitd
  fite = 
    np_glm_b(test_data$outcome ~ 1,
             family = binomial(),
             n_draws = 50)
  fite
  fitf = 
    np_glm_b(outcome ~ x1,
             data = test_data,
             family = binomial(),
             n_draws = 50)
  fitf
  fitg = 
    np_glm_b(outcome ~ 1,
             data = test_data,
             family = binomial(),
             n_draws = 50)
  fitg
  
  
  
  # large sample approach
  ## Make sure CI_level works (and print.lm_b works)
  fitb = 
    np_glm_b(outcome ~ x1 + x2 + x3,
             data = test_data,
             family = binomial(),
             seed = 2025)
  fitc = 
    np_glm_b(outcome ~ x1 + x2 + x3,
             data = test_data,
             family = binomial(),
             seed = 2025,
             CI_level = 0.8)
  fita
  fitb
  fitc
  
  
  ## Test number of inputs
  fitd = 
    np_glm_b(test_data$outcome ~ test_data$x1,
             family = binomial())
  fitd
  fite = 
    np_glm_b(test_data$outcome ~ 1,
             family = binomial())
  fite
  fitf = 
    np_glm_b(outcome ~ x1,
             data = test_data,
             family = binomial())
  fitf
  fitg = 
    np_glm_b(outcome ~ 1,
             data = test_data,
             family = binomial())
  fitg
  
  
  ## Make sure summary.np_glm_b works
  summary(fitb)
  
  ## Make sure coef.lm_b works
  coef(fitb)
  
  ## Make sure credint works
  credint(fitb)
  
  ## Make sure vcov works
  vcov(fitb)
  
  ## Make sure prediction function works
  head(predict(fitb))
  
  
  # Make sure parallelization works
  plan(multisession,
       workers = 10)
  ## Make sure CI_level works (and print.lm_b works)
  fitd = 
    np_glm_b(outcome ~ x1 + x2 + x3,
             data = test_data,
             family = binomial(),
             n_draws = 100,
             seed = 2025)
  fita
  fitd
  summary(fitd)
  coef(fitd)
  preds0d = 
    predict(fitd)
  head(preds0d)
  
  
  # Make sure gls works
  ## Bootstrapping - parallelized
  fite = 
    np_glm_b(outcome ~ x1 + x2 + x3,
             data = test_data,
             family = binomial(),
             loss = "gls",
             n_draws = 50,
             seed = 2025)
  fitd
  fite
  summary(fite)
  coef(fite)
  vcov(fite)
  head(predict(fite))
  
  
  ## Bootstrapping - sequential
  plan(sequential)
  fitf = 
    np_glm_b(outcome ~ x1 + x2 + x3,
             data = test_data,
             family = binomial(),
             loss = "gls",
             n_draws = 250,
             seed = 2025)
  fite
  fitf
  summary(fitf)
  coef(fitf)
  preds0f = 
    predict(fitf)
  head(preds0f)
  
  ## Large sample approx
  fitg = 
    np_glm_b(outcome ~ x1 + x2 + x3,
             data = test_data,
             family = binomial(),
             loss = "gls",
             seed = 2025)
  fitf
  fitg
  summary(fitg)
  credint(fitg)
  coef(fitg)
  head(predict(fitg))
  
  
  # Check if splines and factors work
  library(splines)
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5],
               time = rexp(N))
  test_data$outcome = 
    rbinom(N,1,1.0 / (1.0 + exp(-(-2 + test_data$x1 + test_data$x1^2 + 2 * (test_data$x3 %in% c("d","e")) ))))
  fita = 
    np_glm_b(outcome ~ ns(x1,df = 5) + x2 + x3,
             data = test_data,
             n_draws = 50,
             family = binomial())
  fita
  summary(fita)
  predict(fita,
          newdata = 
            test_data[1,])
  plot(fita,type=c("pdp","ci","pi"))
  plot(fita)
  
  fitb = 
    np_glm_b(outcome ~ ns(x1,df = 5) + x2 + x3,
             data = test_data,
             family = binomial())
  fitb
  summary(fitb)
  predict(fitb,
          newdata = 
            test_data[1,])
  plot(fitb)
  
  
  
  rm(list = ls())
  
  
  # Test loss-likelihood approach (Poisson) ---------------------------------
  
  # Create data
  
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5],
               time = rexp(N))
  test_data$outcome = 
    rpois(N,exp(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e"))) * test_data$time)
  
  
  # Bootstrapping approach - sequential
  plan(sequential)
  ## Make sure CI_level works (and print.lm_b works)
  fita = 
    np_glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
             data = test_data,
             family = poisson(),
             n_draws = 50,
             seed = 2025)
  fitb = 
    np_glm_b(outcome ~ x1 + x2 + x3,
             data = test_data,
             family = poisson(),
             n_draws = 50,
             trials = rep(1,N),
             seed = 2025,
             CI_level = 0.8)
  fita
  fitb
  
  ## Make sure summary.np_glm_b works
  summary(fita)
  summary(fita,
          CI_level = 0.8)
  
  ## Make sure coef.lm_b works
  coef(fita)
  glm(outcome ~ x1 + x2 + x3 + offset(log(time)),
      data = test_data,
      family = poisson())
  
  ## Make sure credint works
  credint(fita)
  
  ## Make sure vcov works
  vcov(fita)
  
  ## Make sure prediction function works
  head(predict(fita))
  predict(fita,
          newdata = fita$data[1,])
  
  ## Test number of inputs
  fitd = 
    np_glm_b(test_data$outcome ~ test_data$x1 + offset(log(test_data$time)),
             family = poisson(),
             n_draws = 50)
  fitd
  fite = 
    np_glm_b(test_data$outcome ~ 1 + offset(log(test_data$time)),
             family = poisson(),
             n_draws = 50)
  fite
  fitf = 
    np_glm_b(outcome ~ x1 + offset(log(time)),
             data = test_data,
             family = poisson(),
             n_draws = 50)
  fitf
  fitg = 
    np_glm_b(outcome ~ 1,
             data = test_data,
             family = poisson(),
             n_draws = 50)
  fitg
  
  
  
  
  
  
  # Bootstrapping approach - parallelized
  plan(multisession, workers = 10)
  ## Make sure CI_level works (and print.lm_b works)
  fitc = 
    np_glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
             data = test_data,
             family = poisson(),
             n_draws = 250,
             seed = 2025)
  fita
  fitc
  summary(fitc)
  coef(fitc)
  head(predict(fitc))
  
  
  
  # large sample approach
  ## Make sure CI_level works (and print.lm_b works)
  fitd = 
    np_glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
             data = test_data,
             family = poisson(),
             seed = 2025)
  fite = 
    np_glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
             data = test_data,
             family = poisson(),
             seed = 2025,
             CI_level = 0.8)
  fita
  fitd
  fite
  
  ## Make sure summary.np_glm_b works
  summary(fitd)
  summary(fitd,
          CI_level = 0.8)
  
  ## Make sure coef.lm_b works
  coef(fitd)
  
  ## Make sure prediction function works
  head(predict(fitd))
  
  ## Make sure credint works
  credint(fitd)
  
  ## Make sure vcov works
  vcov(fitd)
  
  
  
  # Make sure gls works
  ## Bootstrapping - parallelized
  fite = 
    np_glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
             data = test_data,
             family = poisson(),
             loss = "gls",
             n_draws = 250,
             seed = 2025)
  fitd
  fite
  summary(fite)
  coef(fite)
  head(predict(fite))
  
  
  ## Bootstrapping - sequential
  plan(sequential)
  fitf = 
    np_glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
             data = test_data,
             family = poisson(),
             loss = "gls",
             n_draws = 50,
             seed = 2025)
  fite
  fitf
  summary(fitf)
  coef(fitf)
  credint(fitf)
  head(predict(fitf))
  
  
  ## Large sample approx
  fitg = 
    np_glm_b(outcome ~ x1 + x2 + x3,
             data = test_data,
             family = poisson(),
             loss = "gls",
             seed = 2025)
  fitf
  fitg
  summary(fitg)
  credint(fitg)
  coef(fitg)
  head(predict(fitg))
  
  
  
  
  # Check if custom loss function will work.
  fitg = 
    np_glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
             data = test_data,
             family = poisson(),
             loss = function(y,mu) (y-mu)^2,
             n_draws = 50,
             seed = 2025)
  fitg
  summary(fitg)
  coef(fitg)
  head(predict(fitg))
  fith = 
    np_glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
             data = test_data,
             family = poisson(),
             loss = function(y,mu) (y-mu)^2,
             seed = 2025)
  fith
  summary(fith)
  coef(fith)
  head(predict(fith))
  
  
  # Check if splines and factors work
  library(splines)
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5],
               time = rexp(N))
  test_data$outcome = 
    rpois(N,exp(-2 + test_data$x1 + 0.25 * test_data$x1^2 + 2 * (test_data$x3 %in% c("d","e"))) * test_data$time)
  fita = 
    np_glm_b(outcome ~ ns(x1,df = 5) + x2 + x3,
             data = test_data,
             n_draws = 50,
             family = poisson())
  fita
  summary(fita)
  predict(fita,
          newdata = 
            test_data[1,])
  plot(fita,type=c("pdp","ci","pi"))
  plot(fita)
  
  fitb = 
    np_glm_b(outcome ~ ns(x1,df = 5) + x2 + x3,
             data = test_data,
             family = poisson())
  fitb
  summary(fitb)
  predict(fitb,
          newdata = 
            test_data[1,])
  plot(fitb,type=c("pdp","ci","pi"))
  
  
  
  
  rm(list=ls())
  
  
  
  # Test loss-likelihood approach (nbinom) ---------------------------------
  
  # Create data
  
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5],
               time = rexp(N))
  test_data$outcome = 
    rnbinom(N,
            mu = exp(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e"))) * test_data$time,
            size = 0.7)
  
  
  # Bootstrapping approach - sequential
  plan(sequential)
  ## Make sure CI_level works (and print.lm_b works)
  fita = 
    np_glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
             data = test_data,
             family = negbinom(),
             n_draws = 50,
             seed = 2025,
             mc_error = 0.1)
  fitb = 
    np_glm_b(outcome ~ x1 + x2 + x3,
             data = test_data,
             family = negbinom(),
             n_draws = 50,
             ask_before_full_sampling = FALSE,
             seed = 2025,
             CI_level = 0.8,
             mc_error = 0.1)
  fitc = 
    np_glm_b(outcome ~ x1 + x2 + x3,
             data = test_data,
             family = negbinom(),
             n_draws = 50,
             ask_before_full_sampling = FALSE,
             seed = 2025,
             CI_level = 0.95,
             mc_error = 0.1)
  fita
  fitb
  
  ## Make sure summary.np_glm_b works
  summary(fita)
  summary(fita,
          CI_level = 0.8)
  
  ## Make sure coef.lm_b works
  coef(fita)
  
  ## Make sure credint works
  credint(fita)
  
  ## Make sure vcov works
  vcov(fita)
  
  ## Make sure prediction function works
  head(predict(fita))
  predict(fita,
          newdata = fita$data[1,])
  
  ## Test number of inputs
  fitd = 
    np_glm_b(test_data$outcome ~ test_data$x1 + offset(log(test_data$time)),
             family = negbinom(),
             n_draws = 50,
             mc_error = 0.1)
  fitd
  fite = 
    np_glm_b(test_data$outcome ~ 1 + offset(log(test_data$time)),
             family = negbinom(),
             n_draws = 50,
             mc_error = 0.1)
  fite
  fitf = 
    np_glm_b(outcome ~ x1 + offset(log(time)),
             data = test_data,
             family = negbinom(),
             n_draws = 50,
             mc_error = 0.1)
  fitf
  fitg = 
    np_glm_b(outcome ~ 1,
             data = test_data,
             family = negbinom(),
             n_draws = 50,
             mc_error = 0.1)
  fitg
  
  
  
  
  
  
  # Bootstrapping approach - parallelized
  plan(multisession, workers = 10)
  ## Make sure CI_level works (and print.lm_b works)
  fitc = 
    np_glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
             data = test_data,
             family = negbinom(),
             n_draws = 50,
             seed = 2025,
             mc_error = 0.1)
  plan(sequential)
  fita
  fitc
  summary(fitc)
  coef(fitc)
  head(predict(fitc))
  
  
  
  # large sample approach
  ## Make sure CI_level works (and print.lm_b works)
  fitd = 
    np_glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
             data = test_data,
             family = negbinom(),
             seed = 2025)
  fite = 
    np_glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
             data = test_data,
             family = negbinom(),
             seed = 2025,
             CI_level = 0.8)
  fita
  fitd
  fite
  
  ## Make sure summary.np_glm_b works
  summary(fitd)
  summary(fitd,
          CI_level = 0.8)
  
  ## Make sure coef.lm_b works
  coef(fitd)
  
  ## Make sure prediction function works
  head(predict(fitd))
  predict(fitd,
          newdata = fitd$data[1,])
  
  ## Make sure credint works
  credint(fitd)
  
  ## Make sure vcov works
  vcov(fitd)
  
  
  
  
  # Check if custom loss function will work.
  fitg = 
    np_glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
             data = test_data,
             family = negbinom(),
             loss = function(y,mu,phi) 0.5 * (y-mu)^2 / phi + 0.5*log(phi),
             n_draws = 100,
             seed = 2025,
             mc_error = 0.1)
  fith = 
    np_glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
             data = test_data,
             family = negbinom(),
             loss = function(y,mu,phi) 0.5 * (y-mu)^2 / phi + 0.5*log(phi),
             seed = 2025,
             mc_error = 0.1)
  fitg
  summary(fitg)
  coef(fitg)
  head(predict(fitg))
  
  
  # Check if splines and factors work
  library(splines)
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5],
               time = rexp(N))
  test_data$outcome = 
    rnbinom(N,
            mu = exp(-2 + test_data$x1 + 0.25 * test_data$x1^2 + 2 * (test_data$x3 %in% c("d","e"))) * test_data$time,
            size = 0.7)
  
  fita = 
    np_glm_b(outcome ~ ns(x1,df = 5) + x2 + x3,
             data = test_data,
             n_draws = 50,
             family = negbinom(),
             mc_error = 0.1)
  fita
  summary(fita)
  predict(fita,
          newdata = 
            test_data[1,])
  plot(fita,type=c("ci","pi"))
  plot(fita)
  
  fitb = 
    np_glm_b(outcome ~ ns(x1,df = 5) + x2 + x3,
             data = test_data,
             family = negbinom())
  fitb
  summary(fitb)
  predict(fitb,
          newdata = 
            test_data[1,])
  plot(fitb)
  
  
  
  
  rm(list=ls())
  
  
  
  
  
  # Test prop_test_b --------------------------------------------------------
  
  prop_test_b(14,
              19)
  
  prop_test_b(14,
              n_total = 19)
  
  prop_test_b(14,
              19,
              p = 0.45)
  
  
  
  prop_test_b(c(14,22),
              c(19,45))
  
  prop_test_b(c(99,1),
              n_total = c(100,50))
  
  
  prop_test_b(c(14,22),
              c(19,45),
              ROPE = c(1.0 / 1.05, 1.05))
  
  
  
  # Test t_test_b -----------------------------------------------------------
  
  t_test_b(rnorm(50))
  t_test_b(outcome ~ 1,
           data = data.frame(outcome = rnorm(50)))
  t_test_b(rnorm(50),
           rnorm(15,1))
  t_test_b(outcome ~ asdf,
           data = 
             data.frame(outcome = c(rnorm(50),
                                    rnorm(15,1)),
                        asdf = rep(c("a","b"),c(50,15))))
  
  try({
    t_test_b(rnorm(50),
             rnorm(15,1),
             paired = TRUE)
    print("If seeing this message, there's a problem!")
  },silent = T)
  t_test_b(rnorm(50),
           rnorm(50,1),
           paired = TRUE)
  
  
  # Test poisson_test_b -----------------------------------------------------
  
  poisson_test_b(x = 12)
  poisson_test_b(x = 12,
                 r = 11)
  poisson_test_b(x = 12,
                 offset = 2,
                 r = 11)
  
  poisson_test_b(x = 12,
                 offset = 2,
                 r = 11,
                 prior = "flat")
  
  poisson_test_b(x = 12,
                 offset = 2,
                 r = 11,
                 prior_shape_rate = c(1,1))
  
  poisson_test_b(x = c(12,20),
                 offset = c(10,9),
                 prior = "flat")
  
  poisson_test_b(x = c(12,20),
                 offset = c(10,9),
                 prior = "jeff")
  
  poisson_test_b(x = c(12,20),
                 offset = c(10,9),
                 prior = "jeff")
  
  
  
  # sign_test_b -------------------------------------------------------------
  
  # Test input
  sign_test_b(x = rnorm(50))
  sign_test_b(x = rnorm(50,1))
  sign_test_b(x = rbinom(50,1,0.6))
  sign_test_b(x = rnorm(50,1),
              y = rnorm(50,0))
  
  # Test prior
  sign_test_b(x = rnorm(50,0),
              y = rnorm(50,1),
              prior = "uniform")
  sign_test_b(x = rnorm(50,0),
              y = rnorm(50,1),
              prior_shapes = c(2,2))
  
  # Test ROPE
  sign_test_b(x = rnorm(50,1),
              ROPE = 0.1)
  sign_test_b(x = rnorm(50,1),
              ROPE = c(0.45,0.7))
  
  # Test changing p0
  sign_test_b(x = rnorm(50,1),
              p0 = 0.7)
  sign_test_b(x = rnorm(50,1),
              p0 = 0.7,
              ROPE = 0.1)
  
  try({
    sign_test_b(x = rnorm(50,1),
                p0 = 0.7,
                ROPE = 0.3)
    print("If seeing this, you have a problem!")
  },silent=T)
  
  
  
  # Case-control study ------------------------------------------------------
  
  # Test input
  case_control_b(matrix(c(8,47,1,26),2,2))
  case_control_b(c(8,47),
                 c(1,26))
  case_control_b(x = matrix(c(8,47,1,26),2,2))
  
  # Test large sample
  case_control_b(x = 5 + matrix(c(8,47,1,26),2,2))
  
  # Test ROPE
  case_control_b(x = matrix(c(8,47,1,26),2,2),
                 ROPE = 1.05)
  case_control_b(x = 5 + matrix(c(8,47,1,26),2,2),
                 ROPE = 1.05)
  
  # Test prior
  case_control_b(x = 5 + matrix(c(8,47,1,26),2,2),
                 ROPE = 1.05)
  case_control_b(x = 5 + matrix(c(8,47,1,26),2,2),
                 ROPE = 1.05,
                 prior_mean = 10)
  case_control_b(x = 5 + matrix(c(8,47,1,26),2,2),
                 ROPE = 1.05,
                 prior_sd = 0.01)
  
  
  
  # Wilcoxon signed rank test -----------------------------------------------
  
  N = 15
  set.seed(2025)
  
  # Small samples
  ## Test if analysis is done correctly
  test_data_small = 
    data.frame(x = rbeta(N,2,10),
               y = rbeta(N,5,10))
  test_data_small_null = 
    data.frame(x = rbeta(N,2,10),
               y = rbeta(N,2,10))
  hist(test_data_small$x - test_data_small$y)
  wilcoxon_test_b(test_data_small$x - test_data_small$y)
  hist(test_data_small_null$x - test_data_small_null$y)
  wilcoxon_test_b(test_data_small_null$x - test_data_small_null$y)
  
  ## Test input
  wilcoxon_test_b(test_data_small$x,
                  test_data_small$y,
                  paired = TRUE)
  
  ## Test priors
  wilcoxon_test_b(test_data_small$x - test_data_small$y,
                  prior = "centered")
  wilcoxon_test_b(test_data_small$x - test_data_small$y,
                  prior = "uniform")
  wilcoxon_test_b(test_data_small$x - test_data_small$y,
                  prior_shapes = c(5,5))
  
  ## Test ROPE
  wilcoxon_test_b(test_data_small$x - test_data_small$y,
                  ROPE = 0.1)
  wilcoxon_test_b(test_data_small$x - test_data_small$y,
                  ROPE = c(0.4,0.65))
  
  ## Test returned object
  test = 
    wilcoxon_test_b(test_data_small$x - test_data_small$y)
  str(test[-7])
  test$prob_plot
  
  
  # Large samples
  N = 150
  set.seed(2025)
  test_data_big = 
    data.frame(x = rbeta(N,2,10),
               y = rbeta(N,5,10))
  test_data_big_null = 
    data.frame(x = rbeta(N,5,10),
               y = rbeta(N,5,10))
  
  ## Test if analysis is done correctly
  hist(test_data_big$x - test_data_big$y)
  wilcoxon_test_b(test_data_big$x - test_data_big$y)
  hist(test_data_big_null$x - test_data_big_null$y)
  wilcoxon_test_b(test_data_big_null$x - test_data_big_null$y)
  
  ## Test input
  wilcoxon_test_b(test_data_big$x,
                  test_data_big$y,
                  paired = TRUE)
  
  ## Test priors
  wilcoxon_test_b(test_data_big$x - test_data_big$y,
                  prior = "cent")
  wilcoxon_test_b(test_data_big$x - test_data_big$y,
                  prior = "uniform")
  wilcoxon_test_b(test_data_big$x - test_data_big$y,
                  prior_shapes = c(5,5))
  
  ## Test ROPE
  wilcoxon_test_b(test_data_big$x - test_data_big$y,
                  ROPE = 0.1)
  wilcoxon_test_b(test_data_big$x - test_data_big$y,
                  ROPE = c(0.4,0.65))
  
  ## Test returned object
  test = 
    wilcoxon_test_b(test_data_big$x - test_data_big$y)
  str(test[-6])
  test$prob_plot
  
  
  # Wilcoxon rank sum test --------------------------------------------------
  
  
  # Small samples
  N = 15
  x = rbeta(N,2,10)
  y = rbeta(N + 1,5,10)
  y_null = rbeta(N+1,2,10)
  curve(dbeta(x,2,10)); curve(dbeta(x,5,10),add= T,lty = 2,col=2)
  
  ## Test if analysis is done correctly
  wilcoxon_test_b(x,y)
  wilcoxon_test_b(x,y_null)
  
  ## Test priors
  wilcoxon_test_b(x,
                  y,
                  prior = "cent")
  wilcoxon_test_b(x,
                  y,
                  prior = "uniform")
  wilcoxon_test_b(x,
                  y,
                  prior_shapes = c(5,5))
  
  ## Test ROPE
  wilcoxon_test_b(x,
                  y,
                  ROPE = 0.1)
  wilcoxon_test_b(x,
                  y,
                  ROPE = c(0.1,0.8))
  
  
  
  # Large samples
  N = 150
  x = rbeta(N,2,10)
  y = rbeta(N + 1,5,10)
  y_null = rbeta(N+1,2,10)
  curve(dbeta(x,2,10)); curve(dbeta(x,5,10),add= T,lty = 2,col=2)
  
  ## Test if analysis is done correctly
  wilcoxon_test_b(x,y)
  wilcoxon_test_b(x,y_null)
  
  ## Test priors
  wilcoxon_test_b(x,
                  y,
                  prior = "cent")
  wilcoxon_test_b(x,
                  y,
                  prior = "uniform")
  wilcoxon_test_b(x,
                  y,
                  prior_shapes = c(5,5))
  
  ## Test ROPE
  wilcoxon_test_b(x,
                  y,
                  ROPE = 0.1)
  wilcoxon_test_b(x,
                  y,
                  ROPE = c(0.1,0.8))
  
  
  
  # cor_test_b --------------------------------------------------------------
  
  set.seed(2025)
  N = 50
  x = rnorm(N)
  y = x + 4 * rnorm(N)
  
  # Check input
  cor_test_b(x,y)
  cor_test_b(~ asdf + qwer,
             data = data.frame(asdf = x,
                               qwer = y))
  
  # Check tau
  cor_test_b(x,y,
             tau = 0.04)
  
  # Check ROPE
  cor_test_b(x,y,
             ROPE = 0.1)
  cor_test_b(x,y,
             ROPE = c(-0.9,0.02))
  
  # Check prior
  cor_test_b(x,y,
             prior = "uniform")
  cor_test_b(x,y,
             prior = "negati")
  cor_test_b(x,y,
             prior = "pos")
  cor_test_b(x,y,
             prior_shapes = c(10,10))
  
  
  
  
  
  
  # chisq.test --------------------------------------------------------------
  
  set.seed(2025)
  nR = 5
  nC = 3
  dep_probs = 
    extraDistr::rdirichlet(1,rep(2,nR*nC)) |> 
    matrix(nR,nC)
  ind_probs = 
    tcrossprod(rowSums(dep_probs),
               colSums(dep_probs))
  
  # Multinomial sampling design
  ## Test with big N
  N = 500
  ### Check if answers are reasonable wrt the truth
  fita = 
    independence_b(round(N * dep_probs))
  independence_b(round(N * ind_probs))
  str(fita)
  
  ### Try other priors
  independence_b(round(N * dep_probs),
                 prior = "uniform")
  independence_b(round(N * dep_probs),
                 prior_shapes = 2)
  independence_b(round(N * dep_probs),
                 prior_shapes = matrix(1:(nR*nC),nR,nC))
  independence_b(round(N * dep_probs),
                 prior_shapes = rep(2,nR*nC))
  
  ## Test with small N ... not very good until N = 150
  N = 50
  independence_b(round(N * dep_probs))
  independence_b(round(N * ind_probs))
  N = 100
  independence_b(round(N * dep_probs))
  independence_b(round(N * ind_probs))
  N = 150
  independence_b(round(N * dep_probs))
  independence_b(round(N * ind_probs))
  
  
  # Fixed rows sampling design
  ## Test with big N
  N = 500
  ### Check if answers are reasonable wrt the truth
  fita = 
    independence_b(round(N * dep_probs),
                   sampling_design = "rows")
  independence_b(round(N * ind_probs),
                 sampling_design = "rows")
  str(fita)
  
  ### Try other priors
  independence_b(round(N * dep_probs),
                 prior = "uniform",
                 sampling_design = "rows")
  independence_b(round(N * dep_probs),
                 prior_shapes = 2,
                 sampling_design = "rows")
  independence_b(round(N * dep_probs),
                 prior_shapes = matrix(1:(nR*nC),nR,nC),
                 sampling_design = "rows")
  independence_b(round(N * dep_probs),
                 prior_shapes = rep(2,nR*nC),
                 sampling_design = "rows")
  
  ## Test with small N
  N = 50
  independence_b(round(N * dep_probs),
                 sampling_design = "rows")
  independence_b(round(N * ind_probs),
                 sampling_design = "rows")
  N = 100
  independence_b(round(N * dep_probs),
                 sampling_design = "rows")
  independence_b(round(N * ind_probs),
                 sampling_design = "rows")
  N = 150
  independence_b(round(N * dep_probs),
                 sampling_design = "rows")
  independence_b(round(N * ind_probs),
                 sampling_design = "rows")
  
  
  
  # Fixed columns sampling design
  ## Test with big N
  N = 500
  ### Check if answers are reasonable wrt the truth
  fita = 
    independence_b(round(N * dep_probs),
                   sampling_design = "fixed columns")
  independence_b(round(N * ind_probs),
                 sampling_design = "cols")
  str(fita)
  
  ### Try other priors
  independence_b(round(N * dep_probs),
                 prior = "uniform",
                 sampling_design = "cols")
  independence_b(round(N * dep_probs),
                 prior_shapes = 2,
                 sampling_design = "cols")
  independence_b(round(N * dep_probs),
                 prior_shapes = matrix(1:(nR*nC),nR,nC),
                 sampling_design = "cols")
  independence_b(round(N * dep_probs),
                 prior_shapes = rep(2,nR*nC),
                 sampling_design = "cols")
  
  ## Test with small N
  N = 50
  independence_b(round(N * dep_probs),
                 sampling_design = "cols")
  independence_b(round(N * ind_probs),
                 sampling_design = "cols")
  N = 100
  independence_b(round(N * dep_probs),
                 sampling_design = "cols")
  independence_b(round(N * ind_probs),
                 sampling_design = "cols")
  N = 150
  independence_b(round(N * dep_probs),
                 sampling_design = "cols")
  independence_b(round(N * ind_probs),
                 sampling_design = "cols")
  
  # Test to make sure rows and columns on t(x) are equivalent
  test1 = 
    independence_b(round(N * dep_probs),
                   sampling_design = "rows")
  test2 = 
    independence_b(round(N * dep_probs) |> t(),
                   sampling_design = "cols")
  all(near(test1$posterior_shapes,
           t(test2$posterior_shapes)))
  all(near(test1$posterior_mean,
           t(test2$posterior_mean)))
  all(near(test1$lower_bound,
           t(test2$lower_bound)))
  all(near(test1$individual_ROPE,
           t(test2$individual_ROPE)))
  all(near(test1$prob_p_j_given_i_less_than_p_j,
           t(test2$prob_p_j_given_i_less_than_p_j)))
  

  # Semi-parametric survival ------------------------------------------------
  
  
  set.seed(2025)
  N = 300
  test_data = 
    data.frame(outcome = 
                 # rexp(N,1/3))
                 rweibull(N,2,5))
  hist(test_data$outcome)
  test_data$observed = 
    ifelse(test_data$outcome >= 7, 0, 1)
  test_data$outcome =
    ifelse(dplyr::near(test_data$observed,1), test_data$outcome, 7)
  
  Surv_object = survival::Surv(test_data$outcome,
                               test_data$observed)
  
  
  
  set.seed(2025)
  N = 300
  test_data = 
    data.frame(outcome = 
                 c(rweibull(2*N/3,2,5),
                   rweibull(N/3,2,10)),
               x1 = rep(letters[1:3],each = N/3))
  hist(test_data$outcome[which(test_data$x1 == "a")])
  hist(test_data$outcome[which(test_data$x1 == "b")])
  hist(test_data$outcome[which(test_data$x1 == "c")])
  test_data$observed = 
    ifelse(test_data$outcome >= 15, 0, 1)
  test_data$outcome =
    ifelse(dplyr::near(test_data$observed,1), test_data$outcome, 15)
  
  Surv_object = survival::Surv(test_data$outcome,
                               test_data$observed)
  
  
  
  test_data = 
    data.frame(outcome = 
                 c(rexp(2*N/3,2),
                   rexp(N/3,2,10)),
               x1 = rep(letters[1:3],each = N/3))
  hist(test_data$outcome[which(test_data$x1 == "a")])
  hist(test_data$outcome[which(test_data$x1 == "b")])
  hist(test_data$outcome[which(test_data$x1 == "c")])
  test_data$observed = 
    ifelse(test_data$outcome >= 15, 0, 1)
  test_data$outcome =
    ifelse(dplyr::near(test_data$observed,1), test_data$outcome, 15)
  
  Surv_object = survival::Surv(test_data$outcome,
                               test_data$observed)
  
  
  # # LME ---------------------------------------------------------------------
  # 
  # data("sleepstudy", package = "lme4")
  # 
  # library(glmmTMB)
  # set.seed(2025)
  # N = 500
  # test_data = 
  #   data.frame(x1 = rnorm(N),
  #              x2 = rnorm(N),
  #              x3 = letters[1:5])
  # test_re = 
  #   data.frame(x3 = letters[1:5],
  #              z = rnorm(5))
  # test_data =
  #   dplyr::left_join(
  #     test_data,
  #     test_re,
  #     by = "x3")
  # test_data$outcome = 
  #   rbinom(N,1,1.0 / (1.0 + exp(-(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) + test_data$z ))))
  # 
  # fita = glmmTMB(outcome ~ x1 + (1|x3), 
  #                data = test_data, 
  #                family = binomial)
  # 
  # # TMB objective
  # obj <- fita$obj
  # 
  # # Current parameter vector (internal TMB ordering and parameterization)
  # par_hat1 <- obj$par
  # par_hat1[1:2] = unlist(fixef(fita))
  # par_hat1[3] = getME(fita,"theta")#log(sqrt(summary(fita)$varcor$cond$x3[1,1]))
  # par_hat2 <- par_hat1 + rnorm(3,sd = 0.1)
  # 
  # # Evaluate objective (usually: negative log-likelihood up to a constant)
  # ## obj$fn is -log likelihood()
  # 2 * obj$fn(par_hat1)
  # fita
  # gr_hat  <- obj$gr(par_hat)
  # he_hat  <- obj$he(par_hat)
  # 
  # # Evaluate at user-specified parameter vector
  # par_new <- par_hat
  # par_new[1] <- par_new[1] + 0.2
  # 
  # val_new <- obj$fn(par_new)
  # gr_new  <- obj$gr(par_new)
  # 
}




