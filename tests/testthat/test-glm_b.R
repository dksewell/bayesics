

# Binomial ----------------------------------------------------------------

test_that("Test glm_b for binomial data fitting with VB",{
  
  # Generate some data
  set.seed(2025)
  N = 100
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5])
  test_data$outcome = 
    rbinom(N,1,1.0 / (1.0 + exp(-(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) ))))
  
  
  # Test VB fit
  expect_no_error(
    fita <-
      glm_b(outcome ~ x1 + x2 + x3,
            data = test_data,
            family = binomial(),
            seed = 2025)
  )
  
  # Default fit ought to be VB
  expect_identical(fita$algorithm,"VB")
  
  # Make sure print works
  expect_no_error(fita)
  
  # Make sure coef works
  expect_type(coef(fita),"double")
  
  # Make sure credint works
  expect_true(is.matrix(credint(fita)))
  
  # Make sure vcov works
  expect_true(is.matrix(vcov(fita)))
  
  # Make sure summary works
  expect_no_error(
    s <- 
      summary(fita)
  )
  ## Check output format
  expect_s3_class(s,c("tbl_df", "tbl", "data.frame"))
  
  expect_identical(colnames(s),
                   c("Variable","Post Mean","Lower","Upper","Prob Dir",
                     "ROPE","ROPE bounds","BF favoring alternative",
                     "Interpretation"))
  expect_type(s$Variable,"character")
  expect_type(s$`Post Mean`,"double")
  expect_type(s$Lower,"double")
  expect_type(s$Upper,"double")
  expect_type(s$`Prob Dir`,"double")
  expect_type(s$ROPE,"double")
  expect_type(s$`ROPE bounds`,"character")
  expect_type(s$`BF favoring alternative`,"double")
  expect_type(s$Interpretation,"character")
  
  # Make sure prediction function works
  expect_no_error(predict(fita))
  expect_no_error(predict(fita,
                          newdata = fita$data[1,]))
  expect_gte(predict(fita,
                     newdata = fita$data[1,],
                     CI_level = 0.8)$CI_lower[1],
             predict(fita,
                     newdata = fita$data[1,],
                     CI_level = 0.9)$CI_lower[1])
  expect_gte(predict(fita,
                    newdata = fita$data[1,],
                    PI_level = 0.8)$PI_lower[1],
            predict(fita,
                    newdata = fita$data[1,],
                    PI_level = 0.9)$PI_lower[1])
  
  # Make sure savage-dickey ratio works
  expect_s3_class(bayes_factors(fita),
                  c("tbl_df", "tbl", "data.frame"))
  expect_s3_class(bayes_factors(fita,by = "v"),
                  c("tbl_df", "tbl", "data.frame"))
  
  # Make sure information criteria work
  expect_type(AIC(fita),"double")
  expect_type(BIC(fita),"double")
  expect_type(DIC(fita),"double")
  expect_type(WAIC(fita),"double")
  
  # Test number of inputs
  expect_no_error(
    glm_b(test_data$outcome ~ test_data$x1,
          family = binomial(),
          prior = "normal")
  )
  expect_no_error(
    glm_b(test_data$outcome ~ 1,
          family = binomial(),
          prior = "normal")
  )
  expect_no_error(
    glm_b(outcome ~ x1,
          data = test_data,
          family = binomial(),
          prior = "normal")
  )
  expect_no_error(
    glm_b(outcome ~ 1,
          data = test_data,
          family = binomial(),
          prior = "normal")
  )
  
  # Test different priors
  expect_no_error(
    glm_b(outcome ~ x1 + x2 + x3,
          data = test_data,
          family = binomial(),
          prior = "normal",
          seed = 2025)
  )
  expect_no_error(
    fitb <-
      glm_b(outcome ~ x1 + x2 + x3,
            data = test_data,
            family = binomial(),
            prior = "improper")
  )
  ## Make sure bayes_factors doesn't work for improper prior
  expect_error(bayes_factors(fitb))
  
  
  
  # Test plot
  expect_s3_class(plot(fita,
                       type = "diagnostics"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_s3_class(plot(fita,
                       type = "pdp"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       variable = "x1",
                       combine_pi_ci = TRUE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       variable = "x1",
                       combine_pi_ci = FALSE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       combine_pi_ci = TRUE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       combine_pi_ci = TRUE),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       combine_pi_ci = FALSE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_s3_class(plot(fita,
                       type = "cr",
                       variable = "x1"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "cr"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  
  # # Check parallelization
  # plan(multisession,workers = 5)
  # expect_no_error(
  #   glm_b(outcome ~ x1 + x2 + x3,
  #         data = test_data,
  #         family = binomial(),
  #         seed = 2025)
  # )
  # plan(sequential)
  
  
})


test_that("Test glm_b for binomial data fitting with IS",{
  
  # Generate some data
  set.seed(2025)
  N = 100
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5])
  test_data$outcome = 
    rbinom(N,1,1.0 / (1.0 + exp(-(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) ))))
  
  
  # Test IS fit
  expect_no_error(
    fita <-
      glm_b(outcome ~ x1 + x2 + x3,
            data = test_data,
            family = binomial(),
            seed = 2025,
            algorithm = "IS")
  )
  
  # Default fit ought to be VB
  expect_identical(fita$algorithm,"IS")
  
  # Make sure print works
  expect_no_error(fita)
  
  # Make sure coef works
  expect_type(coef(fita),"double")
  
  # Make sure credint works
  expect_true(is.matrix(credint(fita)))
  
  # Make sure vcov works
  expect_true(is.matrix(vcov(fita)))
  
  # Make sure summary works
  expect_no_error(
    s <- 
      summary(fita)
  )
  ## Check output format
  expect_s3_class(s,c("tbl_df", "tbl", "data.frame"))
  
  expect_identical(colnames(s),
                   c("Variable","Post Mean","Lower","Upper","Prob Dir",
                     "ROPE","ROPE bounds","BF favoring alternative",
                     "Interpretation"))
  expect_type(s$Variable,"character")
  expect_type(s$`Post Mean`,"double")
  expect_type(s$Lower,"double")
  expect_type(s$Upper,"double")
  expect_type(s$`Prob Dir`,"double")
  expect_type(s$ROPE,"double")
  expect_type(s$`ROPE bounds`,"character")
  expect_type(s$`BF favoring alternative`,"double")
  expect_type(s$Interpretation,"character")
  
  # Make sure prediction function works
  expect_no_error(predict(fita))
  expect_no_error(predict(fita,
                          newdata = fita$data[1,]))
  expect_gte(predict(fita,
                     newdata = fita$data[1,],
                     CI_level = 0.8)$CI_lower[1],
             predict(fita,
                     newdata = fita$data[1,],
                     CI_level = 0.9)$CI_lower[1])
  expect_gte(predict(fita,
                     newdata = fita$data[1,],
                     PI_level = 0.8)$PI_lower[1],
             predict(fita,
                     newdata = fita$data[1,],
                     PI_level = 0.9)$PI_lower[1])
  
  # Make sure savage-dickey ratio works
  expect_s3_class(bayes_factors(fita),
                  c("tbl_df", "tbl", "data.frame"))
  expect_s3_class(bayes_factors(fita,by = "v"),
                  c("tbl_df", "tbl", "data.frame"))
  
  # Make sure information criteria work
  expect_type(AIC(fita),"double")
  expect_type(BIC(fita),"double")
  expect_type(DIC(fita),"double")
  expect_type(WAIC(fita),"double")
  
  # Test number of inputs
  expect_no_error(
    glm_b(test_data$outcome ~ test_data$x1,
          family = binomial(),
          prior = "normal",
          algorithm = "IS")
  )
  expect_no_error(
    glm_b(test_data$outcome ~ 1,
          family = binomial(),
          prior = "normal",
          algorithm = "IS")
  )
  expect_no_error(
    glm_b(outcome ~ x1,
          data = test_data,
          family = binomial(),
          prior = "normal",
          algorithm = "IS")
  )
  expect_no_error(
    glm_b(outcome ~ 1,
          data = test_data,
          family = binomial(),
          prior = "normal",
          algorithm = "IS")
  )
  
  # Test different priors
  expect_no_error(
    glm_b(outcome ~ x1 + x2 + x3,
          data = test_data,
          family = binomial(),
          prior = "normal",
          seed = 2025,
          algorithm = "IS")
  )
  expect_no_error(
    fitb <-
      glm_b(outcome ~ x1 + x2 + x3,
            data = test_data,
            family = binomial(),
            prior = "improper",
            algorithm = "IS")
  )
  ## Make sure bayes_factors doesn't work for improper prior
  expect_error(bayes_factors(fitb))
  
  
  
  # Test mc_error
  expect_no_error(
    fitc <-
      glm_b(outcome ~ x1 + x2 + x3,
            data = test_data,
            family = binomial(),
            seed = 2025,
            algorithm = "IS",
            mc_error = 0.1)
  )
  expect_no_error(
    fitd <-
      glm_b(outcome ~ x1 + x2 + x3,
            data = test_data,
            family = binomial(),
            seed = 2025,
            algorithm = "IS",
            mc_error = 0.005)
  )
  expect_lt(nrow(fitc$proposal_draws),
            nrow(fitd$proposal_draws))
  
  
  
  # Test plot
  expect_s3_class(plot(fita,
                       type = "diagnostics"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_s3_class(plot(fita,
                       type = "pdp"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       variable = "x1",
                       combine_pi_ci = TRUE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       variable = "x1",
                       combine_pi_ci = FALSE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       combine_pi_ci = TRUE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       combine_pi_ci = TRUE),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       combine_pi_ci = FALSE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "cr",
                       variable = "x1"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "cr"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  
  # # Check parallelization
  # plan(multisession,workers = 5)
  # expect_no_error(
  #   glm_b(outcome ~ x1 + x2 + x3,
  #         data = test_data,
  #         family = binomial(),
  #         seed = 2025,
  #         algorithm = "IS")
  # )
  # plan(sequential)
  
  
})


test_that("Test glm_b for binomial data fitting with LSA",{
  
  # Generate some data
  set.seed(2025)
  N = 100
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5])
  test_data$outcome = 
    rbinom(N,1,1.0 / (1.0 + exp(-(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) ))))
  
  
  # Test IS fit
  expect_no_error(
    fita <-
      glm_b(outcome ~ x1 + x2 + x3,
            data = test_data,
            family = binomial(),
            seed = 2025,
            algorithm = "LSA")
  )
  
  # Default fit ought to be VB
  expect_identical(fita$algorithm,"LSA")
  
  # Make sure print works
  expect_no_error(fita)
  
  # Make sure coef works
  expect_type(coef(fita),"double")
  
  # Make sure credint works
  expect_true(is.matrix(credint(fita)))
  
  # Make sure vcov works
  expect_true(is.matrix(vcov(fita)))
  
  # Make sure summary works
  expect_no_error(
    s <- 
      summary(fita)
  )
  ## Check output format
  expect_s3_class(s,c("tbl_df", "tbl", "data.frame"))
  
  expect_identical(colnames(s),
                   c("Variable","Post Mean","Lower","Upper","Prob Dir",
                     "ROPE","ROPE bounds","BF favoring alternative",
                     "Interpretation"))
  expect_type(s$Variable,"character")
  expect_type(s$`Post Mean`,"double")
  expect_type(s$Lower,"double")
  expect_type(s$Upper,"double")
  expect_type(s$`Prob Dir`,"double")
  expect_type(s$ROPE,"double")
  expect_type(s$`ROPE bounds`,"character")
  expect_type(s$`BF favoring alternative`,"double")
  expect_type(s$Interpretation,"character")
  
  # Make sure prediction function works
  expect_no_error(predict(fita))
  expect_no_error(predict(fita,
                          newdata = fita$data[1,]))
  expect_gte(predict(fita,
                     newdata = fita$data[1,],
                     CI_level = 0.8)$CI_lower[1],
             predict(fita,
                     newdata = fita$data[1,],
                     CI_level = 0.9)$CI_lower[1])
  expect_gte(predict(fita,
                     newdata = fita$data[1,],
                     PI_level = 0.8)$PI_lower[1],
             predict(fita,
                     newdata = fita$data[1,],
                     PI_level = 0.9)$PI_lower[1])
  
  # Make sure savage-dickey ratio works
  expect_s3_class(bayes_factors(fita),
                  c("tbl_df", "tbl", "data.frame"))
  expect_s3_class(bayes_factors(fita,by = "v"),
                  c("tbl_df", "tbl", "data.frame"))
  
  # Make sure information criteria work
  expect_type(AIC(fita),"double")
  expect_type(BIC(fita),"double")
  expect_type(DIC(fita),"double")
  expect_type(WAIC(fita),"double")
  
  # Test number of inputs
  expect_no_error(
    glm_b(test_data$outcome ~ test_data$x1,
          family = binomial(),
          prior = "normal",
          algorithm = "LSA")
  )
  expect_no_error(
    glm_b(test_data$outcome ~ 1,
          family = binomial(),
          prior = "normal",
          algorithm = "LSA")
  )
  expect_no_error(
    glm_b(outcome ~ x1,
          data = test_data,
          family = binomial(),
          prior = "normal",
          algorithm = "LSA")
  )
  expect_no_error(
    glm_b(outcome ~ 1,
          data = test_data,
          family = binomial(),
          prior = "normal",
          algorithm = "LSA")
  )
  
  # Test different priors
  expect_no_error(
    glm_b(outcome ~ x1 + x2 + x3,
          data = test_data,
          family = binomial(),
          prior = "normal",
          seed = 2025,
          algorithm = "LSA")
  )
  expect_no_error(
    fitb <-
      glm_b(outcome ~ x1 + x2 + x3,
            data = test_data,
            family = binomial(),
            prior = "improper",
            algorithm = "LSA")
  )
  ## Make sure bayes_factors doesn't work for improper prior
  expect_error(bayes_factors(fitb))
  
  
  
  # Test plot
  expect_s3_class(plot(fita,
                       type = "diagnostics"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_s3_class(plot(fita,
                       type = "pdp"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       variable = "x1",
                       combine_pi_ci = TRUE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       variable = "x1",
                       combine_pi_ci = FALSE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       combine_pi_ci = TRUE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       combine_pi_ci = TRUE),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       combine_pi_ci = FALSE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "cr",
                       variable = "x1"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "cr"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  
  # # Check parallelization
  # plan(multisession,workers = 5)
  # expect_no_error(
  #   glm_b(outcome ~ x1 + x2 + x3,
  #         data = test_data,
  #         family = binomial(),
  #         seed = 2025,
  #         algorithm = "LSA")
  # )
  # plan(sequential)
  
  
})


test_that("Test glm_b for binomial data with >1 trials",{
  
  # Generate some data
  set.seed(2025)
  N = 100
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5],
               n_trials = rpois(N,20))
  test_data$outcome = 
    rbinom(N,
           test_data$n_trials,
           1.0 / (1.0 + exp(-(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) ))))
  
  
  # Test VB fit
  expect_no_error(
    fita <-
      glm_b(outcome ~ x1 + x2 + x3,
            data = test_data,
            trials = "n_trials",
            family = binomial(),
            seed = 2025)
  )
  expect_no_error(
    fitb <-
      glm_b(outcome ~ x1 + x2 + x3,
            data = test_data,
            trials = test_data$n_trials,
            family = binomial(),
            seed = 2025)
  )
  expect_equal(fita$summary,
               fitb$summary)
  
  expect_no_error(
    preds1 <-
      predict(fita,
              newdata = fita$data[1,],
              trials = "n_trials")
  )
  expect_no_error(
    preds2 <-
      predict(fita,
              newdata = fita$data[1,],
              trials = 2e3)
  )
  expect_lt(preds1$`Post Mean`,
            preds2$`Post Mean`)
  
  
  # Test IS fit
  expect_no_error(
    fita <-
      glm_b(outcome ~ x1 + x2 + x3,
            data = test_data,
            trials = "n_trials",
            family = binomial(),
            seed = 2025,
            algorithm = "IS")
  )
  expect_no_error(
    fitb <-
      glm_b(outcome ~ x1 + x2 + x3,
            data = test_data,
            trials = test_data$n_trials,
            family = binomial(),
            seed = 2025,
            algorithm = "IS")
  )
  
  expect_no_error(
    preds1 <-
      predict(fita,
              newdata = fita$data[1,],
              trials = "n_trials")
  )
  expect_no_error(
    preds2 <-
      predict(fita,
              newdata = fita$data[1,],
              trials = 2e3)
  )
  expect_lt(preds1$`Post Mean`,
            preds2$`Post Mean`)
  
  
  # Test VB fit
  expect_no_error(
    fita <-
      glm_b(outcome ~ x1 + x2 + x3,
            data = test_data,
            trials = "n_trials",
            family = binomial(),
            seed = 2025,
            algorithm = "LSA")
  )
  expect_no_error(
    fitb <-
      glm_b(outcome ~ x1 + x2 + x3,
            data = test_data,
            trials = test_data$n_trials,
            family = binomial(),
            seed = 2025,
            algorithm = "LSA")
  )
  expect_equal(fita$summary,
               fitb$summary)
  
  expect_no_error(
    preds1 <-
      predict(fita,
              newdata = fita$data[1,],
              trials = "n_trials")
  )
  expect_no_error(
    preds2 <-
      predict(fita,
              newdata = fita$data[1,],
              trials = 2e3)
  )
  expect_lt(preds1$`Post Mean`,
            preds2$`Post Mean`)
  
  
})



# Poisson -----------------------------------------------------------------


test_that("Test glm_b for poisson data fitting with VB",{
  
  # Generate some data
  set.seed(2025)
  N = 100
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5],
               time = rexp(N))
  test_data$outcome = 
    rpois(N,exp(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e"))) * test_data$time)
  
  
  # Test VB fit
  expect_no_error(
    fita <-
      glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
            data = test_data,
            family = poisson(),
            seed = 2025)
  )
  
  # Default fit ought to be VB
  expect_identical(fita$algorithm,"VB")
  
  # Make sure print works
  expect_no_error(fita)
  
  # Make sure coef works
  expect_type(coef(fita),"double")
  
  # Make sure credint works
  expect_true(is.matrix(credint(fita)))
  
  # Make sure vcov works
  expect_true(is.matrix(vcov(fita)))
  
  # Make sure summary works
  expect_no_error(
    s <- 
      summary(fita)
  )
  ## Check output format
  expect_s3_class(s,c("tbl_df", "tbl", "data.frame"))
  
  expect_identical(colnames(s),
                   c("Variable","Post Mean","Lower","Upper","Prob Dir",
                     "ROPE","ROPE bounds","BF favoring alternative",
                     "Interpretation"))
  expect_type(s$Variable,"character")
  expect_type(s$`Post Mean`,"double")
  expect_type(s$Lower,"double")
  expect_type(s$Upper,"double")
  expect_type(s$`Prob Dir`,"double")
  expect_type(s$ROPE,"double")
  expect_type(s$`ROPE bounds`,"character")
  expect_type(s$`BF favoring alternative`,"double")
  expect_type(s$Interpretation,"character")
  
  # Make sure prediction function works
  expect_no_error(predict(fita))
  expect_no_error(predict(fita,
                          newdata = fita$data[1,]))
  expect_gte(predict(fita,
                     newdata = fita$data[1,],
                     CI_level = 0.8)$CI_lower[1],
             predict(fita,
                     newdata = fita$data[1,],
                     CI_level = 0.9)$CI_lower[1])
  expect_gte(predict(fita,
                     newdata = fita$data[1,],
                     PI_level = 0.8)$PI_lower[1],
             predict(fita,
                     newdata = fita$data[1,],
                     PI_level = 0.9)$PI_lower[1])
  
  # Make sure savage-dickey ratio works
  expect_s3_class(bayes_factors(fita),
                  c("tbl_df", "tbl", "data.frame"))
  expect_s3_class(bayes_factors(fita,by = "v"),
                  c("tbl_df", "tbl", "data.frame"))
  
  # Make sure information criteria work
  expect_type(AIC(fita),"double")
  expect_type(BIC(fita),"double")
  expect_type(DIC(fita),"double")
  expect_type(WAIC(fita),"double")
  
  # Test number of inputs
  expect_no_error(
    glm_b(test_data$outcome ~ test_data$x1,
          family = poisson(),
          prior = "normal")
  )
  expect_no_error(
    glm_b(test_data$outcome ~ 1,
          family = poisson(),
          prior = "normal")
  )
  expect_no_error(
    glm_b(outcome ~ x1,
          data = test_data,
          family = poisson(),
          prior = "normal")
  )
  expect_no_error(
    glm_b(outcome ~ 1,
          data = test_data,
          family = poisson(),
          prior = "normal")
  )
  
  # Test different priors
  expect_no_error(
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = poisson(),
          prior = "normal",
          seed = 2025)
  )
  expect_no_error(
    fitb <-
      glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
            data = test_data,
            family = poisson(),
            prior = "improper")
  )
  ## Make sure bayes_factors doesn't work for improper prior
  expect_error(bayes_factors(fitb))
  
  
  # Test plot
  expect_s3_class(plot(fita,
                       type = "diagnostics"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_s3_class(plot(fita,
                       type = "pdp"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       variable = "x1",
                       combine_pi_ci = TRUE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       variable = "x1",
                       combine_pi_ci = FALSE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       combine_pi_ci = TRUE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       combine_pi_ci = TRUE),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       combine_pi_ci = FALSE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "pr",
                       variable = "x1"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_s3_class(plot(fita,
                       type = "pr"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "pr",
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "cr",
                       variable = "x1"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "cr"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  
  
  # # Check parallelization
  # plan(multisession,workers = 5)
  # expect_no_error(
  #   glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
  #         data = test_data,
  #         family = poisson(),
  #         seed = 2025)
  # )
  # plan(sequential)
  
  
})


test_that("Test glm_b for poisson data fitting with IS",{
  
  # Generate some data
  set.seed(2025)
  N = 100
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5],
               time = rexp(N))
  test_data$outcome = 
    rpois(N,exp(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e"))) * test_data$time)
  
  
  # Test VB fit
  expect_no_error(
    fita <-
      glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
            data = test_data,
            family = poisson(),
            seed = 2025,
            algorithm = "IS")
  )
  
  # Make sure print works
  expect_no_error(fita)
  
  # Make sure coef works
  expect_type(coef(fita),"double")
  
  # Make sure credint works
  expect_true(is.matrix(credint(fita)))
  
  # Make sure vcov works
  expect_true(is.matrix(vcov(fita)))
  
  # Make sure summary works
  expect_no_error(
    s <- 
      summary(fita)
  )
  ## Check output format
  expect_s3_class(s,c("tbl_df", "tbl", "data.frame"))
  
  expect_identical(colnames(s),
                   c("Variable","Post Mean","Lower","Upper","Prob Dir",
                     "ROPE","ROPE bounds","BF favoring alternative",
                     "Interpretation"))
  expect_type(s$Variable,"character")
  expect_type(s$`Post Mean`,"double")
  expect_type(s$Lower,"double")
  expect_type(s$Upper,"double")
  expect_type(s$`Prob Dir`,"double")
  expect_type(s$ROPE,"double")
  expect_type(s$`ROPE bounds`,"character")
  expect_type(s$`BF favoring alternative`,"double")
  expect_type(s$Interpretation,"character")
  
  # Make sure prediction function works
  expect_no_error(predict(fita))
  expect_no_error(predict(fita,
                          newdata = fita$data[1,]))
  expect_gte(predict(fita,
                     newdata = fita$data[1,],
                     CI_level = 0.8)$CI_lower[1],
             predict(fita,
                     newdata = fita$data[1,],
                     CI_level = 0.9)$CI_lower[1])
  expect_gte(predict(fita,
                     newdata = fita$data[1,],
                     PI_level = 0.8)$PI_lower[1],
             predict(fita,
                     newdata = fita$data[1,],
                     PI_level = 0.9)$PI_lower[1])
  
  # Make sure savage-dickey ratio works
  expect_s3_class(bayes_factors(fita),
                  c("tbl_df", "tbl", "data.frame"))
  expect_s3_class(bayes_factors(fita,by = "v"),
                  c("tbl_df", "tbl", "data.frame"))
  
  # Make sure information criteria work
  expect_type(AIC(fita),"double")
  expect_type(BIC(fita),"double")
  expect_type(DIC(fita),"double")
  expect_type(WAIC(fita),"double")
  
  # Test number of inputs
  expect_no_error(
    glm_b(test_data$outcome ~ test_data$x1,
          family = poisson(),
          prior = "normal",
          algorithm = "IS")
  )
  expect_no_error(
    glm_b(test_data$outcome ~ 1,
          family = poisson(),
          prior = "normal",
          algorithm = "IS")
  )
  expect_no_error(
    glm_b(outcome ~ x1,
          data = test_data,
          family = poisson(),
          prior = "normal",
          algorithm = "IS")
  )
  expect_no_error(
    glm_b(outcome ~ 1,
          data = test_data,
          family = poisson(),
          prior = "normal",
          algorithm = "IS")
  )
  
  # Test different priors
  expect_no_error(
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = poisson(),
          prior = "normal",
          seed = 2025,
          algorithm = "IS")
  )
  expect_no_error(
    fitb <-
      glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
            data = test_data,
            family = poisson(),
            prior = "improper",
            algorithm = "IS")
  )
  ## Make sure bayes_factors doesn't work for improper prior
  expect_error(bayes_factors(fitb))
  
  
  # Test plot
  expect_s3_class(plot(fita,
                       type = "diagnostics"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_s3_class(plot(fita,
                       type = "pdp"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       variable = "x1",
                       combine_pi_ci = TRUE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       variable = "x1",
                       combine_pi_ci = FALSE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       combine_pi_ci = TRUE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       combine_pi_ci = TRUE),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       combine_pi_ci = FALSE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "pr",
                       variable = "x1"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_s3_class(plot(fita,
                       type = "pr"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "pr",
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "cr",
                       variable = "x1"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "cr"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  # # Check parallelization
  # plan(multisession,workers = 5)
  # expect_no_error(
  #   glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
  #         data = test_data,
  #         family = poisson(),
  #         seed = 2025,
  #         algorithm = "IS")
  # )
  # plan(sequential)
  
})


test_that("Test glm_b for poisson data fitting with LSA",{
  
  # Generate some data
  set.seed(2025)
  N = 100
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5],
               time = rexp(N))
  test_data$outcome = 
    rpois(N,exp(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e"))) * test_data$time)
  
  
  # Test VB fit
  expect_no_error(
    fita <-
      glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
            data = test_data,
            family = poisson(),
            seed = 2025)
  )
  
  # Default fit ought to be VB
  expect_identical(fita$algorithm,"VB")
  
  # Make sure print works
  expect_no_error(fita)
  
  # Make sure coef works
  expect_type(coef(fita),"double")
  
  # Make sure credint works
  expect_true(is.matrix(credint(fita)))
  
  # Make sure vcov works
  expect_true(is.matrix(vcov(fita)))
  
  # Make sure summary works
  expect_no_error(
    s <- 
      summary(fita)
  )
  ## Check output format
  expect_s3_class(s,c("tbl_df", "tbl", "data.frame"))
  
  expect_identical(colnames(s),
                   c("Variable","Post Mean","Lower","Upper","Prob Dir",
                     "ROPE","ROPE bounds","BF favoring alternative",
                     "Interpretation"))
  expect_type(s$Variable,"character")
  expect_type(s$`Post Mean`,"double")
  expect_type(s$Lower,"double")
  expect_type(s$Upper,"double")
  expect_type(s$`Prob Dir`,"double")
  expect_type(s$ROPE,"double")
  expect_type(s$`ROPE bounds`,"character")
  expect_type(s$`BF favoring alternative`,"double")
  expect_type(s$Interpretation,"character")
  
  # Make sure prediction function works
  expect_no_error(predict(fita))
  expect_no_error(predict(fita,
                          newdata = fita$data[1,]))
  expect_gte(predict(fita,
                     newdata = fita$data[1,],
                     CI_level = 0.8)$CI_lower[1],
             predict(fita,
                     newdata = fita$data[1,],
                     CI_level = 0.9)$CI_lower[1])
  expect_gte(predict(fita,
                     newdata = fita$data[1,],
                     PI_level = 0.8)$PI_lower[1],
             predict(fita,
                     newdata = fita$data[1,],
                     PI_level = 0.9)$PI_lower[1])
  
  # Make sure savage-dickey ratio works
  expect_s3_class(bayes_factors(fita),
                  c("tbl_df", "tbl", "data.frame"))
  expect_s3_class(bayes_factors(fita,by = "v"),
                  c("tbl_df", "tbl", "data.frame"))
  
  # Make sure information criteria work
  expect_type(AIC(fita),"double")
  expect_type(BIC(fita),"double")
  expect_type(DIC(fita),"double")
  expect_type(WAIC(fita),"double")
  
  # Test number of inputs
  expect_no_error(
    glm_b(test_data$outcome ~ test_data$x1,
          family = poisson(),
          prior = "normal")
  )
  expect_no_error(
    glm_b(test_data$outcome ~ 1,
          family = poisson(),
          prior = "normal")
  )
  expect_no_error(
    glm_b(outcome ~ x1,
          data = test_data,
          family = poisson(),
          prior = "normal")
  )
  expect_no_error(
    glm_b(outcome ~ 1,
          data = test_data,
          family = poisson(),
          prior = "normal")
  )
  
  # Test different priors
  expect_no_error(
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = poisson(),
          prior = "normal",
          seed = 2025)
  )
  expect_no_error(
    fitb <-
      glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
            data = test_data,
            family = poisson(),
            prior = "improper")
  )
  ## Make sure bayes_factors doesn't work for improper prior
  expect_error(bayes_factors(fitb))
  
  
  # Test plot
  expect_s3_class(plot(fita,
                       type = "diagnostics"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_s3_class(plot(fita,
                       type = "pdp"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       variable = "x1",
                       combine_pi_ci = TRUE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       variable = "x1",
                       combine_pi_ci = FALSE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       combine_pi_ci = TRUE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       combine_pi_ci = TRUE),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       combine_pi_ci = FALSE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "pr",
                       variable = "x1"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_s3_class(plot(fita,
                       type = "pr"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "pr",
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "cr",
                       variable = "x1"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "cr"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  
  # # Check parallelization
  # plan(multisession,workers = 5)
  # expect_no_error(
  #   glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
  #         data = test_data,
  #         family = poisson(),
  #         seed = 2025,
  #         algorithm = "LSA")
  # )
  # plan(sequential)
  
})



# Negative Binomial -------------------------------------------------------


test_that("Test glm_b for nbinom data fitting with VB",{
  
  # Generate some data
  set.seed(2025)
  N = 100
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5],
               time = rexp(N))
  test_data$outcome = 
    rnbinom(N,
            mu = exp(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e"))) * test_data$time,
            size = 0.7)
  
  
  # Test VB fit
  expect_no_error(
    fita <-
      glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
            data = test_data,
            family = negbinom(),
            seed = 2025)
  )
  
  # Default fit ought to be VB
  expect_identical(fita$algorithm,"VB")
  
  # Make sure print works
  expect_no_error(fita)
  
  # Make sure coef works
  expect_type(coef(fita),"double")
  
  # Make sure credint works
  expect_true(is.matrix(credint(fita)))
  
  # Make sure vcov works
  expect_true(is.matrix(vcov(fita)))
  
  # Make sure summary works
  expect_no_error(
    s <- 
      summary(fita)
  )
  ## Check output format
  expect_s3_class(s,c("tbl_df", "tbl", "data.frame"))
  
  expect_identical(colnames(s),
                   c("Variable","Post Mean","Lower","Upper","Prob Dir",
                     "ROPE","ROPE bounds","BF favoring alternative",
                     "Interpretation"))
  expect_type(s$Variable,"character")
  expect_type(s$`Post Mean`,"double")
  expect_type(s$Lower,"double")
  expect_type(s$Upper,"double")
  expect_type(s$`Prob Dir`,"double")
  expect_type(s$ROPE,"double")
  expect_type(s$`ROPE bounds`,"character")
  expect_type(s$`BF favoring alternative`,"double")
  expect_type(s$Interpretation,"character")
  
  # Make sure prediction function works
  expect_no_error(predict(fita))
  expect_no_error(predict(fita,
                          newdata = fita$data[1,]))
  expect_gte(predict(fita,
                     newdata = fita$data[1,],
                     CI_level = 0.8)$CI_lower[1],
             predict(fita,
                     newdata = fita$data[1,],
                     CI_level = 0.9)$CI_lower[1])
  expect_gte(predict(fita,
                     newdata = fita$data[1,],
                     PI_level = 0.8)$PI_lower[1],
             predict(fita,
                     newdata = fita$data[1,],
                     PI_level = 0.9)$PI_lower[1])
  
  # Make sure savage-dickey ratio works
  expect_s3_class(bayes_factors(fita),
                  c("tbl_df", "tbl", "data.frame"))
  expect_s3_class(bayes_factors(fita,by = "v"),
                  c("tbl_df", "tbl", "data.frame"))
  
  # Make sure information criteria work
  expect_type(AIC(fita),"double")
  expect_type(BIC(fita),"double")
  expect_type(DIC(fita),"double")
  expect_type(WAIC(fita),"double")
  
  # Test number of inputs
  expect_no_error(
    glm_b(test_data$outcome ~ test_data$x1,
          family = negbinom(),
          prior = "normal")
  )
  expect_no_error(
    glm_b(test_data$outcome ~ 1,
          family = negbinom(),
          prior = "normal")
  )
  expect_no_error(
    glm_b(outcome ~ x1,
          data = test_data,
          family = negbinom(),
          prior = "normal")
  )
  expect_no_error(
    glm_b(outcome ~ 1,
          data = test_data,
          family = negbinom(),
          prior = "normal")
  )
  
  # Test different priors
  expect_no_error(
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = negbinom(),
          prior = "normal",
          seed = 2025)
  )
  expect_no_error(
    fitb <-
      glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
            data = test_data,
            family = negbinom(),
            prior = "improper")
  )
  ## Make sure bayes_factors doesn't work for improper prior
  expect_error(bayes_factors(fitb))
  
  
  # Test plot
  expect_s3_class(plot(fita,
                       type = "diagnostics"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_s3_class(plot(fita,
                       type = "pdp"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       variable = "x1",
                       combine_pi_ci = TRUE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       variable = "x1",
                       combine_pi_ci = FALSE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       combine_pi_ci = TRUE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       combine_pi_ci = TRUE),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       combine_pi_ci = FALSE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "pr",
                       variable = "x1"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_s3_class(plot(fita,
                       type = "pr"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "pr",
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "cr",
                       variable = "x1"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "cr"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  
  
  # # Check parallelization
  # plan(multisession,workers = 5)
  # expect_no_error(
  #   glm_b(outcome ~ x1 + x2 + x3,
  #         data = test_data,
  #         family = negbinom(),
  #         seed = 2025)
  # )
  # plan(sequential)
  
  
})


test_that("Test glm_b for nbinom data fitting with IS",{
  
  # Generate some data
  set.seed(2025)
  N = 100
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5],
               time = rexp(N))
  test_data$outcome = 
    rnbinom(N,
            mu = exp(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e"))) * test_data$time,
            size = 0.7)
  
  
  # Test VB fit
  expect_no_error(
    fita <-
      glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
            data = test_data,
            family = negbinom(),
            seed = 2025,
            algorithm = "IS",
            mc_error = 0.05)
  )
  
  # Make sure print works
  expect_no_error(fita)
  
  # Make sure coef works
  expect_type(coef(fita),"double")
  
  # Make sure credint works
  expect_true(is.matrix(credint(fita)))
  
  # Make sure vcov works
  expect_true(is.matrix(vcov(fita)))
  
  # Make sure summary works
  expect_no_error(
    s <- 
      summary(fita)
  )
  ## Check output format
  expect_s3_class(s,c("tbl_df", "tbl", "data.frame"))
  
  expect_identical(colnames(s),
                   c("Variable","Post Mean","Lower","Upper","Prob Dir",
                     "ROPE","ROPE bounds","BF favoring alternative",
                     "Interpretation"))
  expect_type(s$Variable,"character")
  expect_type(s$`Post Mean`,"double")
  expect_type(s$Lower,"double")
  expect_type(s$Upper,"double")
  expect_type(s$`Prob Dir`,"double")
  expect_type(s$ROPE,"double")
  expect_type(s$`ROPE bounds`,"character")
  expect_type(s$`BF favoring alternative`,"double")
  expect_type(s$Interpretation,"character")
  
  # Make sure prediction function works
  expect_no_error(predict(fita))
  expect_no_error(predict(fita,
                          newdata = fita$data[1,]))
  expect_gte(predict(fita,
                     newdata = fita$data[1,],
                     CI_level = 0.8)$CI_lower[1],
             predict(fita,
                     newdata = fita$data[1,],
                     CI_level = 0.9)$CI_lower[1])
  expect_gte(predict(fita,
                     newdata = fita$data[1,],
                     PI_level = 0.8)$PI_lower[1],
             predict(fita,
                     newdata = fita$data[1,],
                     PI_level = 0.9)$PI_lower[1])
  
  # Make sure savage-dickey ratio works
  expect_s3_class(bayes_factors(fita),
                  c("tbl_df", "tbl", "data.frame"))
  expect_s3_class(bayes_factors(fita,by = "v"),
                  c("tbl_df", "tbl", "data.frame"))
  
  # Make sure information criteria work
  expect_type(AIC(fita),"double")
  expect_type(BIC(fita),"double")
  expect_type(DIC(fita),"double")
  expect_type(WAIC(fita),"double")
  
  # Test number of inputs
  expect_no_error(
    glm_b(test_data$outcome ~ test_data$x1,
          family = negbinom(),
          prior = "normal",
          algorithm = "IS",
          mc_error = 0.05)
  )
  expect_no_error(
    glm_b(test_data$outcome ~ 1,
          family = negbinom(),
          prior = "normal",
          algorithm = "IS",
          mc_error = 0.05)
  )
  expect_no_error(
    glm_b(outcome ~ x1,
          data = test_data,
          family = negbinom(),
          prior = "normal",
          algorithm = "IS",
          mc_error = 0.05)
  )
  expect_no_error(
    glm_b(outcome ~ 1,
          data = test_data,
          family = negbinom(),
          prior = "normal",
          algorithm = "IS",
          mc_error = 0.05)
  )
  
  # Test different priors
  expect_no_error(
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = negbinom(),
          prior = "normal",
          seed = 2025,
          algorithm = "IS",
          mc_error = 0.05)
  )
  expect_no_error(
    fitb <-
      glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
            data = test_data,
            family = negbinom(),
            prior = "improper",
            algorithm = "IS",
            mc_error = 0.05)
  )
  ## Make sure bayes_factors doesn't work for improper prior
  expect_error(bayes_factors(fitb))
  
  
  # Test plot
  expect_s3_class(plot(fita,
                       type = "diagnostics"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_s3_class(plot(fita,
                       type = "pdp"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       variable = "x1",
                       combine_pi_ci = TRUE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       variable = "x1",
                       combine_pi_ci = FALSE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       combine_pi_ci = TRUE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       combine_pi_ci = TRUE),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       combine_pi_ci = FALSE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "pr",
                       variable = "x1"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_s3_class(plot(fita,
                       type = "pr"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "pr",
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "cr",
                       variable = "x1"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "cr"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  # # Check parallelization
  # plan(multisession,workers = 5)
  # expect_no_error(
  #   glm_b(outcome ~ x1 + x2 + x3,
  #         data = test_data,
  #         family = negbinom(),
  #         seed = 2025,
  #         algorithm = "IS",
  #         mc_error = 0.05)
  # )
  # plan(sequential)
  
})


test_that("Test glm_b for nbinom data fitting with LSA",{
  
  # Generate some data
  set.seed(2025)
  N = 100
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5],
               time = rexp(N))
  test_data$outcome = 
    rnbinom(N,
            mu = exp(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e"))) * test_data$time,
            size = 0.7)
  
  
  # Test VB fit
  expect_no_error(
    fita <-
      glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
            data = test_data,
            family = negbinom(),
            seed = 2025,
            algorithm = "LSA")
  )
  
  # Make sure print works
  expect_no_error(fita)
  
  # Make sure coef works
  expect_type(coef(fita),"double")
  
  # Make sure credint works
  expect_true(is.matrix(credint(fita)))
  
  # Make sure vcov works
  expect_true(is.matrix(vcov(fita)))
  
  # Make sure summary works
  expect_no_error(
    s <- 
      summary(fita)
  )
  ## Check output format
  expect_s3_class(s,c("tbl_df", "tbl", "data.frame"))
  
  expect_identical(colnames(s),
                   c("Variable","Post Mean","Lower","Upper","Prob Dir",
                     "ROPE","ROPE bounds","BF favoring alternative",
                     "Interpretation"))
  expect_type(s$Variable,"character")
  expect_type(s$`Post Mean`,"double")
  expect_type(s$Lower,"double")
  expect_type(s$Upper,"double")
  expect_type(s$`Prob Dir`,"double")
  expect_type(s$ROPE,"double")
  expect_type(s$`ROPE bounds`,"character")
  expect_type(s$`BF favoring alternative`,"double")
  expect_type(s$Interpretation,"character")
  
  # Make sure prediction function works
  expect_no_error(predict(fita))
  expect_no_error(predict(fita,
                          newdata = fita$data[1,]))
  expect_gte(predict(fita,
                     newdata = fita$data[1,],
                     CI_level = 0.8)$CI_lower[1],
             predict(fita,
                     newdata = fita$data[1,],
                     CI_level = 0.9)$CI_lower[1])
  expect_gte(predict(fita,
                     newdata = fita$data[1,],
                     PI_level = 0.8)$PI_lower[1],
             predict(fita,
                     newdata = fita$data[1,],
                     PI_level = 0.9)$PI_lower[1])
  
  # Make sure savage-dickey ratio works
  expect_s3_class(bayes_factors(fita),
                  c("tbl_df", "tbl", "data.frame"))
  expect_s3_class(bayes_factors(fita,by = "v"),
                  c("tbl_df", "tbl", "data.frame"))
  
  # Make sure information criteria work
  expect_type(AIC(fita),"double")
  expect_type(BIC(fita),"double")
  expect_type(DIC(fita),"double")
  expect_type(WAIC(fita),"double")
  
  # Test number of inputs
  expect_no_error(
    glm_b(test_data$outcome ~ test_data$x1,
          family = negbinom(),
          prior = "normal",
          algorithm = "LSA")
  )
  expect_no_error(
    glm_b(test_data$outcome ~ 1,
          family = negbinom(),
          prior = "normal",
          algorithm = "LSA")
  )
  expect_no_error(
    glm_b(outcome ~ x1,
          data = test_data,
          family = negbinom(),
          prior = "normal",
          algorithm = "LSA")
  )
  expect_no_error(
    glm_b(outcome ~ 1,
          data = test_data,
          family = negbinom(),
          prior = "normal",
          algorithm = "LSA")
  )
  
  # Test different priors
  expect_no_error(
    glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
          data = test_data,
          family = negbinom(),
          prior = "normal",
          seed = 2025,
          algorithm = "LSA")
  )
  expect_no_error(
    fitb <-
      glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
            data = test_data,
            family = negbinom(),
            prior = "improper",
            algorithm = "LSA")
  )
  ## Make sure bayes_factors doesn't work for improper prior
  expect_error(bayes_factors(fitb))
  
  
  # Test plot
  expect_s3_class(plot(fita,
                       type = "diagnostics"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_s3_class(plot(fita,
                       type = "pdp"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       variable = "x1",
                       combine_pi_ci = TRUE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       variable = "x1",
                       combine_pi_ci = FALSE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       combine_pi_ci = TRUE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       combine_pi_ci = TRUE),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       combine_pi_ci = FALSE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "pr",
                       variable = "x1"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_s3_class(plot(fita,
                       type = "pr"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "pr",
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "cr",
                       variable = "x1"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "cr"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  
  # # Check parallelization
  # plan(multisession,workers = 5)
  # expect_no_error(
  #   glm_b(outcome ~ x1 + x2 + x3,
  #         data = test_data,
  #         family = negbinom(),
  #         seed = 2025,
  #         algorithm = "LSA")
  # )
  # plan(sequential)
  
})



# Gaussian ----------------------------------------------------------------

test_that("Test glm_b for gaussian data. Should pass directly on to lm_b.",{
  
  # Generate some data
  set.seed(2025)
  N = 100
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5])
  test_data$outcome = 
    rnorm(N,-1 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) )
  
  
  expect_no_error(
    fita <- 
      glm_b(outcome ~ x1 + x2 + x3,
            data = test_data,
            prior = "normal",
            family = "gaussian")
  )
  fitb = 
    lm_b(outcome ~ x1 + x2 + x3,
         data = test_data,
         prior = "conjugate")
  expect_equal(fita$summary,
               fitb$summary)
  
  expect_no_error(
    fitc <- 
      glm_b(outcome ~ x1 + x2 + x3,
            data = test_data,
            prior = "zelln",
            family = "gaussian")
  )
  fitd = 
    lm_b(outcome ~ x1 + x2 + x3,
         data = test_data,
         prior = "z")
  expect_equal(fitc$summary,
               fitd$summary)
  
  
})





