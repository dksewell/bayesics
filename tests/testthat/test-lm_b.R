
test_that("Test lm_b with conjugate prior",{
  
  # Generate some data
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5])
  test_data$outcome = 
    rnorm(N,-1 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) )
  
  
  # Test lm_b fit
  expect_no_error(
    fita <-
      lm_b(outcome ~ x1 + x2 + x3,
           data = test_data,
           prior = "conj")
  )
  
  # Test CI_level
  expect_no_error(
    fitb <-
      lm_b(outcome ~ x1 + x2 + x3,
           data = test_data,
           prior = "conj",
           CI_level = 0.9)
  )
  expect_lt(fita$summary$Lower[1],
            fitb$summary$Lower[1])
  
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
  
  # Test number and type of inputs
  expect_no_error(
    lm_b(test_data$outcome ~ test_data$x1,
         prior = "conjugate")
  )
  expect_no_error(
    lm_b(test_data$outcome ~ 1,
         prior = "conju")
  )
  expect_no_error(
    lm_b(outcome ~ x1,
         data = test_data,
         prior = "con")
  )
  expect_no_error(
    lm_b(outcome ~ 1,
         data = test_data,
         prior = "co")
  )
  
  # Test hyperparameters
  expect_no_error(
    fitc <-
      lm_b(outcome ~ x1 + x2 + x3,
           data = test_data,
           prior = "conj",
           prior_beta_mean = rep(2,7),
           prior_beta_precision = diag(7) * 0.1,
           prior_var_shape = 2,
           prior_var_rate = 2)
  )
  expect_equal(fitc$hyperparameters$mu, rep(2,7))
  expect_equal(fitc$hyperparameters$a, 2)
  expect_equal(fitc$hyperparameters$b, 2)
  expect_equal(fitc$hyperparameters$V, diag(7) * 0.1)
  expect_true(!isTRUE(all.equal(fita$summary,
                                fitc$summary)))
  
  
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
                       type = c("ci","pi"),
                       variable = "x1",
                       combine_pi_ci = TRUE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("ci","pi"),
                       variable = "x1",
                       combine_pi_ci = FALSE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("ci","pi"),
                       combine_pi_ci = TRUE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("ci","pi"),
                       combine_pi_ci = TRUE),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("ci","pi"),
                       combine_pi_ci = FALSE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "pi",
                       variable = "x1"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "pi"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "pi",
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "ci",
                       variable = "x1"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "ci"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  
  
  # Test if response transformation works
  test_data$e_outcome = exp(test_data$outcome)
  
  ## Test lm_b fit
  expect_no_error(
    fitc <-
      lm_b(log(e_outcome) ~ x1 + x2 + x3,
           data = test_data,
           prior = "conj")
  )
  expect_equal(fita$summary,
               fitc$summary)
  
  ## Make sure print works
  expect_no_error(fitc)
  
  ## Make sure coef works
  expect_type(coef(fitc),"double")
  
  ## Make sure credint works
  expect_true(is.matrix(credint(fitc)))
  
  ## Make sure vcov works
  expect_true(is.matrix(vcov(fitc)))
  
  ## Make sure summary works
  expect_no_error(
    s <- 
      summary(fitc)
  )
  ### Check output format
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
  
  ## Make sure prediction function works
  expect_no_error(
    fitc_preds <- predict(fitc)
  )
  expect_no_error(predict(fitc,
                          newdata = fitc$data[1,]))
  expect_gt(predict(fitc,
                     newdata = fitc$data[1,],
                     CI_level = 0.8)$CI_lower[1],
             predict(fitc,
                     newdata = fitc$data[1,],
                     CI_level = 0.9)$CI_lower[1])
  expect_gt(predict(fitc,
                     newdata = fitc$data[1,],
                     PI_level = 0.8)$PI_lower[1],
             predict(fitc,
                     newdata = fitc$data[1,],
                     PI_level = 0.9)$PI_lower[1])
  
  ## Check plot
  expect_s3_class(plot(fitc,
                       type = "pdp"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fitc),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
})


test_that("Test lm_b with zellner's g prior",{
  
  
  # Generate some data
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5])
  test_data$outcome = 
    rnorm(N,-1 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) )
  
  
  # Test lm_b fit
  expect_no_error(
    fita <-
      lm_b(outcome ~ x1 + x2 + x3,
           data = test_data,
           prior = "ze")
  )
  
  # Test CI_level
  expect_no_error(
    fitb <-
      lm_b(outcome ~ x1 + x2 + x3,
           data = test_data,
           prior = "ze",
           CI_level = 0.9)
  )
  expect_lt(fita$summary$Lower[1],
            fitb$summary$Lower[1])
  
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
  
  # Test number and type of inputs
  expect_no_error(
    lm_b(test_data$outcome ~ test_data$x1,
         prior = "zellner")
  )
  expect_no_error(
    lm_b(test_data$outcome ~ 1,
         prior = "zellne")
  )
  expect_no_error(
    lm_b(outcome ~ x1,
         data = test_data,
         prior = "zell")
  )
  expect_no_error(
    lm_b(outcome ~ 1,
         data = test_data,
         prior = "z")
  )
  
  # Test hyperparameters
  ## Test beta hyperparameters do NOT affect anything for zellner's g
  expect_no_error(
    fitc <-
      lm_b(outcome ~ x1 + x2 + x3,
           data = test_data,
           prior = "z",
           prior_beta_mean = rep(2,7),
           prior_beta_precision = diag(7) * 0.1)
  )
  expect_true(isTRUE(all.equal(fita$summary,
                               fitc$summary)))
  ## Test g 
  expect_no_error(
    fitd <-
      lm_b(outcome ~ x1 + x2 + x3,
           data = test_data,
           prior = "z",
           zellner_g = 10)
  )
  expect_true(!isTRUE(all.equal(fita$summary,
                                fitd$summary)))
  
  
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
                       type = c("ci","pi"),
                       variable = "x1",
                       combine_pi_ci = TRUE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("ci","pi"),
                       variable = "x1",
                       combine_pi_ci = FALSE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("ci","pi"),
                       combine_pi_ci = TRUE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("ci","pi"),
                       combine_pi_ci = TRUE),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("ci","pi"),
                       combine_pi_ci = FALSE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "pi",
                       variable = "x1"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "pi"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "pi",
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "ci",
                       variable = "x1"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "ci"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  
})


test_that("Test lm_b with improper prior",{
  
  
  # Generate some data
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5])
  test_data$outcome = 
    rnorm(N,-1 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) )
  
  
  # Test lm_b fit
  expect_no_error(
    fita <-
      lm_b(outcome ~ x1 + x2 + x3,
           data = test_data,
           prior = "improper")
  )
  
  # Test CI_level
  expect_no_error(
    fitb <-
      lm_b(outcome ~ x1 + x2 + x3,
           data = test_data,
           prior = "imp",
           CI_level = 0.9)
  )
  expect_lt(fita$summary$Lower[1],
            fitb$summary$Lower[1])
  
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
                     "ROPE","ROPE bounds"))
  expect_type(s$Variable,"character")
  expect_type(s$`Post Mean`,"double")
  expect_type(s$Lower,"double")
  expect_type(s$Upper,"double")
  expect_type(s$`Prob Dir`,"double")
  expect_type(s$ROPE,"double")
  expect_type(s$`ROPE bounds`,"character")
  
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
  
  # Make sure savage-dickey ratio does NOT work
  expect_error(bayes_factors(fita))
  
  # Make sure information criteria work
  expect_type(AIC(fita),"double")
  expect_type(BIC(fita),"double")
  expect_type(DIC(fita),"double")
  expect_type(WAIC(fita),"double")
  
  # Test number and type of inputs
  expect_no_error(
    lm_b(test_data$outcome ~ test_data$x1,
         prior = "improper")
  )
  expect_no_error(
    lm_b(test_data$outcome ~ 1,
         prior = "improp")
  )
  expect_no_error(
    lm_b(outcome ~ x1,
         data = test_data,
         prior = "imp")
  )
  expect_no_error(
    lm_b(outcome ~ 1,
         data = test_data,
         prior = "i")
  )
  
  # Test hyperparameters
  expect_no_error(
    fitc <-
      lm_b(outcome ~ x1 + x2 + x3,
           data = test_data,
           prior = "imp",
           prior_beta_mean = rep(2,7),
           prior_beta_precision = diag(7) * 0.1,
           prior_var_shape = 2,
           prior_var_rate = 2)
  )
  expect_true(is.na(fitc$hyperparameters))
  expect_true(isTRUE(all.equal(fita$summary,
                               fitc$summary)))
  
  
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
                       type = c("ci","pi"),
                       variable = "x1",
                       combine_pi_ci = TRUE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("ci","pi"),
                       variable = "x1",
                       combine_pi_ci = FALSE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("ci","pi"),
                       combine_pi_ci = TRUE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("ci","pi"),
                       combine_pi_ci = TRUE),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("ci","pi"),
                       combine_pi_ci = FALSE,
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "pi",
                       variable = "x1"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "pi"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "pi",
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "ci",
                       variable = "x1"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "ci"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  
})



test_that("Test lm_b with weights",{
  
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
  expect_no_error(
    fita <- 
      lm_b(outcome ~ x1 + x2 + x3,
           data = test_data,
           weights = test_data$w,
           prior = "conj")
  )
  expect_no_error(
    fitb <-
      lm_b(outcome ~ x1 + x2 + x3,
           data = test_data,
           weights = "w",
           prior = "conj")
  )
  expect_no_error(
    fitc <-
      lm_b(outcome ~ x1 + x2 + x3,
           data = test_data,
           prior = "conj")
  )
  expect_equal(fita$summary,
               fitb$summary)
  expect_true(!isTRUE(all.equal(fita$summary,
                                fitc$summary)))
  
  # Zellner's g prior
  expect_no_error(
    fitd <- 
      lm_b(outcome ~ x1 + x2 + x3,
           data = test_data,
           weights = test_data$w,
           prior = "z")
  )
  expect_no_error(
    fite <-
      lm_b(outcome ~ x1 + x2 + x3,
           data = test_data,
           weights = "w",
           prior = "z")
  )
  expect_no_error(
    fitf <-
      lm_b(outcome ~ x1 + x2 + x3,
           data = test_data,
           prior = "z")
  )
  expect_equal(fitd$summary,
               fite$summary)
  expect_true(!isTRUE(all.equal(fitd$summary,
                                fitf$summary)))
  
  # Improper prior
  expect_no_error(
    fitg <- 
      lm_b(outcome ~ x1 + x2 + x3,
           data = test_data,
           weights = test_data$w,
           prior = "i")
  )
  expect_no_error(
    fith <-
      lm_b(outcome ~ x1 + x2 + x3,
           data = test_data,
           weights = "w",
           prior = "i")
  )
  expect_no_error(
    fiti <-
      lm_b(outcome ~ x1 + x2 + x3,
           data = test_data,
           prior = "i")
  )
  expect_equal(fitg$summary,
               fith$summary)
  expect_true(!isTRUE(all.equal(fitg$summary,
                                fiti$summary)))
  
})


test_that("Test complicated terms in lm_b formula",{
  
  
  # I(...) in covariate and interactions
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
            2 * (test_data$x3 %in% c("d","e")) )
  
  ## Test lm_b fit
  expect_no_error(
    fita <-
      lm_b(outcome ~ x1 * x2 + x3 + I(x1^2),
           data = test_data,
           prior = "conj")
  )
  
  ## Make sure print works
  expect_no_error(fita)
  
  ## Make sure coef works
  expect_type(coef(fita),"double")
  
  ## Make sure credint works
  expect_true(is.matrix(credint(fita)))
  
  ## Make sure vcov works
  expect_true(is.matrix(vcov(fita)))
  
  ## Make sure summary works
  expect_no_error(
    s <- 
      summary(fita)
  )
  ### Check output format
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
  
  ## Make sure prediction function works
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
  
  ## Make sure savage-dickey ratio works
  expect_no_error(
    sd1 <-
      bayes_factors(fita)
  )
  expect_s3_class(sd1,
    c("tbl_df", "tbl", "data.frame"))
  expect_identical(sd1$Variable,
                   c("(Intercept)","x1","x2","x3b","x3c","x3d","x3e","I(x1^2)","x1:x2"))
  
  expect_no_error(
    sd2 <-
      bayes_factors(fita, by = "variabl")
  )
  expect_s3_class(sd2,
                  c("tbl_df", "tbl", "data.frame"))
  expect_identical(sd2$Variable,
                   c("x1","x2","x3"))
  
  
  ## Make sure information criteria work
  expect_type(AIC(fita),"double")
  expect_type(BIC(fita),"double")
  expect_type(DIC(fita),"double")
  expect_type(WAIC(fita),"double")
  
  
  ## Test plot
  expect_s3_class(plot(fita,
                       type = "pdp"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  
  
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
  
  ## Test lm_b fit
  expect_no_error(
    fitb <-
      lm_b(outcome ~ ns(x1,df = 5) + x2,
           data = test_data,
           prior = "conj")
  )
  
  ## Make sure print works
  expect_no_error(fitb)
  
  ## Make sure coef works
  expect_type(coef(fitb),"double")
  
  ## Make sure credint works
  expect_true(is.matrix(credint(fitb)))
  
  ## Make sure vcov works
  expect_true(is.matrix(vcov(fitb)))
  
  ## Make sure summary works
  expect_no_error(
    s <- 
      summary(fitb)
  )
  ### Check output format
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
  
  ## Make sure prediction function works
  expect_no_error(predict(fitb))
  expect_no_error(predict(fitb,
                          newdata = fitb$data[1,]))
  expect_gte(predict(fitb,
                     newdata = fitb$data[1,],
                     CI_level = 0.8)$CI_lower[1],
             predict(fitb,
                     newdata = fitb$data[1,],
                     CI_level = 0.9)$CI_lower[1])
  expect_gte(predict(fitb,
                     newdata = fitb$data[1,],
                     PI_level = 0.8)$PI_lower[1],
             predict(fitb,
                     newdata = fitb$data[1,],
                     PI_level = 0.9)$PI_lower[1])
  
  ## Make sure savage-dickey ratio works
  expect_no_error(
    sd1 <-
      bayes_factors(fitb)
  )
  expect_s3_class(sd1,
                  c("tbl_df", "tbl", "data.frame"))
  expect_identical(sd1$Variable,
                   c("(Intercept)",
                     paste("ns(x1, df = 5)",1:5,sep=""),
                     "x2"))
  expect_no_error(
    sd2 <-
      bayes_factors(fitb, by = "variabl")
  )
  expect_s3_class(sd2,
                  c("tbl_df", "tbl", "data.frame"))
  expect_identical(sd2$Variable,
                   c("x1","x2"))
  
  ## Make sure information criteria work
  expect_type(AIC(fitb),"double")
  expect_type(BIC(fitb),"double")
  expect_type(DIC(fitb),"double")
  expect_type(WAIC(fitb),"double")
  
  
  ## Test plot
  expect_s3_class(plot(fitb,
                       type = "pdp"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fitb),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  
  
  # Test if response transformation works
  test_data$e_outcome = exp(test_data$outcome)
  
  ## Test lm_b fit
  expect_no_error(
    fitc <-
      lm_b(log(e_outcome) ~ x1 * x2 + x3 + I(x1^2),
           data = test_data,
           prior = "conj")
  )
  expect_no_error(
    fitc2 <-
      lm_b(outcome ~ x1 * x2 + x3 + I(x1^2),
           data = test_data,
           prior = "conj")
  )
  expect_equal(fitc$summary,
               fitc2$summary)
  
  ## Make sure print works
  expect_no_error(fitc)
  
  ## Make sure coef works
  expect_type(coef(fitc),"double")
  
  ## Make sure credint works
  expect_true(is.matrix(credint(fitc)))
  
  ## Make sure vcov works
  expect_true(is.matrix(vcov(fitc)))
  
  ## Make sure summary works
  expect_no_error(
    s <- 
      summary(fitc)
  )
  ### Check output format
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
  
  ## Make sure prediction function works
  expect_no_error(
    fitc_preds <- predict(fitc)
  )
  expect_no_error(predict(fitc,
                          newdata = fitc$data[1,]))
  expect_gte(predict(fitc,
                     newdata = fitc$data[1,],
                     CI_level = 0.8)$CI_lower[1],
             predict(fitc,
                     newdata = fitc$data[1,],
                     CI_level = 0.9)$CI_lower[1])
  expect_gte(predict(fitc,
                     newdata = fitc$data[1,],
                     PI_level = 0.8)$PI_lower[1],
             predict(fitc,
                     newdata = fitc$data[1,],
                     PI_level = 0.9)$PI_lower[1])
  
  ## Make sure savage-dickey ratio works
  expect_no_error(
    sd1 <-
      bayes_factors(fitc)
  )
  expect_s3_class(sd1,
                  c("tbl_df", "tbl", "data.frame"))
  expect_identical(sd1$Variable,
                   c("(Intercept)","x1","x2","x3b","x3c","x3d","x3e","I(x1^2)","x1:x2"))
  
  expect_no_error(
    sd2 <-
      bayes_factors(fitc, by = "variabl")
  )
  expect_s3_class(sd2,
                  c("tbl_df", "tbl", "data.frame"))
  expect_identical(sd2$Variable,
                   c("x1","x2","x3"))
  
  
  ## Make sure information criteria work
  expect_type(AIC(fitc),"double")
  expect_type(BIC(fitc),"double")
  expect_type(DIC(fitc),"double")
  expect_type(WAIC(fitc),"double")
  
  
  ## Test plot
  expect_s3_class(plot(fitc,
                       type = "pdp"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fitc),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  
  
  # ~ .
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5])
  test_data$outcome = 
    rnorm(N,-1 + 
            test_data$x1 + 
            2 * (test_data$x3 %in% c("d","e")) )
  
  ## Test lm_b fit
  expect_no_error(
    fitd <-
      lm_b(outcome ~ .,
           data = test_data,
           prior = "conj")
  )
  
  ## Make sure print works
  expect_no_error(fitd)
  
  ## Make sure coef works
  expect_type(coef(fitd),"double")
  
  ## Make sure credint works
  expect_true(is.matrix(credint(fitd)))
  
  ## Make sure vcov works
  expect_true(is.matrix(vcov(fitd)))
  
  ## Make sure summary works
  expect_no_error(
    s <- 
      summary(fitd)
  )
  ### Check output format
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
  
  ## Make sure prediction function works
  expect_no_error(predict(fitd))
  expect_no_error(predict(fitd,
                          newdata = fitd$data[1,]))
  expect_gte(predict(fitd,
                     newdata = fitd$data[1,],
                     CI_level = 0.8)$CI_lower[1],
             predict(fitd,
                     newdata = fitd$data[1,],
                     CI_level = 0.9)$CI_lower[1])
  expect_gte(predict(fitd,
                     newdata = fitd$data[1,],
                     PI_level = 0.8)$PI_lower[1],
             predict(fitd,
                     newdata = fitd$data[1,],
                     PI_level = 0.9)$PI_lower[1])
  
  ## Make sure savage-dickey ratio works
  expect_no_error(
    sd1 <-
      bayes_factors(fitd)
  )
  expect_s3_class(sd1,
                  c("tbl_df", "tbl", "data.frame"))
  expect_identical(sd1$Variable,
                   c("(Intercept)","x1","x2","x3b","x3c","x3d","x3e"))
  
  expect_no_error(
    sd2 <-
      bayes_factors(fitd, by = "variabl")
  )
  expect_s3_class(sd2,
                  c("tbl_df", "tbl", "data.frame"))
  expect_identical(sd2$Variable,
                   c("x1","x2","x3"))
  
  
  ## Make sure information criteria work
  expect_type(AIC(fitd),"double")
  expect_type(BIC(fitd),"double")
  expect_type(DIC(fitd),"double")
  expect_type(WAIC(fitd),"double")
  
  
  ## Test plot
  expect_s3_class(plot(fitd,
                       type = "pdp"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fitd),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  
  
  
  
  
})
