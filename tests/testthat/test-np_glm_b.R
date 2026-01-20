
go_fast_for_cran_checks = TRUE

# Binomial ----------------------------------------------------------------

if(!go_fast_for_cran_checks){
  test_that("Test np_glm_b for binomial data fitting with bootstrapping",{
  
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
        np_glm_b(outcome ~ x1 + x2 + x3,
                 data = test_data,
                 family = binomial(),
                 seed = 2025,
                 n_draws = 50,
                 mc_error = 0.1,
                 ask_before_full_sampling = FALSE)
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
    expect_no_error(
      s2 <-
        summary(fita,
                interpretable = FALSE)
    )
    expect_equal(s$`Post Mean`,
                 exp(s2$`Post Mean`[-1]))
  
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
  
  
    # Test number of inputs
    expect_no_error(
      np_glm_b(test_data$outcome ~ test_data$x1,
            family = binomial(),
            n_draws = 50,
            mc_error = 0.2,
            ask_before_full_sampling = FALSE)
    )
    expect_no_error(
      np_glm_b(test_data$outcome ~ 1,
            family = binomial(),
            n_draws = 50,
            mc_error = 0.2,
            ask_before_full_sampling = FALSE)
    )
    expect_no_error(
      np_glm_b(outcome ~ x1,
            data = test_data,
            family = binomial(),
            n_draws = 50,
            mc_error = 0.2,
            ask_before_full_sampling = FALSE)
    )
    expect_no_error(
      np_glm_b(outcome ~ 1,
            data = test_data,
            family = binomial(),
            n_draws = 50,
            mc_error = 0.2,
            ask_before_full_sampling = FALSE)
    )
  
  
    # Test plot
    expect_s3_class(plot(fita,
                         type = "pdp"),
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
    expect_s3_class(plot(fita,
                         type = "cr",
                         exemplar_covariates = fita$data[1,]),
                    c("patchwork","ggplot2::ggplot","ggplot",
                      "ggplot2::gg","S7_object","gg"))
    expect_s3_class(plot(fita),
                    c("patchwork","ggplot2::ggplot","ggplot",
                      "ggplot2::gg","S7_object","gg"))
  
  
    # Check parallelization
    plan(multisession,workers = 5)
    expect_no_error(
      fita <-
        np_glm_b(outcome ~ x1 + x2 + x3,
                 data = test_data,
                 family = binomial(),
                 seed = 2025,
                 n_draws = 100,
                 mc_error = 0.1,
                 ask_before_full_sampling = FALSE)
    )
    plan(sequential)
  
  
  })
}


test_that("Test np_glm_b for binomial data fitting with LSA",{
  
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
      np_glm_b(outcome ~ x1 + x2 + x3,
               data = test_data,
               family = binomial())
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
  expect_no_error(
    s2 <- 
      summary(fita,
              interpretable = FALSE)
  )
  expect_equal(s$`Post Mean`,
               exp(s2$`Post Mean`[-1]))
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
  
  
  # Test number of inputs
  expect_no_error(
    np_glm_b(test_data$outcome ~ test_data$x1,
             family = binomial())
  )
  expect_no_error(
    np_glm_b(test_data$outcome ~ 1,
             family = binomial())
  )
  expect_no_error(
    np_glm_b(outcome ~ x1,
             data = test_data,
             family = binomial())
  )
  expect_no_error(
    np_glm_b(outcome ~ 1,
             data = test_data,
             family = binomial())
  )
  
  
  # Test plot
  expect_s3_class(plot(fita,
                       type = "pdp"),
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
  expect_s3_class(plot(fita,
                       type = "cr",
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  
  
})


test_that("Test np_glm_b for binomial data with >1 trials",{
  
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
      np_glm_b(outcome ~ x1 + x2 + x3,
               data = test_data,
               trials = "n_trials",
               family = binomial())
  )
  expect_no_error(
    fitb <-
      np_glm_b(outcome ~ x1 + x2 + x3,
               data = test_data,
               trials = test_data$n_trials,
               family = binomial())
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

if(!go_fast_for_cran_checks){
  test_that("Test np_glm_b for poisson data fitting with bootstrapping",{
  
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
        np_glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
                 data = test_data,
                 family = poisson(),
                 seed = 2025,
                 n_draws = 50,
                 mc_error = 0.1,
                 ask_before_full_sampling = FALSE)
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
    expect_no_error(
      s2 <-
        summary(fita,
                interpretable = FALSE)
    )
    expect_equal(s$`Post Mean`,
                 exp(s2$`Post Mean`[-1]))
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
  
  
    # Test number of inputs
    expect_no_error(
      np_glm_b(test_data$outcome ~ test_data$x1,
               family = poisson(),
               n_draws = 50,
               mc_error = 0.2,
               ask_before_full_sampling = FALSE)
    )
    expect_no_error(
      np_glm_b(test_data$outcome ~ 1,
               family = poisson(),
               n_draws = 50,
               mc_error = 0.2,
               ask_before_full_sampling = FALSE)
    )
    expect_no_error(
      np_glm_b(outcome ~ x1,
               data = test_data,
               family = poisson(),
               n_draws = 50,
               mc_error = 0.2,
               ask_before_full_sampling = FALSE)
    )
    expect_no_error(
      np_glm_b(outcome ~ 1,
               data = test_data,
               family = poisson(),
               n_draws = 50,
               mc_error = 0.2,
               ask_before_full_sampling = FALSE)
    )
  
  
    # Test plot
    expect_s3_class(plot(fita,
                         type = "pdp"),
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
    expect_s3_class(plot(fita,
                         type = "cr",
                         exemplar_covariates = fita$data[1,]),
                    c("patchwork","ggplot2::ggplot","ggplot",
                      "ggplot2::gg","S7_object","gg"))
    expect_s3_class(plot(fita),
                    c("patchwork","ggplot2::ggplot","ggplot",
                      "ggplot2::gg","S7_object","gg"))
  
  
    # # Check parallelization
    # plan(multisession,workers = 5)
    # expect_no_error(
    #   fita <-
    #     np_glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
    #              data = test_data,
    #              family = poisson(),
    #              seed = 2025,
    #              n_draws = 100,
    #              mc_error = 0.1,
    #              ask_before_full_sampling = FALSE)
    # )
    # plan(sequential)
  
  
  })
}


test_that("Test np_glm_b for poisson data fitting with LSA",{
  
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
      np_glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
               data = test_data,
               family = poisson())
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
  expect_no_error(
    s2 <- 
      summary(fita,
              interpretable = FALSE)
  )
  expect_equal(s$`Post Mean`,
               exp(s2$`Post Mean`[-1]))
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
  
  
  # Test number of inputs
  expect_no_error(
    np_glm_b(test_data$outcome ~ test_data$x1,
             family = poisson())
  )
  expect_no_error(
    np_glm_b(test_data$outcome ~ 1,
             family = poisson())
  )
  expect_no_error(
    np_glm_b(outcome ~ x1,
             data = test_data,
             family = poisson())
  )
  expect_no_error(
    np_glm_b(outcome ~ 1,
             data = test_data,
             family = poisson())
  )
  
  
  # Test plot
  expect_s3_class(plot(fita,
                       type = "pdp"),
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
  expect_s3_class(plot(fita,
                       type = "cr",
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  
  
})



# Negative Binomial -------------------------------------------------------

if(!go_fast_for_cran_checks){
  test_that("Test np_glm_b for negative binomial data fitting with bootstrapping",{
  
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
        np_glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
                 data = test_data,
                 family = negbinom(),
                 seed = 2025,
                 n_draws = 50,
                 mc_error = 0.1,
                 ask_before_full_sampling = FALSE)
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
    expect_no_error(
      s2 <-
        summary(fita,
                interpretable = FALSE)
    )
    expect_equal(s$`Post Mean`,
                 exp(s2$`Post Mean`[-1]))
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
  
  
    # Test number of inputs
    expect_no_error(
      np_glm_b(test_data$outcome ~ test_data$x1,
               family = negbinom(),
               n_draws = 50,
               mc_error = 0.2,
               ask_before_full_sampling = FALSE)
    )
    expect_no_error(
      np_glm_b(test_data$outcome ~ 1,
               family = negbinom(),
               n_draws = 50,
               mc_error = 0.5,
               ask_before_full_sampling = FALSE)
    )
    expect_no_error(
      np_glm_b(outcome ~ x1,
               data = test_data,
               family = negbinom(),
               n_draws = 50,
               mc_error = 0.2,
               ask_before_full_sampling = FALSE)
    )
    expect_no_error(
      np_glm_b(outcome ~ 1,
               data = test_data,
               family = negbinom(),
               n_draws = 50,
               mc_error = 0.5,
               ask_before_full_sampling = FALSE)
    )
  
  
    # Test plot
    expect_s3_class(plot(fita,
                         type = "pdp"),
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
    expect_s3_class(plot(fita,
                         type = "cr",
                         exemplar_covariates = fita$data[1,]),
                    c("patchwork","ggplot2::ggplot","ggplot",
                      "ggplot2::gg","S7_object","gg"))
    expect_s3_class(plot(fita),
                    c("patchwork","ggplot2::ggplot","ggplot",
                      "ggplot2::gg","S7_object","gg"))
  
  
    # # Check parallelization
    # plan(multisession,workers = 5)
    # expect_no_error(
    #   fita <-
    #     np_glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
    #              data = test_data,
    #              family = negbinom(),
    #              seed = 2025,
    #              n_draws = 100,
    #              mc_error = 0.2,
    #              ask_before_full_sampling = FALSE)
    # )
    # plan(sequential)
  
  
  })
}


test_that("Test np_glm_b for negative binomial data fitting with LSA",{
  
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
      np_glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
               data = test_data,
               family = negbinom())
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
  expect_no_error(
    s2 <- 
      summary(fita,
              interpretable = FALSE)
  )
  expect_equal(s$`Post Mean`,
               exp(s2$`Post Mean`[-1]))
  
  
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
  
  
  # Test number of inputs
  expect_no_error(
    np_glm_b(test_data$outcome ~ test_data$x1,
             family = negbinom())
  )
  expect_no_error(
    np_glm_b(test_data$outcome ~ 1,
             family = negbinom())
  )
  expect_no_error(
    np_glm_b(outcome ~ x1,
             data = test_data,
             family = negbinom())
  )
  expect_no_error(
    np_glm_b(outcome ~ 1,
             data = test_data,
             family = negbinom())
  )
  
  
  # Test plot
  expect_s3_class(plot(fita,
                       type = "pdp"),
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
  expect_s3_class(plot(fita,
                       type = "cr",
                       exemplar_covariates = fita$data[1,]),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  
  
})



# Gaussian ----------------------------------------------------------------

if(!go_fast_for_cran_checks){
  test_that("Test np_glm_b for gaussian data fitting with bootstrapping",{
  
    # Generate some data
    set.seed(2025)
    N = 100
    test_data =
      data.frame(x1 = rnorm(N),
                 x2 = rnorm(N),
                 x3 = letters[1:5])
    test_data$outcome =
      rnorm(N,-1 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) )
  
  
    # Test fit
    expect_no_error(
      fita <-
        np_glm_b(outcome ~ x1 + x2 + x3,
                 data = test_data,
                 family = gaussian(),
                 seed = 2025,
                 n_draws = 50,
                 mc_error = 0.1,
                 ask_before_full_sampling = FALSE)
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
    expect_no_error(
      s2 <-
        summary(fita,
                interpretable = FALSE)
    )
    expect_equal(s,
                 s2)
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
  
  
    # Test number of inputs
    expect_no_error(
      np_glm_b(test_data$outcome ~ test_data$x1,
               family = gaussian(),
               n_draws = 50,
               mc_error = 0.2,
               ask_before_full_sampling = FALSE)
    )
    expect_no_error(
      np_glm_b(test_data$outcome ~ 1,
               family = gaussian(),
               n_draws = 50,
               mc_error = 0.2,
               ask_before_full_sampling = FALSE)
    )
    expect_no_error(
      np_glm_b(outcome ~ x1,
               data = test_data,
               family = gaussian(),
               n_draws = 50,
               mc_error = 0.2,
               ask_before_full_sampling = FALSE)
    )
    expect_no_error(
      np_glm_b(outcome ~ 1,
               data = test_data,
               family = gaussian(),
               n_draws = 50,
               mc_error = 0.2,
               ask_before_full_sampling = FALSE)
    )
  
  
    # Test plot
    expect_s3_class(plot(fita,
                         type = "pdp"),
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
    expect_s3_class(plot(fita,
                         type = "cr",
                         exemplar_covariates = fita$data[1,]),
                    c("patchwork","ggplot2::ggplot","ggplot",
                      "ggplot2::gg","S7_object","gg"))
    expect_s3_class(plot(fita),
                    c("patchwork","ggplot2::ggplot","ggplot",
                      "ggplot2::gg","S7_object","gg"))
  
    # Test if response transformation works
    test_data$e_outcome = exp(test_data$outcome)
  
    ## Test np_glm_b with transformed response
    expect_no_error(
      fitb <-
        np_glm_b(log(e_outcome) ~ x1 + x2 + x3,
                 data = test_data,
                 family = gaussian(),
                 seed = 2025,
                 n_draws = 50,
                 mc_error = 0.1,
                 ask_before_full_sampling = FALSE)
    )
    expect_equal(fita$summary,
                 fitb$summary)
  
    ## Make sure prediction function works
    expect_no_error(
      predict(fitb)
    )
    expect_no_error(
      predict(fitb, newdata = fitb$data[1,])
    )
    expect_s3_class(plot(fitb),
                    c("patchwork","ggplot2::ggplot","ggplot",
                      "ggplot2::gg","S7_object","gg"))
  
  
  
    # Check parallelization
    plan(multisession,workers = 5)
    expect_no_error(
      fita <-
        np_glm_b(outcome ~ x1 + x2 + x3,
                 data = test_data,
                 family = gaussian(),
                 seed = 2025,
                 n_draws = 100,
                 mc_error = 0.1,
                 ask_before_full_sampling = FALSE)
    )
    plan(sequential)
  
  
  })
}


test_that("Test np_glm_b for gaussian data fitting with LSA",{
  
  # Generate some data
  set.seed(2025)
  N = 100
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5])
  test_data$outcome = 
    rnorm(N,-1 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) )
  
  
  # Test VB fit
  expect_no_error(
    fita <-
      np_glm_b(outcome ~ x1 + x2 + x3,
               data = test_data,
               family = gaussian())
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
  expect_no_error(
    s2 <- 
      summary(fita,
              interpretable = FALSE)
  )
  expect_equal(s,
               s2)
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
  
  
  # Test number of inputs
  expect_no_error(
    np_glm_b(test_data$outcome ~ test_data$x1,
             family = gaussian())
  )
  expect_no_error(
    np_glm_b(test_data$outcome ~ 1,
             family = gaussian())
  )
  expect_no_error(
    np_glm_b(outcome ~ x1,
             data = test_data,
             family = gaussian())
  )
  expect_no_error(
    np_glm_b(outcome ~ 1,
             data = test_data,
             family = gaussian())
  )
  
  
  # Test plot
  expect_s3_class(plot(fita,
                       type = "pdp"),
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
  
  
  
})


