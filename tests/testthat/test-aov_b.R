

# Proper, heteroscedastic -------------------------------------------------

test_that("Proper prior and heteroscedastic model works", {
  
  # Create data
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rep(letters[1:5],N/5))
  test_data$outcome = 
    rnorm(N,-1 + 2 * (test_data$x1 %in% c("d","e")) )
  
  # No errors upon fitting
  expect_no_error(
    fita <-
      aov_b(outcome ~ x1,
            test_data,
            prior_mean_mu = 2,
            prior_mean_nu = 0.5,
            prior_var_shape = 0.01,
            prior_var_rate = 0.01)
    )
  
  # Make sure print works
  expect_no_error(fita)
  
  # Make sure summary.aov_b works
  expect_no_error(
    s <- 
      summary(fita)
  )
  ## Check output format
  expect_type(s,"list")
  expect_s3_class(s$summary,c("tbl_df", "tbl", "data.frame"))
  
  expect_named(s$summary,c("Variable","Post Mean","Lower","Upper","Prob Dir"))
  expect_type(s$summary$Variable,"character")
  expect_type(s$summary$`Post Mean`,"double")
  expect_type(s$summary$Lower,"double")
  expect_type(s$summary$Upper,"double")
  expect_type(s$summary$`Prob Dir`,"double")
  
  expect_type(s$pairwise$Comparison,"character")
  expect_type(s$pairwise$`Post Mean`,"double")
  expect_type(s$pairwise$Lower,"double")
  expect_type(s$pairwise$Upper,"double")
  expect_type(unlist(s$pairwise[,5]),"double")
  expect_type(s$pairwise$EPR,"double")
  expect_type(s$pairwise$`EPR Lower`,"double")
  expect_type(s$pairwise$`EPR Upper`,"double")
  
  
  ## Make sure coef.aov_b works
  expect_type(coef(fita), "double")
  
  ## Make sure credint works
  expect_true(is.matrix(credint(fita)))
  expect_true(is.matrix(credint(fita,which = "pair")))
  
  ## Make sure vcov works
  expect_true(is.matrix(vcov(fita)))
  
  # Test default hyperparameters and mc_error
  expect_no_error(
    fitb <-
      aov_b(outcome ~ x1,
            test_data,
            mc_error = 0.01)
  )
  expect_equal(fitb$hyperparameters,
               list(mu = mean(test_data$outcome),
                    nu = 0.001,
                    a = 0.001,
                    b = 0.001))
  
  # Make sure prediction function works
  expect_no_error(predict(fita))
  expect_no_error(predict(fita,
                          newdata = fita$data[1,]))
  expect_gt(predict(fita,CI_level = 0.8)$CI_lower[1],
            predict(fita,CI_level = 0.9)$CI_lower[1])
  expect_gt(predict(fita,PI_level = 0.8)$PI_lower[1],
            predict(fita,PI_level = 0.9)$PI_lower[1])
  
  # Make sure information criteria work
  expect_type(AIC(fita),"double")
  expect_type(BIC(fita),"double")
  expect_type(DIC(fita),"double")
  expect_type(WAIC(fita),"double")
  
  
  # Make sure contrasts work
  ## One contrast
  expect_no_error(
    fitc <-
      aov_b(outcome ~ x1,
            test_data,
            mc_error = 0.01,
            contrasts = c(-1/3,-1/3,-1/3,1/2,1/2))
  )
  expect_named(fitc$contrasts,
               c("L","summary"))
  expect_s3_class(fitc$contrasts$summary,c("tbl_df", "tbl", "data.frame"))
  expect_type(fitc$contrasts$summary$`Post Mean`,"double")
  expect_type(fitc$contrasts$summary$Lower,"double")
  expect_type(fitc$contrasts$summary$Upper,"double")
  
  ## Multiple contrasts
  expect_no_error(
    fitd <-
      aov_b(outcome ~ x1,
            test_data,
            mc_error = 0.01,
            contrasts = rbind(c(-1/3,-1/3,-1/3,1/2,1/2),
                              c(-1/3,-1/3,-1/3,1,0)))
  )
  expect_named(fitd$contrasts,
               c("L","summary"))
  expect_s3_class(fitd$contrasts$summary,c("tbl_df", "tbl", "data.frame"))
  expect_type(fitd$contrasts$summary$`Post Mean`,"double")
  expect_type(fitd$contrasts$summary$Lower,"double")
  expect_type(fitd$contrasts$summary$Upper,"double")
  
  ## Invalid contrasts should throw an error
  expect_error(
    aov_b(outcome ~ x1,
          test_data,
          mc_error = 0.01,
          contrasts = c(-1,-1,-1,1,1))
  )
  
  
  
  # Make sure plotting function works
  expect_s3_class(plot(fita,
                       type = "diagnostics"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       combine_pi_ci = TRUE),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "pr"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "pr",
                       PI_level = 0.8),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "cr"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "cr",
                       CI_level = 0.999),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_s3_class(plot(fita),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  
  
  # Test if response transformation works
  test_data$e_outcome = exp(test_data$outcome)
  
  ## Test aov_b fit
  expect_no_error(
    fite <-
      aov_b(log(e_outcome) ~ x1,
            test_data,
            prior_mean_mu = 2,
            prior_mean_nu = 0.5,
            prior_var_shape = 0.01,
            prior_var_rate = 0.01)
  )
  expect_equal(fita$summary,
               fite$summary)
  
  ## Make sure prediction function works
  expect_no_error(
    fite_preds <- predict(fite)
  )
  expect_no_error(predict(fite,
                          newdata = fite$data[1,]))
  expect_gt(predict(fite,
                     newdata = fite$data[1,],
                     CI_level = 0.8)$CI_lower[1],
             predict(fite,
                     newdata = fite$data[1,],
                     CI_level = 0.9)$CI_lower[1])
  expect_gt(predict(fite,
                     newdata = fite$data[1,],
                     PI_level = 0.8)$PI_lower[1],
             predict(fite,
                     newdata = fite$data[1,],
                     PI_level = 0.9)$PI_lower[1])
  
  ## Test plot
  expect_s3_class(plot(fite),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  
  # Make sure parallelization works
  plan(multisession,workers = 5)
  expect_no_error(
    aov_b(outcome ~ x1,
          test_data,
          prior_mean_mu = 2,
          prior_mean_nu = 0.5,
          prior_var_shape = 0.01,
          prior_var_rate = 0.01)
  )
  plan(sequential)
  
  rm(list=ls())
})



# Proper, homoscedastic ---------------------------------------------------

test_that("Proper prior and homoscedastic model works", {
  
  # Create data
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rep(letters[1:5],N/5))
  test_data$outcome = 
    rnorm(N,-1 + 2 * (test_data$x1 %in% c("d","e")) )
  
  # No errors upon fitting
  expect_no_error(
    fita <-
      aov_b(outcome ~ x1,
            test_data,
            heteroscedastic = FALSE,
            prior_mean_mu = 2,
            prior_mean_nu = 0.5,
            prior_var_shape = 0.01,
            prior_var_rate = 0.01)
  )
  
  # Make sure print works
  expect_no_error(fita)
  
  # Make sure summary.aov_b works
  expect_no_error(
    s <- 
      summary(fita)
  )
  ## Check output format
  expect_type(s,"list")
  expect_s3_class(s$summary,c("tbl_df", "tbl", "data.frame"))
  
  expect_named(s$summary,c("Variable","Post Mean","Lower","Upper","Prob Dir"))
  expect_type(s$summary$Variable,"character")
  expect_type(s$summary$`Post Mean`,"double")
  expect_type(s$summary$Lower,"double")
  expect_type(s$summary$Upper,"double")
  expect_type(s$summary$`Prob Dir`,"double")
  
  expect_type(s$pairwise$Comparison,"character")
  expect_type(s$pairwise$`Post Mean`,"double")
  expect_type(s$pairwise$Lower,"double")
  expect_type(s$pairwise$Upper,"double")
  expect_type(unlist(s$pairwise[,5]),"double")
  expect_type(s$pairwise$EPR,"double")
  expect_type(s$pairwise$`EPR Lower`,"double")
  expect_type(s$pairwise$`EPR Upper`,"double")
  
  
  ## Make sure coef.aov_b works
  expect_type(coef(fita), "double")
  
  ## Make sure credint works
  expect_true(is.matrix(credint(fita)))
  expect_true(is.matrix(credint(fita,which = "pair")))
  
  ## Make sure vcov works
  expect_true(is.matrix(vcov(fita)))
  
  # Test default hyperparameters
  expect_no_error(
    fitb <-
      aov_b(outcome ~ x1,
            test_data,
            heteroscedastic = FALSE,
            mc_error = 0.01)
  )
  expect_equal(fitb$hyperparameters,
               list(mu = mean(test_data$outcome),
                    nu = 0.001,
                    a = 0.001,
                    b = 0.001))
  
  # Make sure prediction function works
  expect_no_error(predict(fita))
  expect_no_error(predict(fita,
                          newdata = fita$data[1,]))
  expect_gt(predict(fita,CI_level = 0.8)$CI_lower[1],
            predict(fita,CI_level = 0.9)$CI_lower[1])
  expect_gt(predict(fita,PI_level = 0.8)$PI_lower[1],
            predict(fita,PI_level = 0.9)$PI_lower[1])
  
  # Make sure information criteria work
  expect_type(AIC(fita),"double")
  expect_type(BIC(fita),"double")
  expect_type(DIC(fita),"double")
  expect_type(WAIC(fita),"double")
  
  
  # Make sure contrasts work
  ## One contrast
  expect_no_error(
    fitc <-
      aov_b(outcome ~ x1,
            test_data,
            heteroscedastic = FALSE,
            mc_error = 0.01,
            contrasts = c(-1/3,-1/3,-1/3,1/2,1/2))
  )
  expect_named(fitc$contrasts,
               c("L","summary"))
  expect_s3_class(fitc$contrasts$summary,c("tbl_df", "tbl", "data.frame"))
  expect_type(fitc$contrasts$summary$`Post Mean`,"double")
  expect_type(fitc$contrasts$summary$Lower,"double")
  expect_type(fitc$contrasts$summary$Upper,"double")
  
  ## Multiple contrasts
  expect_no_error(
    fitd <-
      aov_b(outcome ~ x1,
            test_data,
            heteroscedastic = FALSE,
            mc_error = 0.01,
            contrasts = rbind(c(-1/3,-1/3,-1/3,1/2,1/2),
                              c(-1/3,-1/3,-1/3,1,0)))
  )
  expect_named(fitd$contrasts,
               c("L","summary"))
  expect_s3_class(fitd$contrasts$summary,c("tbl_df", "tbl", "data.frame"))
  expect_type(fitd$contrasts$summary$`Post Mean`,"double")
  expect_type(fitd$contrasts$summary$Lower,"double")
  expect_type(fitd$contrasts$summary$Upper,"double")
  
  ## Invalid contrasts should throw an error
  expect_error(
    aov_b(outcome ~ x1,
          test_data,
          heteroscedastic = FALSE,
          mc_error = 0.01,
          contrasts = c(-1,-1,-1,1,1))
  )
  
  
  # Make sure plotting function works
  expect_s3_class(plot(fita,
                       type = "diagnostics"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       combine_pi_ci = TRUE),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "pr"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "pr",
                       PI_level = 0.8),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "cr"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "cr",
                       CI_level = 0.999),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_s3_class(plot(fita),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  # Make sure parallelization works
  plan(multisession,workers = 5)
  expect_no_error(
    aov_b(outcome ~ x1,
          test_data,
          prior_mean_mu = 2,
          prior_mean_nu = 0.5,
          prior_var_shape = 0.01,
          prior_var_rate = 0.01)
  )
  plan(sequential)
  
})





# Improper prior, heteroscedastic -----------------------------------------


test_that("Imroper prior and heteroscedastic model works", {
  
  # Create data
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rep(letters[1:5],N/5))
  test_data$outcome = 
    rnorm(N,-1 + 2 * (test_data$x1 %in% c("d","e")) )
  
  # No errors upon fitting
  expect_no_error(
    fita <-
      aov_b(outcome ~ x1,
            test_data,
            improper = TRUE)
  )
  
  # Make sure print works
  expect_no_error(fita)
  
  # Make sure summary.aov_b works
  expect_no_error(
    s <- 
      summary(fita)
  )
  ## Check output format
  expect_type(s,"list")
  expect_s3_class(s$summary,c("tbl_df", "tbl", "data.frame"))
  
  expect_named(s$summary,c("Variable","Post Mean","Lower","Upper","Prob Dir"))
  expect_type(s$summary$Variable,"character")
  expect_type(s$summary$`Post Mean`,"double")
  expect_type(s$summary$Lower,"double")
  expect_type(s$summary$Upper,"double")
  expect_type(s$summary$`Prob Dir`,"double")
  
  expect_type(s$pairwise$Comparison,"character")
  expect_type(s$pairwise$`Post Mean`,"double")
  expect_type(s$pairwise$Lower,"double")
  expect_type(s$pairwise$Upper,"double")
  expect_type(unlist(s$pairwise[,5]),"double")
  expect_type(s$pairwise$EPR,"double")
  expect_type(s$pairwise$`EPR Lower`,"double")
  expect_type(s$pairwise$`EPR Upper`,"double")
  
  
  ## Make sure coef.aov_b works
  expect_type(coef(fita), "double")
  
  ## Make sure credint works
  expect_true(is.matrix(credint(fita)))
  expect_true(is.matrix(credint(fita,which = "pair")))
  
  ## Make sure vcov works
  expect_true(is.matrix(vcov(fita)))
  
  # Test default hyperparameters DO NOT change anything
  expect_no_error(
    fitb <-
      aov_b(outcome ~ x1,
            test_data,
            improper = TRUE,
            prior_mean_mu = 200,
            prior_mean_nu = 0.5,
            prior_var_shape = 0.01,
            prior_var_rate = 0.01)
  )
  expect_equal(fita$summary,fitb$summary)
  
  # Make sure prediction function works
  expect_no_error(predict(fita))
  expect_no_error(predict(fita,
                          newdata = fita$data[1,]))
  expect_gt(predict(fita,CI_level = 0.8)$CI_lower[1],
            predict(fita,CI_level = 0.9)$CI_lower[1])
  expect_gt(predict(fita,PI_level = 0.8)$PI_lower[1],
            predict(fita,PI_level = 0.9)$PI_lower[1])
  
  # Make sure information criteria work
  expect_type(AIC(fita),"double")
  expect_type(BIC(fita),"double")
  expect_type(DIC(fita),"double")
  expect_type(WAIC(fita),"double")
  
  # Test contrasts
  ## One contrast
  expect_no_error(
    fitc <-
      aov_b(outcome ~ x1,
            test_data,
            improper = TRUE,
            mc_error = 0.01,
            contrasts = c(-1/3,-1/3,-1/3,1/2,1/2))
  )
  expect_named(fitc$contrasts,
               c("L","summary"))
  expect_s3_class(fitc$contrasts$summary,c("tbl_df", "tbl", "data.frame"))
  expect_type(fitc$contrasts$summary$`Post Mean`,"double")
  expect_type(fitc$contrasts$summary$Lower,"double")
  expect_type(fitc$contrasts$summary$Upper,"double")
  
  ## Multiple contrasts
  expect_no_error(
    fitd <-
      aov_b(outcome ~ x1,
            test_data,
            improper = TRUE,
            mc_error = 0.01,
            contrasts = rbind(c(-1/3,-1/3,-1/3,1/2,1/2),
                              c(-1/3,-1/3,-1/3,1,0)))
  )
  expect_named(fitd$contrasts,
               c("L","summary"))
  expect_s3_class(fitd$contrasts$summary,c("tbl_df", "tbl", "data.frame"))
  expect_type(fitd$contrasts$summary$`Post Mean`,"double")
  expect_type(fitd$contrasts$summary$Lower,"double")
  expect_type(fitd$contrasts$summary$Upper,"double")
  
  ## Invalid contrasts should throw an error
  expect_error(
    aov_b(outcome ~ x1,
          test_data,
          improper = TRUE,
          mc_error = 0.01,
          contrasts = c(-1,-1,-1,1,1))
  )
  
  
  # Make sure plotting function works
  expect_s3_class(plot(fita,
                       type = "diagnostics"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       combine_pi_ci = TRUE),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "pr"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "pr",
                       PI_level = 0.8),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "cr"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "cr",
                       CI_level = 0.999),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_s3_class(plot(fita),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  # Make sure parallelization works
  plan(multisession,workers = 5)
  expect_no_error(
    aov_b(outcome ~ x1,
          test_data,
          improper = TRUE)
  )
  plan(sequential)
  
  rm(list=ls())
})



# Improper prior, homoscedastic -----------------------------------------


test_that("Imroper prior and homoscedastic model works", {
  
  # Create data
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rep(letters[1:5],N/5))
  test_data$outcome = 
    rnorm(N,-1 + 2 * (test_data$x1 %in% c("d","e")) )
  
  # No errors upon fitting
  expect_no_error(
    fita <-
      aov_b(outcome ~ x1,
            test_data,
            improper = TRUE,
            heteroscedastic = FALSE)
  )
  
  # Make sure print works
  expect_no_error(fita)
  
  # Make sure summary.aov_b works
  expect_no_error(
    s <- 
      summary(fita)
  )
  ## Check output format
  expect_type(s,"list")
  expect_s3_class(s$summary,c("tbl_df", "tbl", "data.frame"))
  
  expect_named(s$summary,c("Variable","Post Mean","Lower","Upper","Prob Dir"))
  expect_type(s$summary$Variable,"character")
  expect_type(s$summary$`Post Mean`,"double")
  expect_type(s$summary$Lower,"double")
  expect_type(s$summary$Upper,"double")
  expect_type(s$summary$`Prob Dir`,"double")
  
  expect_type(s$pairwise$Comparison,"character")
  expect_type(s$pairwise$`Post Mean`,"double")
  expect_type(s$pairwise$Lower,"double")
  expect_type(s$pairwise$Upper,"double")
  expect_type(unlist(s$pairwise[,5]),"double")
  expect_type(s$pairwise$EPR,"double")
  expect_type(s$pairwise$`EPR Lower`,"double")
  expect_type(s$pairwise$`EPR Upper`,"double")
  
  
  ## Make sure coef.aov_b works
  expect_type(coef(fita), "double")
  
  ## Make sure credint works
  expect_true(is.matrix(credint(fita)))
  expect_true(is.matrix(credint(fita,which = "pair")))
  
  ## Make sure vcov works
  expect_true(is.matrix(vcov(fita)))
  
  # Test default hyperparameters DO NOT change anything
  expect_no_error(
    fitb <-
      aov_b(outcome ~ x1,
            test_data,
            improper = TRUE,
            prior_mean_mu = 200,
            prior_mean_nu = 0.5,
            prior_var_shape = 0.01,
            prior_var_rate = 0.01,
            heteroscedastic = FALSE)
  )
  expect_equal(fita$summary,fitb$summary)
  
  # Make sure prediction function works
  expect_no_error(predict(fita))
  expect_no_error(predict(fita,
                          newdata = fita$data[1,]))
  expect_gt(predict(fita,CI_level = 0.8)$CI_lower[1],
            predict(fita,CI_level = 0.9)$CI_lower[1])
  expect_gt(predict(fita,PI_level = 0.8)$PI_lower[1],
            predict(fita,PI_level = 0.9)$PI_lower[1])
  
  # Make sure information criteria work
  expect_type(AIC(fita),"double")
  expect_type(BIC(fita),"double")
  expect_type(DIC(fita),"double")
  expect_type(WAIC(fita),"double")
  
  # Test contrasts
  ## One contrast
  expect_no_error(
    fitc <-
      aov_b(outcome ~ x1,
            test_data,
            improper = TRUE,
            heteroscedastic = FALSE,
            mc_error = 0.01,
            contrasts = c(-1/3,-1/3,-1/3,1/2,1/2))
  )
  expect_named(fitc$contrasts,
               c("L","summary"))
  expect_s3_class(fitc$contrasts$summary,c("tbl_df", "tbl", "data.frame"))
  expect_type(fitc$contrasts$summary$`Post Mean`,"double")
  expect_type(fitc$contrasts$summary$Lower,"double")
  expect_type(fitc$contrasts$summary$Upper,"double")
  
  ## Multiple contrasts
  expect_no_error(
    fitd <-
      aov_b(outcome ~ x1,
            test_data,
            improper = TRUE,
            heteroscedastic = FALSE,
            mc_error = 0.01,
            contrasts = rbind(c(-1/3,-1/3,-1/3,1/2,1/2),
                              c(-1/3,-1/3,-1/3,1,0)))
  )
  expect_named(fitd$contrasts,
               c("L","summary"))
  expect_s3_class(fitd$contrasts$summary,c("tbl_df", "tbl", "data.frame"))
  expect_type(fitd$contrasts$summary$`Post Mean`,"double")
  expect_type(fitd$contrasts$summary$Lower,"double")
  expect_type(fitd$contrasts$summary$Upper,"double")
  
  ## Invalid contrasts should throw an error
  expect_error(
    aov_b(outcome ~ x1,
          test_data,
          improper = TRUE,
          heteroscedastic = FALSE,
          mc_error = 0.01,
          contrasts = c(-1,-1,-1,1,1))
  )
  
  # Make sure plotting function works
  expect_s3_class(plot(fita,
                       type = "diagnostics"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = c("cr","pr"),
                       combine_pi_ci = TRUE),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "pr"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "pr",
                       PI_level = 0.8),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "cr"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       type = "cr",
                       CI_level = 0.999),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_s3_class(plot(fita),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  # Make sure parallelization works
  plan(multisession,workers = 5)
  expect_no_error(
    aov_b(outcome ~ x1,
          test_data,
          improper = TRUE,
          heteroscedastic = FALSE)
  )
  plan(sequential)
  
  rm(list=ls())
})


