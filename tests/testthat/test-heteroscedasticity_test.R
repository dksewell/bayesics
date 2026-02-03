
test_that("Testing ANOVA heteroscedastity BF",{
  
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
  
  expect_no_error(
    htest1 <- 
      heteroscedasticity_test(hetero_model,
                              homo_model)
  )
  expect_s3_class(htest1,c("tbl_df", "tbl", "data.frame"))
  expect_type(htest1$BF,"double")
  expect_type(htest1$Interpretation,"character")
  
  
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
  expect_no_error(
    htest2 <- 
      heteroscedasticity_test(hetero_model,
                              homo_model)
  )
  
  # Make sure homoscedastic data gives larger BF than heteroscedastic data
  expect_gte(htest1$BF,
             htest2$BF)
  
})