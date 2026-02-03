

test_that("Test poisson_test_b for a single population",{
  
  # No offset
  expect_no_error(
    fita <- 
      poisson_test_b(x = 12)
  )
  expect_equal(fita$x,12)
  expect_equal(fita$offset,1)
  expect_type(fita$posterior_mean,"double")
  expect_type(fita$posterior_parameters,"double")
  expect_type(fita$CI,"double")
  expect_s3_class(fita$rate_plot,
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  # with offset
  expect_no_error(
    fitb <- 
      poisson_test_b(x = 12,
                     offset = 2)
  )
  expect_equal(fitb$x,12)
  expect_equal(fitb$offset,2)
  expect_type(fitb$posterior_mean,"double")
  expect_type(fitb$posterior_parameters,"double")
  expect_type(fitb$CI,"double")
  expect_s3_class(fitb$rate_plot,
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  # With reference value
  expect_no_error(
    fitc <- 
      poisson_test_b(x = 12,
                     offset = 2,
                     r = 10)
  )
  expect_type(fitc$Pr_less_than_r,"double")
  
  # With different prior
  expect_no_error(
    poisson_test_b(x = 12,
                   offset = 2,
                   r = 11,
                   prior = "flat")
  )
  expect_no_error(
    poisson_test_b(x = 12,
                   offset = 2,
                   r = 11,
                   prior_shape_rate = c(1,1))
  )
  
   
})



test_that("Test poisson_test_b for two populations",{
  
  # No offset
  expect_no_error(
    fita <- 
      poisson_test_b(x = c(12,20))
  )
  expect_equal(fita$x,c(12,20))
  expect_equal(fita$offset,c(1,1))
  expect_type(fita$posterior_mean_pop1,"double")
  expect_type(fita$posterior_mean_pop2,"double")
  expect_type(fita$CI_pop1,"double")
  expect_type(fita$CI_pop2,"double")
  expect_type(fita$CI_lambda1_over_lambda2,"double")
  expect_type(fita$Pr_rateratio_in_ROPE,"double")
  expect_type(fita$posterior_parameters$population_1,"double")
  expect_type(fita$posterior_parameters$population_2,"double")
  expect_s3_class(fita$rate_plot,
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  # with offset
  expect_no_error(
    fitb <- 
      poisson_test_b(x = c(12,20),
                     offset = c(10,9))
  )
  expect_equal(fitb$x,c(12,20))
  expect_equal(fitb$offset,c(10,9))
  expect_type(fitb$posterior_mean_pop1,"double")
  expect_type(fitb$posterior_mean_pop2,"double")
  expect_type(fitb$CI_pop1,"double")
  expect_type(fitb$CI_pop2,"double")
  expect_type(fitb$CI_lambda1_over_lambda2,"double")
  expect_type(fitb$Pr_rateratio_in_ROPE,"double")
  expect_type(fitb$posterior_parameters$population_1,"double")
  expect_type(fitb$posterior_parameters$population_2,"double")
  expect_s3_class(fitb$rate_plot,
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  # reference value shouldn't do anything here
  expect_no_error(
    fitc <- 
      poisson_test_b(x = c(12,20),
                     offset = c(10,9),
                     r = 10)
  )
  expect_equal(fitb[c(1:9,11)],
               fitc[c(1:9,11)])
  
  # With different prior
  expect_no_error(
    poisson_test_b(x = c(12,20),
                   offset = c(10,9),
                   prior = "flat")
  )
  expect_no_error(
    poisson_test_b(x = c(12,20),
                   offset = c(10,9),
                   prior_shape_rate = c(1,1))
  )
  
  
})

