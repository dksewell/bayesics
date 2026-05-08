

test_that("Test poisson_test_b for a single population",{
  
  # No offset
  expect_no_error(
    fita <- 
      poisson_test_b(x = 12)
  )
  expect_no_error(fita)
  expect_s3_class(plot(fita),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  # with offset
  expect_no_error(
    fitb <- 
      poisson_test_b(x = 12,
                     offset = 2)
  )
  expect_no_error(fitb)
  expect_s3_class(plot(fitb),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  # With reference value
  expect_no_error(
    fitc <- 
      poisson_test_b(x = 12,
                     offset = 2,
                     r = 10)
  )
  expect_no_error(fitc)
  expect_s3_class(plot(fitc),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
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
  expect_no_error(fita)
  expect_s3_class(plot(fita),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  
  # with offset
  expect_no_error(
    fitb <- 
      poisson_test_b(x = c(12,20),
                     offset = c(10,9))
  )
  expect_no_error(fitb)
  expect_s3_class(plot(fitb),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  # reference value shouldn't do anything here
  expect_no_error(
    fitc <- 
      poisson_test_b(x = c(12,20),
                     offset = c(10,9),
                     r = 10)
  )
  expect_equal(fitb$results,
               fitc$results)
  
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

