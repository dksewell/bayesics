

test_that("Test prop_test_b for a single population",{
  
  # Test input methods
  expect_no_error(
    fita <-
      prop_test_b(14,
                  19)
  )
  expect_no_error(fita)
  expect_s3_class(plot(fita),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_no_error(
    fitb <-
      prop_test_b(14,
                  n_total = 14 + 19)
  )
  expect_equal(fita$results,
               fitb$results)
  
  # Test probability of comparison
  expect_no_error(
    fitc <-
      prop_test_b(14,
                  19,
                  p = 0.45)
  )
  expect_no_error(fitc)
  expect_s3_class(plot(fitc),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
})



test_that("Test prop_test_b for two populations",{
  
  # Test input methods
  expect_no_error(
    fita <-
      prop_test_b(c(14,22),
                  c(19,45))
  )
  expect_no_error(
    fitb <-
      prop_test_b(c(14,22),
                  n_total = c(14,22) + c(19,45))
  )
  expect_equal(fita$results,
               fitb$results)
  
  # Test output
  expect_no_error(fita)
  expect_s3_class(plot(fita),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  # Test probability of comparison doesn't do anything
  expect_no_error(
    fitc <-
      prop_test_b(c(14,22),
                  c(19,45),
                  p = 0.45)
  )
  expect_equal(fita$results,
               fitc$results)
  
})