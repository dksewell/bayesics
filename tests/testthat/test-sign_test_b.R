
test_that("Test sign_test_b",{
  
  # Test input
  expect_no_error(
    fita <-
      sign_test_b(x = rnorm(50))
  )
  expect_no_error(fita)
  expect_s3_class(plot(fita),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_no_error(
    fitb <-
      sign_test_b(x = rnorm(50,1),
                  y = rnorm(50,0))
  )
  expect_no_error(fitb)
  expect_s3_class(plot(fitb),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  
  # Test prior
  set.seed(2025)
  x = rnorm(50,0)
  y = rnorm(50,1)
  expect_no_error(
    fitb <-
      sign_test_b(x = x,
                  y = y,
                  prior = "uniform")
  )
  expect_no_error(
    fitc <-
      sign_test_b(x = x,
                  y = y,
                  prior_shapes = c(1,1))
  )
  expect_no_error(
    fitd <-
      sign_test_b(x = x,
                  y = y,
                  prior_shapes = c(2,2))
  )
  expect_equal(fitb$results,
               fitc$results)
  expect_true(!isTRUE(all.equal(fitb$results,
                                fitd$results)))
  
  # Test ROPE
  expect_no_error(
    fite <-
      sign_test_b(x = x,
                  y = y,
                  ROPE = 0.1)
  )
  expect_no_error(
    fitf <-
      sign_test_b(x = x,
                  y = y,
                  ROPE = 0.15)
  )
  expect_no_error(
    fitg <-
      sign_test_b(x = x,
                  y = y,
                  ROPE = c(0.4,0.6))
  )
  expect_gt(fite$results$ROPE_lower_bound,
            fitf$results$ROPE_lower_bound)
  expect_equal(fite$results,
               fitg$results)
  
  # Test changing reference probability
  expect_no_error(
    fith <-
      sign_test_b(x = x,
                  y = y,
                  p0 = 0.7)
  )
  expect_lt(fite$pdir$pdir,
            fith$pdir$pdir)
  expect_error(
    sign_test_b(x = rnorm(50,1),
                p0 = 0.71,
                ROPE = 0.3)
  )
  
  
})