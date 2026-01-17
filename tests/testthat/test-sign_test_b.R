
test_that("Test sign_test_b",{
  
  # Test input
  expect_no_error(
    fita <-
      sign_test_b(x = rnorm(50))
  )
  expect_no_error(
    sign_test_b(x = rnorm(50,1),
                y = rnorm(50,0))
  )
  
  # Test output
  expect_type(fita$posterior_mean,"double")
  expect_type(fita$CI,"double")
  expect_type(fita$Pr_less_than_p,"double")
  expect_type(fita$ROPE,"double")
  expect_type(fita$ROPE_bounds,"double")
  expect_type(fita$posterior_parameters,"double")
  expect_s3_class(fita$prop_plot,
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
  expect_equal(fitb[c(1:5,7)],
               fitc[c(1:5,7)])
  expect_true(!isTRUE(all.equal(fitb[c(1:5,7)],
                                fitd[c(1:5,7)])))
  
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
  expect_lt(fite$ROPE,
            fitf$ROPE)
  expect_equal(fite[c(1:5,7)],
               fitg[c(1:5,7)])
  
  # Test changing reference probability
  expect_no_error(
    fith <-
      sign_test_b(x = x,
                  y = y,
                  p0 = 0.7)
  )
  expect_lt(fite$Pr_less_than_p,
            fith$Pr_less_than_p)
  expect_error(
    sign_test_b(x = rnorm(50,1),
                p0 = 0.71,
                ROPE = 0.3)
  )
  
  
})