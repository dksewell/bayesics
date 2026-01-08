
test_that("Bayesian model averaging", {
  
  # Create data
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5],
               x4 = rnorm(N),
               x5 = rnorm(N),
               x6 = rnorm(N),
               x7 = rnorm(N),
               x8 = rnorm(N),
               x9 = rnorm(N),
               x10 = rnorm(N))
  test_data$outcome = 
    rnorm(N,-1 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) )
  
  
  expect_no_error(
    fita <-
      bma_inference(outcome ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10,
                    test_data,
                    user.int = FALSE)
  )
  
  # Make sure print works
  expect_no_error(fita)
  
  # Make sure coef works
  expect_type(coef(fita), "double")
  
  # Make sure credint works
  expect_true(is.matrix(credint(fita)))
  
  # Make sure summary.aov_b works
  expect_no_error(
    s <- 
      summary(fita)
  )
  ## Check output format
  expect_s3_class(s,c("tbl_df", "tbl", "data.frame"))
  expect_named(s,
               c("Variable","Post Mean","Lower","Upper","ROPE","Prob Dir"))
  expect_type(s$Variable,"character")
  expect_type(s$`Post Mean`,"double")
  expect_type(s$Lower,"double")
  expect_type(s$Upper,"double")
  expect_type(s$ROPE,"double")
  expect_type(s$`Prob Dir`,"double")
  
  # Make sure predict works
  expect_no_error(predict(fita))
  expect_no_error(predict(fita,
                          newdata = fita$data[1,]))
  expect_gt(predict(fita,CI_level = 0.8)$newdata$CI_lower[1],
            predict(fita,CI_level = 0.9)$newdata$CI_lower[1])
  expect_gt(predict(fita,PI_level = 0.8)$newdata$PI_lower[1],
            predict(fita,PI_level = 0.9)$newdata$PI_lower[1])
  
  # Test plot
  expect_s3_class(plot(fita,
                       "dx"),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       "pdp",
                       variable = c("x1","x2","x3")),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       "ci",
                       variable = c("x1","x2","x3")),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       "pi",
                       variable = c("x1","x2","x3")),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       c("ci","pi"),
                       variable = c("x1","x2","x3")),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(plot(fita,
                       c("ci","pi"),
                       combine_pi_ci = FALSE,
                       variable = c("x1","x2","x3")),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  
  
  # # Make sure parallelization works. 
  # plan(multisession, workers = 5)
  # expect_no_error(
  #   fitb <-
  #     bma_inference(outcome ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10,
  #                   test_data,
  #                   mc_draws = 5e2,
  #                   user.int = FALSE)
  # )
  # plan(sequential)
  
  # Check if splines and factors work
  library(splines)
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5],
               x4 = rnorm(N),
               x5 = rnorm(N),
               x6 = rnorm(N),
               x7 = rnorm(N),
               x8 = rnorm(N),
               x9 = rnorm(N),
               x10 = rnorm(N))
  test_data$outcome = 
    rnorm(N,-1 + test_data$x1 + test_data$x1^2 + 2 * (test_data$x3 %in% c("d","e")) )
  expect_no_error(
    fitc <- 
      bma_inference(outcome ~ ns(x1,df = 5) + x2 + x3,
                    data = test_data,
                    mc_draws = 5e3,
                    user.int = FALSE)
  )
  
  rm(list = ls())
  
})