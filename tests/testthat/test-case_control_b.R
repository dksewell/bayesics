
test_that("Simple case-control analysis", {
  
  # Test input
  expect_no_error(
    fita <- case_control_b(matrix(c(8,47,1,26),2,2))
  )
  expect_no_error(fita)
  expect_s3_class(plot(fita),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_no_error(
    fitb <- case_control_b(c(8,47),
                           c(1,26))
  )
  expect_no_error(fitb)
  expect_s3_class(plot(fitb),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_no_error(
    fitc <- case_control_b(x = matrix(c(8,47,1,26),2,2))
  )
  expect_no_error(fitc)
  expect_s3_class(plot(fitc),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  # Test large sample
  expect_no_error(
    fitd <- case_control_b(x = 5 + matrix(c(8,47,1,26),2,2))
  )
  expect_no_error(fitd)
  expect_s3_class(plot(fitd),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  # Test ROPE
  expect_no_error(
    fite <- case_control_b(x = matrix(c(8,47,1,26),2,2),
                           ROPE = 1.05)
  )
  expect_no_error(fite)
  expect_s3_class(plot(fite),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_no_error(
    fitf <- case_control_b(x = 5 + matrix(c(8,47,1,26),2,2),
                           ROPE = 1.05)
  )
  expect_no_error(fitf)
  expect_s3_class(plot(fitf),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  # Test prior
  expect_no_error(
    fitg <- case_control_b(x = 5 + matrix(c(8,47,1,26),2,2),
                           ROPE = 1.05)
  )
  expect_no_error(fitg)
  expect_s3_class(plot(fitg),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_no_error(
    fith <- case_control_b(x = 5 + matrix(c(8,47,1,26),2,2),
                           ROPE = 1.05,
                           prior_mean = 10)
  )
  expect_no_error(fith)
  expect_s3_class(plot(fith),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_no_error(
    fiti <- case_control_b(x = 5 + matrix(c(8,47,1,26),2,2),
                           ROPE = 1.05,
                           prior_sd = 0.01)
  )
  expect_no_error(fiti)
  expect_s3_class(plot(fiti),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
})