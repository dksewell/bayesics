
test_that("Test single population semi-parametric survival curve fitting",{
  set.seed(2025)
  N = 300
  test_data = 
    data.frame(outcome = 
                 rweibull(N,2,5))
  test_data$observed = 
    ifelse(test_data$outcome >= 7, 0, 1)
  test_data$outcome =
    ifelse(dplyr::near(test_data$observed,1), test_data$outcome, 7)
  
  # Check fitting function
  expect_no_error(
    fita <-
      survfit_b(Surv(test_data$outcome,
                     test_data$observed) ~ 1)
  )
  expect_error(
    fita <-
      survfit_b(Surv(test_data$outcome,
                     test_data$observed))
  )
  expect_no_error(
    fitb <-
      survfit_b(Surv(outcome,
                     observed) ~ 1,
                data = test_data)
  )
  
  # Check print
  expect_no_error(
    fita
  )
  
  # Check plotting
  expect_no_error(
    survplot <- 
      plot(fita)
  )
  expect_type(survplot,"list")
  expect_s3_class(survplot$plot,
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(survplot$data,c("tbl_df", "tbl", "data.frame"))
  
  
})


test_that("Test multiple populations via semi-parametric survival curve fitting",{
  set.seed(2025)
  N = 300
  test_data = 
    data.frame(outcome = 
                 c(rweibull(2*N/3,2,5),
                   rweibull(N/3,2,10)),
               x1 = rep(letters[1:3],each = N/3))
  test_data$observed = 
    ifelse(test_data$outcome >= 9, 0, 1)
  test_data$outcome =
    ifelse(dplyr::near(test_data$observed,1), test_data$outcome, 9)
  
  # Check fitting function
  expect_no_error(
    fita <-
      survfit_b(Surv(outcome,
                     observed) ~ x1,
                data = test_data)
  )
  expect_no_error(
    fitb <-
      survfit_b(Surv(test_data$outcome,
                     test_data$observed) ~ test_data$x1)
  )
  
  # Check print
  expect_no_error(
    fita
  )
  
  # Check plotting
  expect_no_error(
    survplot <- 
      plot(fita)
  )
  expect_type(survplot,"list")
  expect_s3_class(survplot$plot,
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_s3_class(survplot$data,c("tbl_df", "tbl", "data.frame"))
  
})
