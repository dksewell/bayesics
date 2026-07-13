do_full_testing = FALSE

# avoid the automatic warning from future
suppressWarnings({
  future.apply::future_sapply(1:2,sum)
})


test_that("Binary treatment, GLM response mediation",{
  
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(tr = rep(0:1,N/2),
               x1 = rnorm(N))
  test_data$m = 
    rnorm(N, 0.4 * test_data$tr - 0.25 * test_data$x1)
  test_data$outcome = 
    rpois(N,exp(-1 + 0.6 * test_data$tr + 1.5 * test_data$m + 0.25 * test_data$x1))
  
  m1 = 
    lm_b(m ~ tr + x1,
         data = test_data)
  m2 = 
    glm_b(outcome ~ m + tr + x1,
          data = test_data,
          family = poisson())
  expect_no_error(
    m3 <-
      mediate_b(m1,m2,
                treat = "tr",
                control_value = 0,
                treat_value = 1,
                n_draws = 500,
                mc_error = 0.05,
                ask_before_full_sampling = FALSE)
  )
  
  # Make sure print works
  expect_no_error(m3)
  
  # Check output
  expect_s3_class(m3$summary,c("tbl_df", "tbl", "data.frame"))
  expect_s3_class(m3$posterior_draws,c("tbl_df", "tbl", "data.frame"))
  expect_type(m3$treat_value,"double")
  expect_type(m3$control_value,"double")
  expect_type(m3$CI_level,"double")
  expect_type(m3$mc_error,"double")
  expect_s3_class(m3$model_m,c("lm_b","glm_b"))
  expect_s3_class(m3$model_y,c("lm_b","glm_b"))
  expect_identical(m3$summary$Estimand,
                   c("ACME (Control)",
                     "ACME (Treatment)",
                     "ADE (Control)",
                     "ADE (Treatment)",
                     "Total Effect",
                     "ACME (Average)",
                     "ADE (Average)",
                     "Prop. Mediated (Average)"))
  expect_type(m3$summary$Estimate,"double")
  expect_type(m3$summary$Lower,"double")
  expect_type(m3$summary$Upper,"double")
  expect_type(m3$summary$`Prob Dir`,"double")
  
  # Check summary 
  expect_no_error(
    s <- 
      summary(m3)
  )
  expect_silent(
    summary(m3,print_results=FALSE)
  )
  expect_identical(s$Estimand,
                   c("ACME (Control)",
                     "ACME (Treatment)",
                     "ADE (Control)",
                     "ADE (Treatment)",
                     "Total Effect",
                     "ACME (Average)",
                     "ADE (Average)",
                     "Prop. Mediated (Average)"))
  expect_type(s$Estimate,"double")
  expect_type(s$Lower,"double")
  expect_type(s$Upper,"double")
  expect_type(s$`Prob Dir`,"double")
  
  # Check plot
  if(do_full_testing){
    expect_s3_class(plot(m3),
                    c("patchwork","ggplot2::ggplot","ggplot",
                      "ggplot2::gg","S7_object","gg"))
    expect_s3_class(plot(m3,type = "diagnostics"),
                    c("patchwork","ggplot2::ggplot","ggplot",
                      "ggplot2::gg","S7_object","gg"))
    expect_s3_class(plot(m3,type = "acme"),
                    c("patchwork","ggplot2::ggplot","ggplot",
                      "ggplot2::gg","S7_object","gg"))
    expect_s3_class(plot(m3,type = "ade"),
                    c("patchwork","ggplot2::ggplot","ggplot",
                      "ggplot2::gg","S7_object","gg"))
    expect_s3_class(plot(m3,type = c("acme","ade")),
                    c("patchwork","ggplot2::ggplot","ggplot",
                      "ggplot2::gg","S7_object","gg"))
  }
  
})


test_that("Continuous treatment, GLM response mediation",{
  
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(tr = rnorm(N),
               x1 = rnorm(N))
  test_data$m = 
    rnorm(N, 0.4 * test_data$tr - 0.25 * test_data$x1)
  test_data$outcome = 
    rpois(N,exp(-1 + 0.6 * test_data$tr + 1.5 * test_data$m + 0.25 * test_data$x1))
  
  m1 = 
    lm_b(m ~ tr + x1,
         data = test_data)
  m2 = 
    glm_b(outcome ~ m + tr + x1,
          data = test_data,
          family = poisson())
  expect_no_error(
    m3 <-
      mediate_b(m1,m2,
                treat = "tr",
                control_value = -1,
                treat_value = 1,
                n_draws = 500,
                mc_error = 0.05,
                ask_before_full_sampling = FALSE)
  )
  
  # Check automatic selection of control and treat values
  expect_no_error(
    m4 <-
      mediate_b(m1,m2,
                treat = "tr",
                n_draws = 500,
                mc_error = 0.05,
                ask_before_full_sampling = FALSE)
  )
  expect_no_error(
    m5 <-
      mediate_b(m1,m2,
                treat = "tr",
                control_value = quantile(test_data$tr,0.2),
                treat_value = quantile(test_data$tr,0.8),
                n_draws = 500,
                mc_error = 0.05,
                ask_before_full_sampling = FALSE)
  )
  expect_equal(m4$summary,
               m5$summary)
  
  # Make sure print works
  expect_no_error(m3)
  
  # Check output
  expect_s3_class(m3$summary,c("tbl_df", "tbl", "data.frame"))
  expect_s3_class(m3$posterior_draws,c("tbl_df", "tbl", "data.frame"))
  expect_type(m3$treat_value,"double")
  expect_type(m3$control_value,"double")
  expect_type(m3$CI_level,"double")
  expect_type(m3$mc_error,"double")
  expect_s3_class(m3$model_m,c("lm_b","glm_b"))
  expect_s3_class(m3$model_y,c("lm_b","glm_b"))
  expect_identical(m3$summary$Estimand,
                   c("ACME (Control)",
                     "ACME (Treatment)",
                     "ADE (Control)",
                     "ADE (Treatment)",
                     "Total Effect",
                     "ACME (Average)",
                     "ADE (Average)",
                     "Prop. Mediated (Average)"))
  expect_type(m3$summary$Estimate,"double")
  expect_type(m3$summary$Lower,"double")
  expect_type(m3$summary$Upper,"double")
  expect_type(m3$summary$`Prob Dir`,"double")
  
  # Check summary 
  expect_no_error(
    s <- 
      summary(m3)
  )
  expect_silent(
    summary(m3,print_results=F)
  )
  expect_identical(s$Estimand,
                   c("ACME (Control)",
                     "ACME (Treatment)",
                     "ADE (Control)",
                     "ADE (Treatment)",
                     "Total Effect",
                     "ACME (Average)",
                     "ADE (Average)",
                     "Prop. Mediated (Average)"))
  expect_type(s$Estimate,"double")
  expect_type(s$Lower,"double")
  expect_type(s$Upper,"double")
  expect_type(s$`Prob Dir`,"double")
  
  # Check plot
  if(do_full_testing){
    expect_s3_class(plot(m3),
                    c("patchwork","ggplot2::ggplot","ggplot",
                      "ggplot2::gg","S7_object","gg"))
    expect_s3_class(plot(m3,type = "diagnostics"),
                    c("patchwork","ggplot2::ggplot","ggplot",
                      "ggplot2::gg","S7_object","gg"))
    expect_s3_class(plot(m3,type = "acme"),
                    c("patchwork","ggplot2::ggplot","ggplot",
                      "ggplot2::gg","S7_object","gg"))
    expect_s3_class(plot(m3,type = "ade"),
                    c("patchwork","ggplot2::ggplot","ggplot",
                      "ggplot2::gg","S7_object","gg"))
    expect_s3_class(plot(m3,type = c("acme","ade")),
                    c("patchwork","ggplot2::ggplot","ggplot",
                      "ggplot2::gg","S7_object","gg"))
  }
  
})


test_that("Simple case mediation",{
  
  set.seed(2025)
  N = 500
  test_data = 
    data.frame(tr = rnorm(N),
               x1 = rnorm(N))
  test_data$m = 
    rnorm(N, 0.4 * test_data$tr - 0.25 * test_data$x1)
  test_data$outcome = 
    rnorm(N,-1 + 0.6 * test_data$tr + 1.5 * test_data$m + 0.25 * test_data$x1)
  
  m1 = 
    lm_b(m ~ tr + x1,
         data = test_data)
  m2 = 
    lm_b(outcome ~ m + tr + x1,
         data = test_data)
  expect_no_error(
    m3 <-
      mediate_b(m1,m2,
                treat = "tr",
                control_value = -2,
                treat_value = 2,
                n_draws = 500,
                mc_error = 0.05,
                ask_before_full_sampling = FALSE)
  )
  
  # Make sure print works
  expect_no_error(m3)
  
  # Check output
  expect_s3_class(m3$summary,c("tbl_df", "tbl", "data.frame"))
  expect_s3_class(m3$posterior_draws,c("tbl_df", "tbl", "data.frame"))
  expect_type(m3$treat_value,"double")
  expect_type(m3$control_value,"double")
  expect_type(m3$CI_level,"double")
  expect_type(m3$mc_error,"double")
  expect_s3_class(m3$model_m,c("lm_b","glm_b"))
  expect_s3_class(m3$model_y,c("lm_b","glm_b"))
  expect_identical(m3$summary$Estimand,
                   c("ACME",
                     "ADE",
                     "Total Effect",
                     "Prop. Mediated"))
  expect_type(m3$summary$Estimate,"double")
  expect_type(m3$summary$Lower,"double")
  expect_type(m3$summary$Upper,"double")
  expect_type(m3$summary$`Prob Dir`,"double")
  
  # Check summary 
  expect_no_error(
    s <- 
      summary(m3)
  )
  expect_silent(
    summary(m3,print_results=FALSE)
  )
  expect_identical(s$Estimand,
                   c("ACME",
                     "ADE",
                     "Total Effect",
                     "Prop. Mediated"))
  expect_type(s$Estimate,"double")
  expect_type(s$Lower,"double")
  expect_type(s$Upper,"double")
  expect_type(s$`Prob Dir`,"double")
  
  # Check plot
  if(do_full_testing){
    expect_s3_class(plot(m3),
                    c("patchwork","ggplot2::ggplot","ggplot",
                      "ggplot2::gg","S7_object","gg"))
    expect_s3_class(plot(m3,type = "diagnostics"),
                    c("patchwork","ggplot2::ggplot","ggplot",
                      "ggplot2::gg","S7_object","gg"))
    expect_s3_class(plot(m3,type = "acme"),
                    c("patchwork","ggplot2::ggplot","ggplot",
                      "ggplot2::gg","S7_object","gg"))
    expect_s3_class(plot(m3,type = "ade"),
                    c("patchwork","ggplot2::ggplot","ggplot",
                      "ggplot2::gg","S7_object","gg"))
    expect_s3_class(plot(m3,type = c("acme","ade")),
                    c("patchwork","ggplot2::ggplot","ggplot",
                      "ggplot2::gg","S7_object","gg"))
  }
  
})
