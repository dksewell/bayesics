
go_fast_for_cran_checks = TRUE

test_that("Test Wilcoxon signed rank analysis",{
  
  # Test small sample
  if(!go_fast_for_cran_checks){
    N = 15
    test_data_small = 
      data.frame(x = rbeta(N,2,10),
                 y = rbeta(N,5,10))
    
    ## Test input
    expect_no_error(
      fita <- 
        wilcoxon_test_b(test_data_small$x - test_data_small$y)
    )
    expect_no_error(
      fitb <- 
        wilcoxon_test_b(test_data_small$x,
                        test_data_small$y,
                        paired = TRUE)
    )
    
    ## Test output
    expect_equal(fita[c(1:5)],
                 fitb[c(1:5)],
                 tolerance = 0.05)
    expect_s3_class(fita$posterior_distribution,
                    c("tbl_df", "tbl", "data.frame"))
    expect_type(fita$posterior_mean,"double")
    expect_type(fita$CI,"double")
    expect_type(fita$Pr_less_than_p,"double")
    expect_type(fita$Pr_in_ROPE,"double")
    expect_type(fita$BF_for_phi_gr_onehalf_vs_phi_less_onehalf,"double")
    
    ## Test priors
    expect_no_error(
      wilcoxon_test_b(test_data_small$x - test_data_small$y,
                      prior = "uniform")
    )
    expect_no_error(
      wilcoxon_test_b(test_data_small$x - test_data_small$y,
                      prior_shapes = c(5,5))
    )
    
    ## Test ROPE
    expect_no_error(
      wilcoxon_test_b(test_data_small$x - test_data_small$y,
                      ROPE = 0.1)
    )
    expect_no_error(
      wilcoxon_test_b(test_data_small$x - test_data_small$y,
                      ROPE = c(0.4,0.65))
    )
  }
  
  
  # Large samples
  N = 150
  set.seed(2025)
  test_data_big = 
    data.frame(x = rbeta(N,2,10),
               y = rbeta(N,5,10))
  
  ## Test input
  expect_no_error(
    fitc <- 
      wilcoxon_test_b(test_data_big$x - test_data_big$y)
  )
  expect_no_error(
    fitd <- 
      wilcoxon_test_b(test_data_big$x,
                      test_data_big$y,
                      paired = TRUE)
  )
  
  ## Test output
  expect_equal(fitc[c(1:5)],
               fitd[c(1:5)],
               tolerance = 0.05)
  expect_type(fitc$posterior_parameters,"double")
  expect_type(fitc$posterior_mean,"double")
  expect_type(fitc$CI,"double")
  expect_type(fitc$Pr_less_than_p,"double")
  expect_type(fitc$Pr_in_ROPE,"double")
  expect_type(fitc$BF_for_phi_gr_onehalf_vs_phi_less_onehalf,"double")
  
  ## Test priors
  expect_no_error(
    wilcoxon_test_b(test_data_big$x - test_data_big$y,
                    prior = "uniform")
  )
  expect_no_error(
    wilcoxon_test_b(test_data_big$x - test_data_big$y,
                    prior_shapes = c(5,5))
  )
  
  ## Test ROPE
  expect_no_error(
    wilcoxon_test_b(test_data_big$x - test_data_big$y,
                    ROPE = 0.1)
  )
  expect_no_error(
    wilcoxon_test_b(test_data_big$x - test_data_big$y,
                    ROPE = c(0.4,0.65))
  )
  
  
})


test_that("Test Wilcoxon rank sum analysis",{
  
  # Small samples
  if(!go_fast_for_cran_checks){
    set.seed(2025)
    N = 15
    x = rbeta(N,2,10)
    y = rbeta(N + 1,5,10)
    
    ## Test input
    expect_no_error(
      fita <-
        wilcoxon_test_b(x,y)
    )
    expect_s3_class(fita$posterior_distribution,
                    c("tbl_df", "tbl", "data.frame"))
    expect_type(fita$posterior_mean,"double")
    expect_type(fita$CI,"double")
    expect_type(fita$Pr_less_than_p,"double")
    expect_type(fita$Pr_in_ROPE,"double")
    expect_type(fita$BF_for_Omegax_gr_onehalf_vs_Omegax_less_onehalf,"double")
    expect_s3_class(fita$prob_plot,
                    c("patchwork","ggplot2::ggplot","ggplot",
                      "ggplot2::gg","S7_object","gg"))
    
    ## Test priors
    expect_no_error(
      wilcoxon_test_b(x,
                      y,
                      prior = "uniform")
    )
    expect_no_error(
      wilcoxon_test_b(x,
                      y,
                      prior_shapes = c(5,5))
    )
    
    ## Test ROPE
    expect_no_error(
      wilcoxon_test_b(x,
                      y,
                      ROPE = 0.1)
    )
    expect_no_error(
      wilcoxon_test_b(x,
                      y,
                      ROPE = c(0.1,0.8))
    )
  }
  
  
  
  # Large samples
  set.seed(2025)
  N = 150
  x = rbeta(N,2,10)
  y = rbeta(N + 1,5,10)
  
  
  ## Test input
  expect_no_error(
    fitb <-
      wilcoxon_test_b(x,y)
  )
  expect_type(fitb$posterior_parameters,"double")
  expect_type(fitb$posterior_mean,"double")
  expect_type(fitb$CI,"double")
  expect_type(fitb$Pr_less_than_p,"double")
  expect_type(fitb$Pr_in_ROPE,"double")
  expect_type(fitb$BF_for_Omegax_gr_onehalf_vs_Omegax_less_onehalf,"double")
  expect_s3_class(fitb$prob_plot,
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  ## Test priors
  expect_no_error(
    wilcoxon_test_b(x,
                    y,
                    prior = "uniform")
  )
  expect_no_error(
    wilcoxon_test_b(x,
                    y,
                    prior_shapes = c(5,5))
  )
  
  ## Test ROPE
  expect_no_error(
    wilcoxon_test_b(x,
                    y,
                    ROPE = 0.1)
  )
  expect_no_error(
    wilcoxon_test_b(x,
                    y,
                    ROPE = c(0.1,0.8))
  )
  
})
