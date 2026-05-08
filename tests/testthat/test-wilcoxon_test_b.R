
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
    expect_no_error(fita)
    expect_s3_class(plot(fita),
                    c("patchwork","ggplot2::ggplot","ggplot",
                      "ggplot2::gg","S7_object","gg"))
    expect_no_error(
      fitb <- 
        wilcoxon_test_b(test_data_small$x,
                        test_data_small$y,
                        paired = TRUE)
    )
    expect_no_error(fitb)
    expect_s3_class(plot(fitb),
                    c("patchwork","ggplot2::ggplot","ggplot",
                      "ggplot2::gg","S7_object","gg"))
    
    ## Test output
    expect_equal(fita$results,
                 fitb$results,
                 tolerance = 0.05)
    
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
  expect_no_error(fitc)
  expect_s3_class(plot(fitc),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_no_error(
    fitd <- 
      wilcoxon_test_b(test_data_big$x,
                      test_data_big$y,
                      paired = TRUE)
  )
  expect_equal(fitc$results,
               fitd$results,
               tolerance = 0.05)
  
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
    expect_no_error(fita)
    expect_s3_class(plot(fita),
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
  expect_no_error(fitb)
  expect_s3_class(plot(fitb),
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
