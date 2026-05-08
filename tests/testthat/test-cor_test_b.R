
test_that("Testing Kendall's tau analysis",{
  
  # Generate data
  set.seed(2025)
  N = 50
  x = rnorm(N)
  y = x + 4 * rnorm(N)
  
  # Check input
  expect_no_error(
    cor_test_b(x,y)
  )
  expect_no_error(
    cor_test_b(~ asdf + qwer,
               data = data.frame(asdf = x,
                                 qwer = y))
  )
    
  # Check tau
  test1 = 
    cor_test_b(x,y,
               tau = 0.04)
  test2 = 
    cor_test_b(x,y,
               tau = 0.1)
  expect_gt(
    test1$pdir$pdir,
    test2$pdir$pdir
  )
  
  # Check ROPE
  test3 = 
    cor_test_b(x,y,
               ROPE = 0.1)
  test4 = 
    cor_test_b(x,y,
               ROPE = 0.2)
  expect_lt(
    test3$results$Pr_in_ROPE,
    test4$results$Pr_in_ROPE
  )
  
  # Check prior
  expect_no_error(
    cor_test_b(x,y,
               prior = "uniform")
  )
  expect_no_error(
    cor_test_b(x,y,
               prior = "negati")
  )
  expect_no_error(
    cor_test_b(x,y,
               prior = "pos")
  )
  expect_no_error(
    cor_test_b(x,y,
               prior_shapes = c(10,10))
  )
  
  
})