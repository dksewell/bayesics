
test_that("Simple case-control analysis", {
  
  # Test input
  expect_no_error(
    case_control_b(matrix(c(8,47,1,26),2,2))
  )
  
  expect_no_error(
    case_control_b(c(8,47),
                   c(1,26))
  )
  
  expect_no_error(
    case_control_b(x = matrix(c(8,47,1,26),2,2))
  )
  
  # Test large sample
  expect_no_error(
    case_control_b(x = 5 + matrix(c(8,47,1,26),2,2))
  )
  
  # Test ROPE
  expect_no_error(
    case_control_b(x = matrix(c(8,47,1,26),2,2),
                   ROPE = 1.05)
  )
  
  expect_no_error(
    case_control_b(x = 5 + matrix(c(8,47,1,26),2,2),
                   ROPE = 1.05)
  )
  
  # Test prior
  expect_no_error(
    case_control_b(x = 5 + matrix(c(8,47,1,26),2,2),
                   ROPE = 1.05)
  )
  
  expect_no_error(
    case_control_b(x = 5 + matrix(c(8,47,1,26),2,2),
                   ROPE = 1.05,
                   prior_mean = 10)
  )
  
  expect_no_error(
    case_control_b(x = 5 + matrix(c(8,47,1,26),2,2),
                   ROPE = 1.05,
                   prior_sd = 0.01)
  )
  
})