
test_that("Test independence analysis for 2-way tables",{
  
  # Generate data
  set.seed(2025)
  N = 500
  nR = 5
  nC = 3
  dep_probs = 
    extraDistr::rdirichlet(1,rep(2,nR*nC)) |> 
    matrix(nR,nC)
  ind_probs = 
    tcrossprod(rowSums(dep_probs),
               colSums(dep_probs))
  
  # Test with big N
  expect_no_error(
    independence_b(round(N * dep_probs))
  )
  expect_no_error(
    independence_b(round(N * ind_probs))
  )
  ## Try other priors
  expect_no_error(
    independence_b(round(N * dep_probs),
                   prior = "uniform")
  )
  expect_no_error(
    independence_b(round(N * dep_probs),
                   prior_shapes = 2)
  )
  expect_no_error(
    independence_b(round(N * dep_probs),
                   prior_shapes = matrix(1:(nR*nC),nR,nC))
  )
  expect_no_error(
    independence_b(round(N * dep_probs),
                   prior_shapes = rep(2,nR*nC))
  )
  
  # Fixed rows sampling design
  expect_no_error(
    independence_b(round(N * dep_probs),
                   sampling_design = "rows")
  )
  expect_no_error(
    independence_b(round(N * ind_probs),
                   sampling_design = "rows")
  )
  ## Try other priors
  expect_no_error(
    independence_b(round(N * dep_probs),
                   prior = "uniform",
                   sampling_design = "rows")
  )
  expect_no_error(
    independence_b(round(N * dep_probs),
                   prior_shapes = 2,
                   sampling_design = "rows")
  )
  expect_no_error(
    independence_b(round(N * dep_probs),
                   prior_shapes = matrix(1:(nR*nC),nR,nC),
                   sampling_design = "rows")
  )
  expect_no_error(
    independence_b(round(N * dep_probs),
                   prior_shapes = rep(2,nR*nC),
                   sampling_design = "rows")
  )
  
  # Fixed columns sampling design
  expect_no_error(
    independence_b(round(N * dep_probs),
                   sampling_design = "fixed columns")
  )
  expect_no_error(
    independence_b(round(N * ind_probs),
                   sampling_design = "cols")
  )
  ## Try other priors
  expect_no_error(
    independence_b(round(N * dep_probs),
                   prior = "uniform",
                   sampling_design = "cols")
  )
  expect_no_error(
    independence_b(round(N * dep_probs),
                   prior_shapes = 2,
                   sampling_design = "cols")
  )
  expect_no_error(
    independence_b(round(N * dep_probs),
                   prior_shapes = matrix(1:(nR*nC),nR,nC),
                   sampling_design = "cols")
  )
  expect_no_error(
    independence_b(round(N * dep_probs),
                   prior_shapes = rep(2,nR*nC),
                   sampling_design = "cols")
  )
  
  # Test to make sure rows and columns on t(x) are equivalent
  test1 = 
    independence_b(round(N * dep_probs),
                   sampling_design = "rows")
  test2 = 
    independence_b(round(N * dep_probs) |> t(),
                   sampling_design = "cols")
  expect_true(all(near(test1$posterior_shapes,
                       t(test2$posterior_shapes))))
  expect_true(all(near(test1$posterior_mean,
                       t(test2$posterior_mean))))
  expect_true(all(near(test1$lower_bound,
                       t(test2$lower_bound))))
  expect_true(all(near(test1$individual_ROPE,
                       t(test2$individual_ROPE))))
  expect_true(all(near(test1$prob_p_j_given_i_less_than_p_j,
                       t(test2$prob_p_j_given_i_less_than_p_j))))
  
  
  # expect_no_error(
  # )
  # expect_no_error(
  # )
  # expect_no_error(
  # )
  # expect_no_error(
  # )
  # expect_no_error(
  # )
  # expect_no_error(
  # )
  # expect_no_error(
  # )
  # expect_no_error(
  # )
  # 
  
  
  
})