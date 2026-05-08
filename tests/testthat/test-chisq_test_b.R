
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
    fit <- independence_b(round(N * dep_probs))
  )
  expect_no_error(fit)
  expect_warning(plot(fit))
  
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
                   sampling_design = "fixed rows")
  )
  expect_no_error(
    independence_b(round(N * ind_probs),
                   sampling_design = "fixed rows")
  )
  ## Try other priors
  expect_no_error(
    independence_b(round(N * dep_probs),
                   prior = "uniform",
                   sampling_design = "fixed rows")
  )
  expect_no_error(
    independence_b(round(N * dep_probs),
                   prior_shapes = 2,
                   sampling_design = "fixed rows")
  )
  expect_no_error(
    independence_b(round(N * dep_probs),
                   prior_shapes = matrix(1:(nR*nC),nR,nC),
                   sampling_design = "fixed rows")
  )
  expect_no_error(
    independence_b(round(N * dep_probs),
                   prior_shapes = rep(2,nR*nC),
                   sampling_design = "fixed rows")
  )
  
  # Fixed columns sampling design
  expect_no_error(
    independence_b(round(N * dep_probs),
                   sampling_design = "fixed columns")
  )
  expect_no_error(
    independence_b(round(N * ind_probs),
                   sampling_design = "fixed c")
  )
  ## Try other priors
  expect_no_error(
    independence_b(round(N * dep_probs),
                   prior = "uniform",
                   sampling_design = "fixed c")
  )
  expect_no_error(
    independence_b(round(N * dep_probs),
                   prior_shapes = 2,
                   sampling_design = "fixed c")
  )
  expect_no_error(
    independence_b(round(N * dep_probs),
                   prior_shapes = matrix(1:(nR*nC),nR,nC),
                   sampling_design = "fixed c")
  )
  expect_no_error(
    independence_b(round(N * dep_probs),
                   prior_shapes = rep(2,nR*nC),
                   sampling_design = "fixed c")
  )
  
  # Test to make sure rows and columns on t(x) are equivalent
  test1 = 
    independence_b(round(N * dep_probs),
                   sampling_design = "fixed r")
  test1$results = 
    test1$results |>
    mutate(row = 
             as.integer(stringr::str_extract(Quantity, "(?<=Row )\\d+")),
           col = 
             as.integer(stringr::str_extract(Quantity, "(?<=Col )\\d+"))
    )
  test2_m = test1_m = 
    matrix(0.0,5,3)
  test1_m[cbind(test1$results$row,
                test1$results$col)] = 
    test1$results$`Post Mean`
  test2 = 
    independence_b(round(N * dep_probs) |> t(),
                   sampling_design = "fixed c")
  test2$results = 
    test2$results |>
    mutate(row = 
             as.integer(stringr::str_extract(Quantity, "(?<=Row )\\d+")),
           col = 
             as.integer(stringr::str_extract(Quantity, "(?<=Col )\\d+"))
    )
  test2_m[cbind(test2$results$col,
                test2$results$row)] = 
    test2$results$`Post Mean`
  
  expect_true(all(near(test1_m,test2_m)))
  
  
  
})