
test_that("Test t_test_b",{
  
  # Test inputs
  expect_no_error(
    fita <- 
      t_test_b(rnorm(50))
  )
  expect_no_error(
    fitb <-
      t_test_b(outcome ~ 1,
               data = data.frame(outcome = rnorm(50)))
  )
  expect_no_error(
    fitc <-
      t_test_b(rnorm(50),
               rnorm(15,1))
  )
  expect_no_error(
    fitd <-
      t_test_b(outcome ~ asdf,
               data = 
                 data.frame(outcome = c(rnorm(50),
                                        rnorm(15,1)),
                            asdf = rep(c("a","b"),c(50,15))))
  )
  expect_no_error(
    fite <-
      t_test_b(rnorm(50),
               rnorm(50,1),
               paired = TRUE)
  )
  expect_error(
    t_test_b(rnorm(50),
             rnorm(15,1), # Different length should throw an error if paired = TRUE
             paired = TRUE)
  )
  
  
  # Test output
  ## Single population
  expect_s3_class(fita,c("tbl_df", "tbl", "data.frame"))
  expect_type(fita$Variable,"character")
  expect_type(fita$`Post Mean`,"double")
  expect_type(fita$Lower,"double")
  expect_type(fita$Upper,"double")
  expect_type(fita$`Prob Dir`,"double")
  
  expect_s3_class(fitb,c("tbl_df", "tbl", "data.frame"))
  expect_type(fitb$Variable,"character")
  expect_type(fitb$`Post Mean`,"double")
  expect_type(fitb$Lower,"double")
  expect_type(fitb$Upper,"double")
  expect_type(fitb$`Prob Dir`,"double")
  
  expect_s3_class(fite,c("tbl_df", "tbl", "data.frame"))
  expect_type(fite$Variable,"character")
  expect_type(fite$`Post Mean`,"double")
  expect_type(fite$Lower,"double")
  expect_type(fite$Upper,"double")
  expect_type(fite$`Prob Dir`,"double")
  
  ## Two populations
  expect_s3_class(fitc$summary,c("tbl_df", "tbl", "data.frame"))
  expect_identical(colnames(fitc$summary),
                   c("Variable","Post Mean","Lower","Upper","Prob Dir"))
  expect_type(fitc$summary$Variable,"character")
  expect_type(fitc$summary$`Post Mean`,"double")
  expect_type(fitc$summary$Lower,"double")
  expect_type(fitc$summary$Upper,"double")
  expect_type(fitc$summary$`Prob Dir`,"double")
  
  expect_s3_class(fitc$pairwise,c("tbl_df", "tbl", "data.frame"))
  expect_type(fitc$pairwise$Comparison,"character")
  expect_type(fitc$pairwise$`Post Mean`,"double")
  expect_type(fitc$pairwise$Lower,"double")
  expect_type(fitc$pairwise$Upper,"double")
  expect_type(fitc$pairwise$`Prob Dir`,"double")
  expect_type(unlist(fitc$pairwise[,6]),"double")
  expect_type(fitc$pairwise$EPR,"double")
  expect_type(fitc$pairwise$`EPR Lower`,"double")
  expect_type(fitc$pairwise$`EPR Upper`,"double")
  
  
  
})