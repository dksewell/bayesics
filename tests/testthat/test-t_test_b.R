
test_that("Test t_test_b",{
  
  # Test inputs
  expect_no_error(
    fita <- 
      t_test_b(rnorm(50))
  )
  expect_no_error(fita)
  expect_s3_class(plot(fita),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_no_error(
    fitb <-
      t_test_b(outcome ~ 1,
               data = data.frame(outcome = rnorm(50)))
  )
  expect_no_error(fitb)
  expect_s3_class(plot(fitb),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_no_error(
    fitc <-
      t_test_b(rnorm(50),
               rnorm(15,1))
  )
  expect_no_error(fitc)
  expect_s3_class(plot(fitc),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  
  expect_no_error(
    fitd <-
      t_test_b(outcome ~ asdf,
               data = 
                 data.frame(outcome = c(rnorm(50),
                                        rnorm(15,1)),
                            asdf = rep(c("a","b"),c(50,15))))
  )
  expect_no_error(fitd)
  expect_s3_class(plot(fitd),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_no_error(
    fite <-
      t_test_b(rnorm(50),
               rnorm(50,1),
               paired = TRUE)
  )
  expect_no_error(fite)
  expect_s3_class(plot(fite),
                  c("patchwork","ggplot2::ggplot","ggplot",
                    "ggplot2::gg","S7_object","gg"))
  expect_error(
    t_test_b(rnorm(50),
             rnorm(15,1), # Different length should throw an error if paired = TRUE
             paired = TRUE)
  )
  
  
  
})