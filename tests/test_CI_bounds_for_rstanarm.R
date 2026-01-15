


run <- FALSE

if (run) {
  
  library(janitor)
  library(lubridate)
  library(tidyverse)
  library(magrittr)
  library(coda)
  library(rstanarm)
  library(future)
  
  # x = as.mcmc(cbind(rbeta(600,1,2),rnorm(600)))
  # 
  # test = raftery.diag(x)
  # 
  # qnorm(0.5 * (1.95) )^2 * 0.025 * 0.975 / 0.005^2
  
  
  # n = 100:1e4
  # k = round(0.025 * n)
  # test0 = (k/n - 0.025) * sqrt(n)
  # plot(test0,type='l')
  
  N = 100
  test_data = 
    data.frame(x1 = rnorm(N),
               x2 = rnorm(N),
               x3 = letters[1:5])
  test_data$outcome = 
    rbinom(N,1,1.0 / (1.0 + exp(-(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) ))))
  
  test_quantiles = function(i){
    set.seed(i)
    
    fit1 = 
      stan_glm(outcome ~ x1 + x2 + x3,
               data = test_data,
               family = "binomial")
    fit2 = 
      stan_glm(outcome ~ x1 + x2 + x3,
               data = test_data,
               family = "binomial",
               algorithm = "optimizing",
               importance_resampling = TRUE)
    
    bounds = c(posterior_interval(fit1,
                                  prob = 0.95)[expand.grid(2:3,1:2) |> 
                                                 as.matrix()],
               posterior_interval(fit2,
                                  prob = 0.95)[expand.grid(2:3,1:2) |> 
                                                 as.matrix()])
    return(bounds)
  }
  
  plan(multisession,
       workers = 20)
  
  quants_from_rstanarm = 
    future_sapply(1:1e3,test_quantiles)
  
  expand.grid(2:3,1:2) |> 
    as.matrix()                
  
  for(i in 1:nrow(quants_from_rstanarm)){
    cat("\n\n")
    cat(c("LB MCMC 1","LB MCMC 2",
          "UB MCMC 1","UB MCMC 2",
          "LB IS 1","LB IS 2",
          "UB IS 1","UB IS 2")[i])
    hist(quants_from_rstanarm[i,])
    abline(v = mean(quants_from_rstanarm[i,]))
    abline(v = c(0.99,1.01) * mean(quants_from_rstanarm[i,]))
    legend("topleft",
           legend = paste0("Mean outside of error bounds = ",
                           round(mean(
                             (quants_from_rstanarm[i,] < 
                                min(mean(quants_from_rstanarm[i,]) * c(0.99,1.01))) |
                               (quants_from_rstanarm[i,] >
                                  max(mean(quants_from_rstanarm[i,])*c(0.99,1.01)))
                           ), 3)),
           bty = 'n')
    if(i < nrow(quants_from_rstanarm)) readline() 
  }
  
  
}












