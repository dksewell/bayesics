
run <- FALSE

if (run) {
  
# Look into reasonable/achievable mc_errors
pacman::p_load(bayesics,
               future)

# aov_b -------------------------------------------------------------------

# small data
set.seed(2025)
N = 100
test_data = 
  data.frame(x1 = rep(letters[1:5],N/5))
test_data$outcome = 
  rnorm(N,-1 + 2 * (test_data$x1 %in% c("d","e")) )

test = 
  aov_b(outcome ~ x1,
        test_data,
        prior_mean_mu = 2,
        prior_mean_nu = 0.5,
        prior_var_shape = 0.01,
        prior_var_rate = 0.01,
        mc_error = 0.002)
test
test$mc_error
nrow(test$posterior_draws)

# medium
set.seed(2025)
N = 500
test_data = 
  data.frame(x1 = rep(letters[1:5],N/5))
test_data$outcome = 
  rnorm(N,-1 + 2 * (test_data$x1 %in% c("d","e")) )

test = 
  aov_b(outcome ~ x1,
        test_data,
        prior_mean_mu = 2,
        prior_mean_nu = 0.5,
        prior_var_shape = 0.01,
        prior_var_rate = 0.01,
        mc_error = 0.002)
test
test$mc_error
nrow(test$posterior_draws)

# large
set.seed(2025)
N = 5000
test_data = 
  data.frame(x1 = rep(letters[1:5],N/5))
test_data$outcome = 
  rnorm(N,-1 + 2 * (test_data$x1 %in% c("d","e")) )

test = 
  aov_b(outcome ~ x1,
        test_data,
        prior_mean_mu = 2,
        prior_mean_nu = 0.5,
        prior_var_shape = 0.01,
        prior_var_rate = 0.01,
        mc_error = 0.002,
        compute_bayes_factor = FALSE)
test
test$mc_error
nrow(test$posterior_draws)


rm(list=ls())


# case_control_b ----------------------------------------------------------

# small
case_control_b(matrix(c(8,47,1,26),2,2),
               large_sample_approx = FALSE,
               mc_error = 0.005)
case_control_b(c(8,47),
               c(1,26),
               large_sample_approx = FALSE,
               mc_error = 0.005)
case_control_b(x = matrix(c(8,47,1,26),2,2),
               large_sample_approx = FALSE,
               mc_error = 0.005)


# chisq_test --------------------------------------------------------------

set.seed(2025)
nR = 5
nC = 3
dep_probs = 
  extraDistr::rdirichlet(1,rep(2,nR*nC)) |> 
  matrix(nR,nC)
ind_probs = 
  tcrossprod(rowSums(dep_probs),
             colSums(dep_probs))

# Multinomial sampling design
## Test with big N
N = 500
### Check if answers are reasonable wrt the truth
independence_b(round(N * dep_probs))
independence_b(round(N * ind_probs))


rm(list=ls())


# glm_b -------------------------------------------------------------------

# BINARY
## small
set.seed(2025)
N = 100
test_data = 
  data.frame(x1 = rnorm(N),
             x2 = rnorm(N),
             x3 = letters[1:5])
test_data$outcome = 
  rbinom(N,1,1.0 / (1.0 + exp(-(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) ))))

test = 
  glm_b(outcome ~ x1 + x2 + x3,
        data = test_data,
        family = binomial(),
        seed = 2025,
        algorithm = "IS",
        mc_error = 0.01)
test
test$mc_error
length(test$importance_sampling_weights)

## medium
set.seed(2025)
N = 500
test_data = 
  data.frame(x1 = rnorm(N),
             x2 = rnorm(N),
             x3 = letters[1:5])
test_data$outcome = 
  rbinom(N,1,1.0 / (1.0 + exp(-(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) ))))

test = 
  glm_b(outcome ~ x1 + x2 + x3,
        data = test_data,
        family = binomial(),
        seed = 2025,
        algorithm = "IS",
        mc_error = 0.01)
test
test$mc_error
length(test$importance_sampling_weights)


## large
set.seed(2025)
N = 5000
test_data = 
  data.frame(x1 = rnorm(N),
             x2 = rnorm(N),
             x3 = letters[1:5])
test_data$outcome = 
  rbinom(N,1,1.0 / (1.0 + exp(-(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e")) ))))

test = 
  glm_b(outcome ~ x1 + x2 + x3,
        data = test_data,
        family = binomial(),
        seed = 2025,
        algorithm = "IS",
        mc_error = 0.01)
test
test$mc_error
length(test$importance_sampling_weights)

# POISSON
## small
set.seed(2025)
N = 100
test_data = 
  data.frame(x1 = rnorm(N),
             x2 = rnorm(N),
             x3 = letters[1:5],
             time = rexp(N))
test_data$outcome = 
  rpois(N,exp(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e"))) * test_data$time)

test = 
  glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
        data = test_data,
        family = poisson(),
        seed = 2025,
        algorithm = "IS",
        mc_error = 0.01)
test
test$mc_error
length(test$importance_sampling_weights)

## medium
set.seed(2025)
N = 500
test_data = 
  data.frame(x1 = rnorm(N),
             x2 = rnorm(N),
             x3 = letters[1:5],
             time = rexp(N))
test_data$outcome = 
  rpois(N,exp(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e"))) * test_data$time)

test = 
  glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
        data = test_data,
        family = poisson(),
        seed = 2025,
        algorithm = "IS",
        mc_error = 0.01)
test
test$mc_error
length(test$importance_sampling_weights)

## large
set.seed(2025)
N = 5000
test_data = 
  data.frame(x1 = rnorm(N),
             x2 = rnorm(N),
             x3 = letters[1:5],
             time = rexp(N))
test_data$outcome = 
  rpois(N,exp(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e"))) * test_data$time)

test = 
  glm_b(outcome ~ x1 + x2 + x3 + offset(log(time)),
        data = test_data,
        family = poisson(),
        seed = 2025,
        algorithm = "IS",
        mc_error = 0.01)
test
test$mc_error
length(test$importance_sampling_weights)


# NBINOM
# small
set.seed(2025)
N = 100
test_data = 
  data.frame(x1 = rnorm(N),
             x2 = rnorm(N),
             x3 = letters[1:5],
             time = rexp(N))
test_data$outcome = 
  rnbinom(N,
          mu = exp(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e"))) * test_data$time,
          size = 0.7)

test = 
  glm_b(formula = outcome ~ x1 + x2 + x3 + offset(log(time)),
        data = test_data,
        family = negbinom(),
        seed = 2025,
        algorithm = "IS",
        mc_error = 0.01)
test
test$mc_error
length(test$importance_sampling_weights)

# medium
set.seed(2025)
N = 500
test_data = 
  data.frame(x1 = rnorm(N),
             x2 = rnorm(N),
             x3 = letters[1:5],
             time = rexp(N))
test_data$outcome = 
  rnbinom(N,
          mu = exp(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e"))) * test_data$time,
          size = 0.7)

test = 
  glm_b(formula = outcome ~ x1 + x2 + x3 + offset(log(time)),
        data = test_data,
        family = negbinom(),
        seed = 2025,
        algorithm = "IS",
        mc_error = 0.01)
test
test$mc_error
length(test$importance_sampling_weights)

# large
set.seed(2025)
N = 5000
test_data = 
  data.frame(x1 = rnorm(N),
             x2 = rnorm(N),
             x3 = letters[1:5],
             time = rexp(N))
test_data$outcome = 
  rnbinom(N,
          mu = exp(-2 + test_data$x1 + 2 * (test_data$x3 %in% c("d","e"))) * test_data$time,
          size = 0.7)

test = 
  glm_b(formula = outcome ~ x1 + x2 + x3 + offset(log(time)),
        data = test_data,
        family = negbinom(),
        seed = 2025,
        algorithm = "IS",
        mc_error = 0.01)
test
test$mc_error
length(test$importance_sampling_weights)


rm(list=ls())



# poisson_test_b ----------------------------------------------------------

poisson_test_b(x = c(12,20),
               offset = c(10,9),
               mc_error = 0.002)



# prop_test_b -------------------------------------------------------------

prop_test_b(c(14,22),
            c(19,45),
            mc_error = 0.002)


}



