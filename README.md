# bayesics 2.0.0

Bayesian statistics is an integral part of contemporary applied science.  **bayesics** provides a single framework, unified in syntax and output, for performing the most commonly used statistical procedures, ranging from one- and two-sample inference to general mediation analysis.  **bayesics** leans hard away from the requirement that users be familiar with algorithms by using closed-form solutions whenever possible, and automatically selecting the number of posterior samples required for accurate inference when such solutions are not possible. **bayesics** focuses instead on providing key inferential quantities: point estimates, credible intervals, probability of direction, region of practical equivalance (ROPE), and, when applicable, Bayes factors.  While *algorithmic* assessment is not required in **bayesics**, *model* assessment is still critical; towards that, **bayesics** provides diagnostic plots for parametric inference, including Bayesian p-values. Finally, **bayesics** provides extensions to models implemented within **bayesics** or in alternative \R packages and, in the case of mediation analysis, correction to existing implementations.


# Installation
``` r
# Development version from GitHub
# install.packages("devtools")
devtools::install_github("a1arakkal/bayesics")
```

# Basic usage
``` r
# Load in an example dataset
data(indo_rct,
     package = "medicaldata")

# Run two-sample difference in mean analysis
t_test_b(age ~ rx,
         data = indo_rct)

# Run two-sample difference in proportion analysis
## Create the contingency table
gender_table =
  table(indo_rct$gender,
        indo_rct$rx)
## Perform the analysis
prop_test_b(gender_table[1,],
            gender_table[2,],
            CI_level = 0.99)

# Run the non-parametric Bayesian Wilcoxon rank sum analysis
wilcoxon_test_b(
  indo_rct$risk[which(indo_rct$rx == "1_indomethacin")],
  indo_rct$risk[which(indo_rct$rx == "0_placebo")]
)

# Run a Bernoulli GLM
## Fit the model
pep_fit =
  glm_b(outcome ~ age + gender + risk + rx,
        data = indo_rct,
        family = binomial())

## Plot results, including diagnostic plots
plot(pep_fit)

## Look at summary
summary(pep_fit)

## We could also have used a non-parametric loss-likelihood bootstrap
pep_np_fit =
  np_glm_b(outcome ~ age + gender + risk + rx,
           data = indo_rct,
           family = binomial())
```
