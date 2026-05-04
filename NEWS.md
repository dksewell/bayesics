# bayesics 3.0.0

* Major restructuring of the S3 structure in bayesics.  All regression objects now inherit the lm_b class, and methods now rely heavily on those, better enabling expandability of the bayesics functionality.  Objects have been standardized so in most cases generics will work similarly on any object with the same estimation algorithm. Added plot_dx, plot_bands functions, as well as get_posterior_draws function for aov_b, glm_b, and np_glm_b objects.


# bayesics 2.1.1

* Fixed bug relating to response variable transformations and improper prior
* Greatly improved conjugate prior defaults in lm_b for the intercept
* Fixed default `fractional_proportion` in `frac_bayes_factors()` to be symmetric and not dependent on model order.

# bayesics 2.1.0

* Added fractional Bayes factors for lm_b objects
* Added check in lm_b when finding default hyperparameters for sigma^2 in case optim fails

# bayesics 2.0.2

* Updated licence
* Updated documentation for aov_b to list residuals and standardized residuals in Values


# bayesics 2.0.1

Updated help for stats-like functions, and changed cat()/print() to message().


# bayesics 2.0.0

* Initial CRAN submission.
