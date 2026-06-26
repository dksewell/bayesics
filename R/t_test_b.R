#' t-test
#' 
#' One and two sample t-tests on vectors of data
#' 
#' @details
#' A one and two sample t-test is nothing more than a special case of 
#' one-way anova.  See \code{\link{aov_b}} for details.
#' 
#' 
#' @param x Either a (non-empty) numeric vector of data values, or a formula 
#' of the form outcome ~ grouping variable.
#' @param y an optional (non-empty) numeric vector of data values
#' @param mu optional.  If supplied, \code{t_test_b} will return the 
#' posterior probabilty that the population mean (ignored in 2 sample inference) 
#' is less than this value.
#' @param paired logical.  If TRUE, provide both x and y as vectors.
#' @param data logical.  Only used if x is a formula.
#' @param heteroscedastic logical.  Set to FALSE to assume all groups have 
#' equal variance.
#' @param prior_mean_mu numeric. Hyperparameter for the a priori mean of the 
#' group means.
#' @param prior_mean_nu numeric. Hyperparameter which scales the precision of 
#' the group means.
#' @param prior_var_shape numeric. Twice the shape parameter for the inverse gamma prior on
#' the residual variance(s).  I.e., \eqn{\sigma^2\sim IG}(prior_var_shape/2,prior_var_rate/2).
#' @param prior_var_rate  numeric. Twice the rate parameter for the inverse gamma prior on
#' the residual variance(s).  I.e., \eqn{\sigma^2\sim IG}(prior_var_shape/2,prior_var_rate/2).
#' @param CI_level numeric. Credible interval level.
#' @param ROPE numeric.  Used to compute posterior probability that Cohen's D +/- ROPE
#' @param mc_error The number of posterior draws will ensure that with 99% 
#' probability the bounds of the credible intervals will be within \eqn{\pm} 
#' \code{mc_error}\eqn{\times 4s_y}, that is, within 100\code{mc_error}% of the 
#' trimmed range of y. (Ignored for single population inference.)
#' @param improper logical.  Should we use an improper prior that is proportional 
#' to the inverse of the variance?
#' @param seed integer.  Always set your seed!!!
#' @param plot logical. Should the resulting inverse gamma distribution be plotted?
#' 
#' @returns An object of class \code{\link{b_procedure-class}}.
#' 
#' 
#' @examples
#' \donttest{
#' # Single population
#' t_test_b(rnorm(50))
#' # or an alternative input format
#' t_test_b(outcome ~ 1,
#'          data = data.frame(outcome = rnorm(50)))
#' 
#' # Two populations
#' t_test_b(rnorm(50),
#'          rnorm(15,1))
#' 
#' # or an alternative input format
#' t_test_b(outcome ~ group_variable,
#'          data = 
#'            data.frame(outcome = c(rnorm(50),
#'                                   rnorm(15,1)),
#'                       group_variable = rep(c("a","b"),
#'                                            c(50,15))))
#' }
#' 
#' 
#' @export


t_test_b = function(x,
                    y,
                    mu,
                    paired = FALSE,
                    data,
                    heteroscedastic = TRUE,
                    prior_mean_mu,
                    prior_mean_nu = 0.001,
                    prior_var_shape = 0.001,
                    prior_var_rate = 0.001,
                    CI_level = 0.95,
                    ROPE = 0.1,
                    improper = FALSE,
                    plot = TRUE,
                    seed = 1,
                    mc_error = 0.002){
  
  outcome_name = NULL
  if(rlang::is_formula(x)){
    outcome_name = all.vars(x)[1]
    if(length(all.vars(x)) == 1){#Intercept only model, i.e., one population
      x = data[[outcome_name]]
    }
  }
  
  if(rlang::is_formula(x) & missing(data)) stop("If formula is given, data must also be given.")
  
  
  if(is.numeric(x)){
    
    # One sample inference
    if(missing(y) | paired){
      if(!missing(y) && paired && (length(x) != length(y)) ) stop("Length of x must equal that of y.")
      if(paired & !missing(y)){
        x = x - y
        outcome_name = "x minus y"
      }
      
      if(is.null(outcome_name)) outcome_name = "x"
      
      # Set alpha lv
      a = 1 - CI_level
      
      # Set prior_mean_mu if missing
      if(missing(prior_mean_mu)) prior_mean_mu = 0.0
      
      # Check if improper prior \propto 1/\sigma^2 is requested
      if(improper){
        prior_mean_mu = 0.0
        prior_mean_nu = 0.0
        prior_var_shape = -1.0
        prior_var_rate = 0.0
      }
      
      # Get summary stats
      data_quants = 
        tibble::tibble(n = NROW(x),
                       ybar = mean(x),
                       y2 = sum(x^2),
                       sample_var = var(x)) |> 
        mutate(s2 = (n - 1.0) / n * .data$sample_var)
      
      # Get posterior parameters
      nu_g = 
        prior_mean_nu + data_quants$n
      mu_g =
        (prior_mean_nu * prior_mean_mu + data_quants$n * data_quants$ybar) /
        nu_g
      a_G =
        prior_var_shape + sum(data_quants$n)
      b_G =
        prior_var_rate +
        sum(
          data_quants$n * data_quants$s2 +
            prior_mean_nu * data_quants$n / (nu_g + data_quants$n) * (prior_mean_mu - data_quants$ybar)^2
        )
      
      # Construct results
      results = 
        list(name = 
               "One sample population mean analysis",
             data = x,
             print_data = FALSE,
             CI_level = CI_level,
             prior = 
               paste0(
                 "Prior: mu ~ N(",
                 format(signif(prior_mean_mu, 3), 
                        scientific = FALSE),
                 ", sigma^2/",
                 format(signif(prior_mean_nu, 3), 
                        scientific = FALSE),
                 "), sigma^2 ~ IG(shape=",
                 format(signif(prior_var_shape, 3), 
                        scientific = FALSE),
                 "/2, rate=",
                 format(signif(prior_var_rate, 3), 
                        scientific = FALSE),
                 "/2)"
               )
        )
      
      # posterior summary
      results$results = 
        tibble(
          Quantity = 
            c("Population mean",
              "Population variance"),
          `Post Mean` = 
            c(mu_g, 
              b_G/2 / (a_G/2 - 1.0)),
          Lower = c(extraDistr::qlst(a/2, 
                                     df = a_G,
                                     mu = mu_g,
                                     sigma = sqrt(b_G / nu_g / a_G)),
                    extraDistr::qinvgamma(a/2, alpha = a_G/2, beta = b_G/2)),
          Upper = c(extraDistr::qlst(1 - a/2, 
                                     df = a_G,
                                     mu = mu_g,
                                     sigma = sqrt(b_G / nu_g / a_G)),
                    extraDistr::qinvgamma(1 - a/2, alpha = a_G/2, beta = b_G/2))
        )
      
      # Compute pdir
      if(paired){
        results$pdir = 
          list(pdir = extraDistr::plst(0, 
                                       df = a_G,
                                       mu = mu_g,
                                       sigma = sqrt(b_G / nu_g / a_G)))
        results$pdir$description = 
          paste0("Probability that the difference in means (x - y) is ",
                 ifelse(results$pdir$pdir > 0.5,
                        "less",
                        "greater"),
                 " than 0")
        results$pdir$pdir = 
          max(results$pdir$pdir,
              1.0 - results$pdir$pdir)
      }
      
      
      if(plot){
        
        results$plot =
          tibble::tibble(x = 
                           seq(
                             extraDistr::qlst(0.005,
                                              df = a_G,
                                              mu = mu_g,
                                              sigma = 
                                                sqrt( b_G / a_G / nu_g)
                             ),
                             extraDistr::qlst(0.995,
                                              df = a_G,
                                              mu = mu_g,
                                              sigma = 
                                                sqrt( b_G / a_G / nu_g)
                             ),
                             l = 50)) |> 
          ggplot(aes(x=x)) +
          stat_function(fun = 
                          function(x){
                            extraDistr::dlst(x,
                                             df = a_G,
                                             mu = mu_g,
                                             sigma = 
                                               sqrt( b_G / a_G / nu_g))
                          },
                        linewidth = 2) +
          theme_classic(base_size = 15) +
          xlab(expression(mu)) + 
          ylab("") + 
          ggtitle(ifelse(paired,
                         "Difference in means (x - y)",
                         "Population mean"))
      }
      
      results = 
        structure(results,
                  class = "b_procedure")
      
      return(results)
      
    }else{#End: one sample inference
      
      # Create data tibble
      ttest_data = 
        tibble::tibble(group = rep(c("x","y"),
                                   c(length(x),
                                     length(y)))) |> 
        dplyr::mutate(y = c(x,y))
      
      # Set prior_mean_mu if missing
      if(missing(prior_mean_mu))
        prior_mean_mu = mean(ttest_data$y)
      
      
      ret = 
        aov_b(y ~ group,
              data = ttest_data,
              heteroscedastic = heteroscedastic,
              prior_mean_mu = prior_mean_mu,
              prior_mean_nu = prior_mean_nu,
              prior_var_shape = prior_var_shape,
              prior_var_rate = prior_var_rate,
              CI_level = CI_level,
              ROPE = ROPE,
              improper = improper,
              seed = seed,
              mc_error = mc_error)
      
      # Create results object
      results = 
        list(
          name = "Two sample population means analysis",
          data = ttest_data,
          print_data = FALSE,
          CI_level = CI_level,
          prior =
            paste0(
              "Prior: mu ~ N(",
              format(signif(prior_mean_mu, 3), 
                     scientific = FALSE),
              ", sigma^2/",
              format(signif(prior_mean_nu, 3), 
                     scientific = FALSE),
              "), sigma^2 ~ IG(shape=",
              format(signif(prior_var_shape, 3), 
                     scientific = FALSE),
              "/2, rate=",
              format(signif(prior_var_rate, 3), 
                     scientific = FALSE),
              "/2)"
            ),
          notes = "ROPE for the difference in means is given in terms of Cohen's D."
        )
      
      # Get posterior summary
      ## Get each population's parameters
      if(heteroscedastic){
        results$results = 
          tibble(Quantity = 
                   c("Population 1 mean",
                     "Population 2 mean",
                     "Population 1 variance",
                     "Population 2 variance"))
      }else{
          results$results = 
            tibble(Quantity = 
                     c("Population 1 mean",
                       "Population 2 mean",
                       "Population 1 and 2 variance"))
      }
      
      for(j in c("Post Mean","Lower","Upper")) results$results[[j]] = ret$summary[[j]]
      
      ## Get difference in means
      results$results = 
        results$results |> 
        dplyr::bind_rows(
          tibble(
            Quantity = "Difference in population means (Pop 1 - Pop 2)",
            `Post Mean` = ret$pairwise_summary$`Post Mean`[1],
            Lower = ret$pairwise_summary$Lower[1],
            Upper = ret$pairwise_summary$Upper[1],
            Pr_in_ROPE = 
              ret$pairwise_summary |> 
              dplyr::pull(dplyr::contains("ROPE")),
            ROPE_lower_bound = -ROPE,
            ROPE_upper_bound = ROPE
          )
        )
      
      
      # Get pdir
      results$pdir = 
        list(
          pdir = ret$pairwise_summary$`Prob Dir`[1],
          description = 
            paste0("Probability that the difference in means (x - y) is ",
                   ifelse(ret$pairwise_summary$`Post Mean` < 0,
                          "less",
                          "greater"),
                   " than 0")
        )
      
      
      # Get BF
      results$BF = 
        list(
          description = "Bayes factor in favor of unequal group means",
          BF = ret$BF_for_different_vs_same_means
        )
      bf_max = max(results$BF$BF,
                   1.0 / results$BF$BF)
      results$BF$interpretation =
        ifelse(bf_max <= 3.2,
               "Not worth more than a bare mention",
               ifelse(bf_max <= 10,
                      "Substantial",
                      ifelse(bf_max <= 100,
                             "Strong",
                             "Decisive")))
      
      
      
      if(plot){
        results$plot = 
          tibble::tibble(x = 
                           seq(
                             min(
                               extraDistr::qlst(0.005,
                                                df = ret$posterior_parameters$a_g,
                                                mu = ret$posterior_parameters$mu_g,
                                                sigma = sqrt(ret$posterior_parameters$b_g / 
                                                               ret$posterior_parameters$nu_g / 
                                                               ret$posterior_parameters$a_g))
                             ),
                             max(
                               extraDistr::qlst(0.995,
                                                df = ret$posterior_parameters$a_g,
                                                mu = ret$posterior_parameters$mu_g,
                                                sigma = sqrt(ret$posterior_parameters$b_g / 
                                                               ret$posterior_parameters$nu_g / 
                                                               ret$posterior_parameters$a_g))
                             ),
                             l = 50)) |> 
          ggplot(aes(x=x)) +
          stat_function(fun = 
                          function(x){
                            extraDistr::dlst(x,
                                             df = ret$posterior_parameters$a_g[1],
                                             mu = ret$posterior_parameters$mu_g[1],
                                             sigma = sqrt(ret$posterior_parameters$b_g[1] / 
                                                            ret$posterior_parameters$nu_g[1] / 
                                                            ret$posterior_parameters$a_g[1]))
                          },
                        aes(color = "Posterior (Pop1)"),
                        linewidth = 2) +
          stat_function(fun = 
                          function(x){
                            extraDistr::dlst(x,
                                             df = ret$posterior_parameters$a_g[1 + heteroscedastic],
                                             mu = ret$posterior_parameters$mu_g[2],
                                             sigma = sqrt(ret$posterior_parameters$b_g[1 + heteroscedastic] / 
                                                            ret$posterior_parameters$nu_g[2] / 
                                                            ret$posterior_parameters$a_g[1 + heteroscedastic]))
                          },
                        aes(color = "Posterior (Pop2)"),
                        linewidth = 2)
        if(improper){
          post_modes = 
            extraDistr::dlst(0,
                             df = ret$posterior_parameters$a_g,
                             mu = 0,
                             sigma = sqrt(ret$posterior_parameters$b_g / 
                                            ret$posterior_parameters$nu_g / 
                                            ret$posterior_parameters$a_g)) |> 
            max()
          results$plot = 
            results$plot +
            stat_function(fun = 
                            function(x){
                              post_modes / 10
                            },
                          aes(color = "Prior"),
                          linewidth = 2) 
          # geom_hline(yintercept = post_modes / 10,
          #            aes(color = "Prior"),
          #            linewidth = 2)
        }else{
          results$plot =
            results$plot  +
            stat_function(fun = 
                            function(x){
                              extraDistr::dlst(x,
                                               df = ret$hyperparameters$a,
                                               mu = ret$hyperparameters$mu,
                                               sigma = 
                                                 ret$hyperparameters$b / 
                                                 ret$hyperparameters$a / 
                                                 ret$hyperparameters$nu)
                            },
                          aes(color = "Prior"),
                          linewidth = 2) 
        }
        results$plot = 
          results$plot + 
          scale_color_manual(values = c("Prior" = "#440154FF", 
                                        "Posterior (Pop1)" = "#21908CFF", 
                                        "Posterior (Pop2)" = "#FDE725FF")) +
          theme_classic(base_size = 15) +
          xlab("") + 
          ylab("") + 
          labs(color = "Distribution") + 
          ggtitle("Population means")
      }
      
      
      # attach aov_b object
      results$object_fit = ret
      
      results = 
        structure(results,
                  class = "b_procedure")
      
      return(results)
    }
  }else{
    
    # Set prior_mean_mu if missing
    if(missing(prior_mean_mu))
      prior_mean_mu = mean(data[[outcome_name]])
    
    # If formula (which implies it must be two sample inference):
    ret = 
      aov_b(x,
            data = data,
            heteroscedastic = heteroscedastic,
            prior_mean_mu = prior_mean_mu,
            prior_mean_nu = prior_mean_nu,
            prior_var_shape = prior_var_shape,
            prior_var_rate = prior_var_rate,
            CI_level = CI_level,
            ROPE = ROPE,
            improper = improper,
            seed = seed,
            mc_error = mc_error)
    
    # Get factor levels
    factor_levels = 
      sapply(ret$summary$Variable,
             function(z) trimws(strsplit(z,":")[[1]][3])) |> 
      unique() |> 
      na.omit()
    
    
    # Create results object
    results = 
      list(
        name = "Two sample population means analysis",
        data = data,
        print_data = FALSE,
        CI_level = CI_level,
        prior =
          paste0(
            "Prior: mu ~ N(",
            format(signif(prior_mean_mu, 3), 
                   scientific = FALSE),
            ", sigma^2/",
            format(signif(prior_mean_nu, 3), 
                   scientific = FALSE),
            "), sigma^2 ~ IG(shape=",
            format(signif(prior_var_shape, 3), 
                   scientific = FALSE),
            "/2, rate=",
            format(signif(prior_var_rate, 3), 
                   scientific = FALSE),
            "/2)"
          ),
        notes = "ROPE for the difference in means is given in terms of Cohen's D."
      )
    
    # Get posterior summary
    ## Get each population's parameters
    if(heteroscedastic){
      results$results = 
        tibble(Quantity = 
                 c(paste0("Population ",
                          factor_levels[1],
                          " mean"),
                   paste0("Population ",
                          factor_levels[2],
                          " mean"),
                   paste0("Population ",
                          factor_levels[1],
                          " variance"),
                   paste0("Population ",
                          factor_levels[2],
                          " variance")))
    }else{
      results$results = 
        tibble(Quantity = 
                 c(paste0("Population ",
                          factor_levels[1],
                          " mean"),
                   paste0("Population ",
                          factor_levels[2],
                          " mean"),
                   paste0("Population ",
                          factor_levels[1],
                          " and ",
                          factor_levels[2],
                          " variance")))
    }
    
    for(j in c("Post Mean","Lower","Upper")) results$results[[j]] = ret$summary[[j]]
    
    ## Get difference in means
    results$results = 
      results$results |> 
      dplyr::bind_rows(
        tibble(
          Quantity = 
            paste0("Difference in population means (Pop ",
                   factor_levels[1]," - Pop ",
                   factor_levels[2],")"),
          `Post Mean` = ret$pairwise_summary$`Post Mean`[1],
          Lower = ret$pairwise_summary$Lower[1],
          Upper = ret$pairwise_summary$Upper[1],
          Pr_in_ROPE = 
            ret$pairwise_summary |> 
            dplyr::pull(dplyr::contains("ROPE")),
          ROPE_lower_bound = -ROPE,
          ROPE_upper_bound = ROPE
        )
      )
    
    
    # Get pdir
    results$pdir = 
      list(
        pdir = ret$pairwise_summary$`Prob Dir`[1],
        description = 
          paste0("Probability that the difference in means (",
                 factor_levels[1],
                 " - ",
                 factor_levels[2],
                 ") is ",
                 ifelse(ret$pairwise_summary$`Post Mean` < 0,
                        "less",
                        "greater"),
                 " than 0")
      )
    
    
    if(plot){
      
      color_labels = 
        paste("Posterior (",
              factor_levels,
              ")",
              sep = "")
      color_values = 
        c("#440154FF",
          "#21908CFF",
          "#FDE725FF")
      names(color_values) = 
        c("Prior", color_labels)
      
        
      post_means = 
        ret$summary |> 
        dplyr::filter(grepl("Mean : ",ret$summary$Variable)) |> 
        dplyr::pull(.data$`Post Mean`)
      post_sds = 
        sqrt(ret$summary |> 
               dplyr::filter(!grepl("Mean :",ret$summary$Variable)) |> 
               dplyr::pull(.data$`Post Mean`))
      results$plot = 
        tibble::tibble(x = 
                         seq(
                           min(
                             qnorm(0.005,
                                   post_means,
                                   post_sds)
                           ),
                           max(
                             qnorm(0.995,
                                   ret$summary |> 
                                     dplyr::filter(grepl("Mean : ",
                                                         ret$summary$Variable)) |> 
                                     dplyr::pull(.data$`Post Mean`),
                                   sqrt(ret$summary |> 
                                          dplyr::filter(!grepl("Mean : ",
                                                               ret$summary$Variable)) |> 
                                          dplyr::pull(.data$`Post Mean`)))
                           ),
                           l = 50)) |> 
        ggplot(aes(x=x)) +
        stat_function(fun = 
                        function(x){
                          dnorm(x,
                                post_means[1],
                                post_sds[1])
                        },
                      aes(color = color_labels[1]),
                      linewidth = 2) +
        stat_function(fun = 
                        function(x){
                          dnorm(x,
                                post_means[2],
                                post_sds[1 + heteroscedastic])
                        },
                      aes(color = color_labels[2]),
                      linewidth = 2)
      if(improper){
        post_modes = 
          dnorm(post_means,
                post_means,
                post_sds) |> 
          max()
        results$plot = 
          results$plot +
          geom_hline(yintercept = post_modes / 10,
                     aes(color = "Prior"),
                     linewidth = 2)
      }else{
        results$plot =
          results$plot  +
          stat_function(fun = 
                          function(x){
                            dlst(x,
                                 df = ret$hyperparameters$a,
                                 mu = ret$hyperparameters$mu,
                                 sigma = 
                                   ret$hyperparameters$b / 
                                   ret$hyperparameters$a / 
                                   ret$hyperparameters$nu)
                          },
                        aes(color = "Prior"),
                        linewidth = 2) 
      }
      results$plot = 
        results$plot + 
        scale_color_manual(values = color_values) +
        theme_classic(base_size = 15) +
        xlab("") + 
        ylab("") + 
        labs(color = "Distribution") + 
        ggtitle("Population means")
    }
    
    # attach aov_b object
    results$object_fit = ret
    
    results = 
      structure(results,
                class = "b_procedure")
    
    return(results)
  }
  
}