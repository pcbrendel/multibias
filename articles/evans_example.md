# Multibias example: Evans

This vignette demonstrates how to use `multibias` to adjust for
uncontrolled confounding in a real-world dataset from the Evans County
Heart Study. It specifically showcases how to reason and derive bias
parameters that can be used to adjust for the uncontrolled confounding.
We’ll examine the relationship between smoking and coronary heart
disease (CHD), showing how failing to account for age as a confounder
can bias our estimates and how `multibias` can be used to arrive at the
unbiased effect estimate despite missing data on the confounder.

``` r
library(multibias)
library(dplyr)
```

The Evans County Heart Study was a prospective cohort study conducted in
Evans County, Georgia, from 1960 to 1969. The study aimed to investigate
risk factors for cardiovascular disease in a rural population. For this
example, we’ll use a subset of the data focusing on 609 participants
aged 40 and older.

The key variables in our analysis are:

- `SMK`: Smoking status (1 = smoker, 0 = non-smoker)
- `CHD`: Coronary heart disease status (1 = present, 0 = absent)
- `HPT`: Hypertension status (1 = present, 0 = absent)
- `AGE`: Age in years

Let’s load and examine the data:

``` r
evans <- read.csv("evans.csv")

summary_stats <- evans %>%
  summarise(
    n = n(),
    mean_age = mean(AGE),
    sd_age = sd(AGE),
    prop_smokers = mean(SMK),
    prop_chd = mean(CHD),
    prop_hpt = mean(HPT)
  )

print(summary_stats)
#>     n mean_age   sd_age prop_smokers  prop_chd  prop_hpt
#> 1 609 53.70608 9.258388     0.635468 0.1165846 0.4187192
```

For purposes of demonstrating `multibias`, let’s pretend that our data
was missing information on the confounding variable, AGE. Let’s create
our `data_observed` object with uncontrolled confounding (bias = “uc”)
and inspect the biased SMK-CHD effect estimate, adjusted for
hypertension (HPT).

``` r
df_obs <- data_observed(
  data = evans,
  bias = "uc",
  exposure = "SMK",
  outcome = "CHD",
  confounders = "HPT"
)

print(df_obs)
#> Observed Data
#> ---------------------------------
#> The following biases are present: 
#> Uncontrolled Confounding 
#> ---------------------------------
#> Exposure: SMK 
#> Outcome: CHD 
#> Confounders: HPT 
#> ---------------------------------
#> Data head: 
#>   SMK CHD HPT
#> 1   0   0   0
#> 2   1   0   0
#> 3   1   1   1
#> 4   1   0   1
#> 5   1   0   0
summary(df_obs)
#> Note: Estimates are exponentiated (odds ratios) for binary outcomes
#> # A tibble: 3 × 7
#>   term        estimate std.error statistic  p.value conf.low conf.high
#>   <chr>          <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
#> 1 (Intercept)   0.0529     0.299     -9.83 7.97e-23   0.0284    0.0920
#> 2 SMK           1.99       0.294      2.34 1.92e- 2   1.14      3.64  
#> 3 HPT           2.39       0.260      3.36 7.86e- 4   1.45      4.02
```

Can we anticipate whether this odds ratio without age-adjustment is
biased towards or away from the null? Let’s consider the association of
the uncontrolled confounder with the exposure and outcome.

``` r
cor(evans$SMK, evans$AGE)
#> [1] -0.1391298
cor(evans$CHD, evans$AGE)
#> [1] 0.1393077
```

In our data, AGE has a negative association with SMK (older people are
**less** likely to be smokers) and a positive association with CHD
(older people are **more** likely to have CHD). These opposite
associations must be biasing the odds ratio towards the null, creating a
distortion where those who are less likely to smoke are more likely to
experience the outcome.

We’ll treat AGE as a binary indicator of over (1) or under (0) age 60.
To adjust for the uncontrolled confounding from AGE, let’s refer to the
appropriate bias model for a binary uncontrolled confounder:
logit(P(U=1)) = α₀ + α₁X + α₂Y + α_(2+j)C_(j).

To derive the necessary bias parameters, let’s make the following
assumption:

- the odds of an age\>60 is half as likely in smokers than non-smokers
- the odds of an age\>60 is 2.5x in those with CHD than those without
  CHD
- the odds of an age\>60 is 2x in those with HPT than those without HPT

To convert these relationships as parameters in the model, we’ll
log-transform them from odds ratios. For the model intercept, we can use
the following reasoning: what is the probability that a non-smoker (X=0)
without CHD (Y=0) and HPT (C=0) is over age 60 in this population? We’ll
assume this is a 25% probability. We’ll use the inverse logit function
[`qlogis()`](https://rdrr.io/r/stats/Logistic.html) from the `stats`
package to convert this from a probability to the intercept coefficient
of the logistic regression model.

``` r
u_0 <- qlogis(0.25)
u_x <- log(0.5)
u_y <- log(2.5)
u_c <- log(2)

u_coefs <- list(u = c(u_0, u_x, u_y, u_c))
```

Now let’s plug these bias parameters into
[`multibias_adjust()`](http://www.paulbrendel.com/multibias/reference/multibias_adjust.md)
along with our `data_observed` object to obtain a bias-adjusted effect
estimate.

``` r
set.seed(1234)
multibias_adjust(
  data_observed = df_obs,
  bias_params = bias_params(coef_list = u_coefs),
  bootstrap = TRUE,
  bootstrap_reps = 100
)
#> $estimate
#> [1] 2.229563
#> 
#> $std.error
#> [1] 0.6849506
#> 
#> $ci
#> [1] 1.364743 3.962324
```

We get an odds ratio of 2.2. This matches our expectation that the bias
adjustment would pull the odds ratio away from the null. How does this
result compare to the result we would get if age **wasn’t** missing in
the data and was incorporated in the outcome regression?

``` r
full_model <- glm(CHD ~ SMK + HPT + AGE,
  family = binomial(link = "logit"),
  data = evans
)
or <- round(exp(coef(full_model)[2]), 2)
or_ci_low <- round(
  exp(coef(full_model)[2] - 1.96 * summary(full_model)$coef[2, 2]), 2
)
or_ci_high <- round(
  exp(coef(full_model)[2] + 1.96 * summary(full_model)$coef[2, 2]), 2
)

print(paste0("Odds Ratio: ", or))
#> [1] "Odds Ratio: 2.31"
print(paste0("95% CI: (", or_ci_low, ", ", or_ci_high, ")"))
#> [1] "95% CI: (1.28, 4.16)"
```

Based on these results, it appears that the bias-adjusted odds ratio
obtained via `multibias` is close to this complete-data odds ratio of
2.3.
