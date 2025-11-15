# Data source for `df_uc`

Data with complete information on one source of bias, three known
confounders, and 100,000 observations. This data is used to derive
`df_uc` and can be used to obtain bias parameters for purposes of
validating the simultaneous multi-bias adjustment method with `df_uc`.
With this source data, the fitted regression g(P(Y=1)) = α₀ + α₁X +
α₂C1 + α₃C2 + α₄C3 + α₅U shows that the true, unbiased exposure-outcome
effect estimate = 2 when:

1.  g = logit, Y = *Y_bi*, and X = *X_bi* or

2.  g = identity, Y = *Y_cont*, X = *X_cont*.

## Usage

``` r
df_uc_source
```

## Format

A dataframe with 100,000 rows and 8 columns:

- X_bi:

  binary exposure, 1 = present and 0 = absent

- X_cont:

  continuous exposure

- Y_bi:

  binary outcome corresponding to exposure *X_bi*, 1 = present and 0 =
  absent

- Y_cont:

  continuous outcome corresponding to exposure *X_cont*

- C1:

  1st confounder, 1 = present and 0 = absent

- C2:

  2nd confounder, 1 = present and 0 = absent

- C3:

  3rd confounder, 1 = present and 0 = absent

- U:

  uncontrolled confounder, 1 = present and 0 = absent
