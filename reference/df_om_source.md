# Data source for `df_om`

Data with complete information on one sources of bias, three known
confounders, and 100,000 observations. This data is used to derive
`df_om` and can be used to obtain bias parameters for purposes of
validating the simultaneous multi-bias adjustment method with `df_om`.
With this source data, the fitted regression logit(P(Y=1)) = α₀ + α₁X +
α₂C1 + α₃C2 + α₄C3 shows that the true, unbiased exposure-outcome odds
ratio = 2.

## Usage

``` r
df_om_source
```

## Format

A dataframe with 100,000 rows and 6 columns:

- X:

  exposure, 1 = present and 0 = absent

- Y:

  true outcome, 1 = present and 0 = absent

- C1:

  1st confounder, 1 = present and 0 = absent

- C2:

  2nd confounder, 1 = present and 0 = absent

- C3:

  3rd confounder, 1 = present and 0 = absent

- Ystar:

  misclassified outcome, 1 = present and 0 = absent
