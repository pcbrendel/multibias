# Simulated data with exposure misclassification and outcome misclassification

Data containing two sources of bias, three known confounders, and
100,000 observations. This data is obtained from `df_emc_omc_source` by
removing the columns *X* and *Y*. The resulting data corresponds to what
a researcher would see in the real-world: a misclassified exposure,
*Xstar*, and a misclassified outcome, *Ystar*. As seen in
`df_em_om_source`, the true, unbiased exposure-outcome odds ratio = 2.

## Usage

``` r
df_em_om
```

## Format

A dataframe with 100,000 rows and 5 columns:

- Xstar:

  misclassified exposure, 1 = present and 0 = absent

- Ystar:

  misclassified outcome, 1 = present and 0 = absent

- C1:

  1st confounder, 1 = present and 0 = absent

- C2:

  2nd confounder, 1 = present and 0 = absent

- C3:

  3rd confounder, 1 = present and 0 = absent
