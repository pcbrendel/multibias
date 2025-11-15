# Simulated data with exposure misclassification and selection bias

Data containing two sources of bias, three known confounders, and
100,000 observations. This data is obtained by sampling with replacement
with probability = *S* from `df_em_sel_source` then removing the columns
*X* and *S*. The resulting data corresponds to what a researcher would
see in the real-world: a misclassified exposure, *Xstar*, and missing
data for those not selected into the study (*S*=0). As seen in
`df_em_sel_source`, the true, unbiased exposure-outcome odds ratio = 2.

## Usage

``` r
df_em_sel
```

## Format

A dataframe with 100,000 rows and 5 columns:

- Xstar:

  misclassified exposure, 1 = present and 0 = absent

- Y:

  outcome, 1 = present and 0 = absent

- C1:

  1st confounder, 1 = present and 0 = absent

- C2:

  2nd confounder, 1 = present and 0 = absent

- C3:

  3rd confounder, 1 = present and 0 = absent
