# Simulated data with exposure misclassification

Data containing one source of bias, three known confounders, and 100,000
observations. This data is obtained from `df_emc_source` by removing the
column *X*. The resulting data corresponds to what a researcher would
see in the real-world: a misclassified exposure, *Xstar*, and no data on
the true exposure. As seen in `df_emc_source`, the true, unbiased
exposure-outcome odds ratio = 2.

## Usage

``` r
df_em
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
