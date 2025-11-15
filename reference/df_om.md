# Simulated data with outcome misclassification

Data containing one source of bias, three known confounders, and 100,000
observations. This data is obtained from `df_om_source` by removing the
column *Y*. The resulting data corresponds to what a researcher would
see in the real-world: a misclassified outcome, *Ystar*, and no data on
the true outcome. As seen in `df_om_source`, the true, unbiased
exposure-outcome odds ratio = 2.

## Usage

``` r
df_om
```

## Format

A dataframe with 100,000 rows and 5 columns:

- X:

  exposure, 1 = present and 0 = absent

- Ystar:

  misclassified outcome, 1 = present and 0 = absent

- C1:

  1st confounder, 1 = present and 0 = absent

- C2:

  2nd confounder, 1 = present and 0 = absent

- C3:

  3rd confounder, 1 = present and 0 = absent
