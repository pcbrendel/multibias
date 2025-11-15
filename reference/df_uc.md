# Simulated data with uncontrolled confounding

Data containing one source of bias, three known confounders, and 100,000
observations. This data is obtained from `df_uc_source` by removing the
column *U*. The resulting data corresponds to what a researcher would
see in the real-world: information on known confounders (*C1*, *C2*, and
*C3*), but not for confounder *U*. As seen in `df_uc_source`, the true,
unbiased exposure-outcome effect estimate = 2.

## Usage

``` r
df_uc
```

## Format

A dataframe with 100,000 rows and 7 columns:

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
