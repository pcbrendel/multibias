# Represent validation causal data

`data_validation` is one of two different options to represent bias
assumptions for bias adjustment. It combines the validation dataframe
with specific identification of the appropriate columns for bias
adjustment, including: true exposure, true outcome, confounders,
misclassified exposure, misclassified outcome, and selection. The
purpose of validation data is to use an external data source to
transport the necessary causal relationships that are missing in the
observed data.

## Usage

``` r
data_validation(
  data,
  true_exposure,
  true_outcome,
  confounders = NULL,
  misclassified_exposure = NULL,
  misclassified_outcome = NULL,
  selection = NULL
)
```

## Arguments

- data:

  Dataframe of validation data

- true_exposure:

  String name of the column in `data` corresponding to the true
  exposure.

- true_outcome:

  String name of the column in `data` corresponding to the true outcome.

- confounders:

  String name(s) of the column(s) in `data` corresponding to the
  confounding variable(s).

- misclassified_exposure:

  String name of the column in `data` corresponding to the misclassified
  exposure.

- misclassified_outcome:

  String name of the column in `data` corresponding to the misclassified
  outcome.

- selection:

  String name of the column in `data` corresponding to the selection
  indicator.

## Value

An object of class `data_validation` containing:

- data:

  A dataframe with the selected columns

- true_exposure:

  The name of the true exposure variable

- true_outcome:

  The name of the true outcome variable

- confounders:

  The name(s) of the confounder variable(s)

- misclassified_exposure:

  The name of the misclassified exposure variable

- misclassified_outcome:

  The name of the misclassified outcome variable

- selection:

  The name of the selection indicator variable

## Examples

``` r
df <- data_validation(
  data = df_sel_source,
  true_exposure = "X",
  true_outcome = "Y",
  confounders = c("C1", "C2", "C3"),
  selection = "S"
)
```
