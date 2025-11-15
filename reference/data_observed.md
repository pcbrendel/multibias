# Represent observed causal data

`data_observed` combines the observed dataframe with specific
identification of the columns corresponding to the exposure, outcome,
and confounders. It is an essential input of the
[`multibias_adjust()`](http://www.paulbrendel.com/multibias/reference/multibias_adjust.md)
function.

## Usage

``` r
data_observed(data, bias, exposure, outcome, confounders = NULL)
```

## Arguments

- data:

  Dataframe for bias analysis.

- bias:

  String type(s) of bias distorting the effect of the exposure on the
  outcome. Can choose from a subset of the following: "uc", "em", "om",
  "sel". These correspond to uncontrolled confounding, exposure
  misclassification, outcome misclassification, and selection bias,
  respectively.

- exposure:

  String name of the column in `data` corresponding to the exposure
  variable.

- outcome:

  String name of the column in `data` corresponding to the outcome
  variable.

- confounders:

  String name(s) of the column(s) in `data` corresponding to the
  confounding variable(s).

## Value

An object of class `data_observed` containing:

- data:

  A dataframe with the selected columns

- bias:

  The type(s) of bias present

- exposure:

  The name of the exposure variable

- outcome:

  The name of the outcome variable

- confounders:

  The name(s) of the confounder variable(s)

## Examples

``` r
df <- data_observed(
  data = df_sel,
  bias = "uc",
  exposure = "X",
  outcome = "Y",
  confounders = c("C1", "C2", "C3")
)
```
