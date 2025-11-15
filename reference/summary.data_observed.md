# Summary method for data_observed objects

Provides a statistical summary of the observed data by fitting either:

- A logistic regression model for binary outcomes

- A linear regression model for continuous outcomes

The model includes the exposure and all confounders as predictors. For
binary outcomes, estimates are exponentiated to show odds ratios.

## Usage

``` r
# S3 method for class 'data_observed'
summary(object, ...)
```

## Arguments

- object:

  A `data_observed` object

- ...:

  Additional arguments passed to summary

## Value

A data frame containing model coefficients, standard errors, confidence
intervals, and p-values. For binary outcomes, coefficients are
exponentiated to show odds ratios.
