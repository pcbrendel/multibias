
<!-- README.md is generated from README.Rmd. Please edit that file -->

# multibias <img src="man/figures/logo.png" align="right" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/pcbrendel/multibias/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pcbrendel/multibias/actions/workflows/R-CMD-check.yaml)
[![cranlogs](https://cranlogs.r-pkg.org/badges/multibias)](https://cran.r-project.org/package=multibias)
<!-- badges: end -->

## Overview

The multibias package is used to adjust for multiple biases in causal
inference when working with observational data. Bias here refers to the
case when the associational estimate of effect does not equal the causal
estimate of effect:

$$(P(Y=1|X=1,C=0) / P(Y=1|X=0,C=0)) \neq (P(Y^{X=1}=1) / P(Y^{X=0}=1))$$

The `multibias_adjust()` function outputs odds ratio estimates adjusted
for any combination of: uncontrolled confounding (**uc**), exposure
misclassification (**em**), outcome misclassification (**om**), and
selection bias (**sel**).

The package also includes several dataframes that are useful for
validating the bias adjustment methods. Each dataframe contains
different combinations of bias as identified by the same prefixing
system. For each bias combination, there is a dataframe with incomplete
information (as would be encountered in the real world) (e.g., `df_uc`)
and a dataframe with complete information that was used to derive the
biased data (e.g., `df_uc_source`).

## Installation

``` r
# install from CRAN
install.packages("multibias")

# install from github using devtools
# library("devtools")
devtools::install_github("pcbrendel/multibias")
```

## Getting started

1.  Represent the observed causal data as a `data_observed` object. Here
    you provide the data, specify the key variables, and list the biases
    present in the data. See list below for the different bias
    combinations that multibias can handle.
2.  Obtain one of the two sources for bias adjustment:
    1.  Bias parameters - via the `bias_params` object. Values for these
        parameters could come from the literature, validation data, or
        expert opinion. Each parameter can be represented as a single
        value or as a probability distribution. See the `bias_params`
        documentation for the full bias models.
    2.  Validation dataframe - via the `data_validation` object. The
        purpose of validation data is to use an external data source to
        transport the necessary causal relationships that are missing in
        the observed data.
3.  Run `multibias_adjust()` using the above inputs to obtain the
    bias-adjusted exposure-outcome odds ratio and confidence interval.
4.  Visualize a Forest Plot of the observed effect estimate against
    various bias-adjusted estimates via `multibias_plot()`.

### Possible bias adjustments

**Single Bias**

- exposure misclassification
- outcome misclassification
- selection bias
- uncontrolled confounding

**Multiple Biases**

- exposure misclassification & selection bias
- exposure misclassification & outcome misclassification
- outcome misclassification & selection bias
- uncontrolled confounding & exposure misclassificaiton
- uncontrolled confounding & outcome misclassification
- uncontrolled confounding & selection bias
- uncontrolled confounding, exposure misclassification, & selection bias
- uncontrolled confounding, outcome misclassification, & selection bias

## Resources

- Brendel PB, Torres AZ, Arah OA, Simultaneous adjustment of
  uncontrolled confounding, selection bias and misclassification in
  multiple-bias modelling, *International Journal of Epidemiology*,
  Volume 52, Issue 4, Pages 1220–1230.
  <https://doi.org/10.1093/ije/dyad001>
- [Applying Quantitative Bias Analysis to Epidemiologic
  Data](https://link.springer.com/book/10.1007/978-0-387-87959-8)
