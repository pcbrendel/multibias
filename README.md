# multibias <img src="man/figures/logo.png" align="right" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/pcbrendel/multibias/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pcbrendel/multibias/actions/workflows/R-CMD-check.yaml)
[![cranlogs](https://cranlogs.r-pkg.org/badges/multibias)](https://cran.r-project.org/package=multibias)
<!-- badges: end -->

## Overview

The multibias package is used to adjust for multiple biases in causal inference when working with observational data. Bias here refers to the case when the associational estimate of effect (e.g., $`P(Y=1|X=1,C=0) / P(Y=1|X=0,C=0)`$) does not equal the causal estimate of effect (e.g., $`P(Y^{X=1}=1) / P(Y^{X=0}=1)`$). The underlying methods are explained in the [article](https://doi.org/10.1093/ije/dyad001):

Brendel PB, Torres AZ, Arah OA, Simultaneous adjustment of uncontrolled confounding, selection bias and misclassification in multiple-bias modelling, *International Journal of Epidemiology*, Volume 52, Issue 4, Pages 1220–1230

The functions provide odds ratio estimates adjusted for any combination of: uncontrolled confounding (**uc**), exposure misclassification (**em**), outcome misclassification (**om**), and selection bias (**sel**).

Single bias adjustments:

| Function | Adjusts for |
| -------- | ----------- |
| `adjust_em()` | exposure misclassification |
| `adjust_om()` | outcome misclassification |
| `adjust_sel()` | selection bias |
| `adjust_uc()` | uncontrolled confounding |

Multiple bias adjustments:

| Function | Adjusts for |
| -------- | ----------- |
| `adjust_em_sel()` | exposure misclassification & selection bias |
| `adjust_em_om` | exposure misclassification & outcome misclassification |
| `adjust_om_sel()` | outcome misclassification & selection bias |
| `adjust_uc_em()` | uncontrolled confounding & exposure misclassificaiton |
| `adjust_uc_om()` | uncontrolled confounding & outcome misclassification |
| `adjust_uc_sel()` | uncontrolled confounding & selection bias |
| `adjust_uc_em_sel()` | uncontrolled confounding, exposure misclassification, & selection bias |
| `adjust_uc_om_sel()` | uncontrolled confounding, outcome misclassification, & selection bias |

The package also includes several dataframes that are useful for demonstrating and validating the bias adjustment methods. Each dataframe contains different combinations of bias as identified by the same prefixing system (e.g., **uc** for uncontrolled confounding). For each bias combination, there is a dataframe with incomplete information (as would be encountered in the real world) (e.g., `df_uc`) and a dataframe with complete information that was used to derive the biased data (e.g., `df_uc_source`).

If you are new to bias analysis, check out [Applying Quantitative Bias Analysis to Epidemiologic Data](https://link.springer.com/book/10.1007/978-0-387-87959-8) or visit my [website](https://www.paulbrendel.com). For examples, see the vignette.

## Installation

```{r, eval = FALSE}
# install from CRAN
install.packages("multibias")

# install from github using devtools
# library("devtools")
devtools::install_github("pcbrendel/multibias")
```

## How it works

1. Determine the desired biases to adjust for in your observational data for a given exposure-outcome effect and identify the corresponding `adjust` function.
2. Obtain the necessary bias parameters (see function documentation) for bias adjustment. These values could come from the literature, validation data, or expert opinion. Each parameter can be represented as a single value or as a probability distribution.
3. Run the `adjust` function after inputting:
   1. The observational data
   2. Column names of the exposure, outcome, and measured confounders
   3. The bias parameters
   4. Level for outupt confidence interval
4. The `adjust` function will output the bias-adjusted exposure-outcome odds ratio and confidence interval.