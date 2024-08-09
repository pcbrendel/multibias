# multibias <img src="man/figures/logo.png" align="right" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/pcbrendel/multibias/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pcbrendel/multibias/actions/workflows/R-CMD-check.yaml)
[![cranlogs](https://cranlogs.r-pkg.org/badges/multibias)](https://cran.r-project.org/package=multibias)
<!-- badges: end -->

## Overview

The multibias package is used to adjust for multiple biases in causal inference when working with observational data. Bias here refers to the case when the associational estimate of effect (e.g., $P(Y=1|X=1,C=0) / P(Y=1|X=0,C=0)$) does not equal the causal estimate of effect (e.g., $P(Y^{X=1}=1) / P(Y^{X=0}=1)$). The underlying methods are explained in the [article](https://doi.org/10.1093/ije/dyad001):

Brendel PB, Torres AZ, Arah OA, Simultaneous adjustment of uncontrolled confounding, selection bias and misclassification in multiple-bias modelling, *International Journal of Epidemiology*, Volume 52, Issue 4, Pages 1220–1230

The functions provide odds ratio estimates adjusted for any combination of: uncontrolled confounding (**uc**), exposure misclassification (**emc**), outcome misclassification (**omc**), and selection bias (**sel**):

| Function | Adjusts for |
| -------- | ----------- |
| `adjust_emc()` | exposure misclassification |
| `adjust_omc()` | outcome misclassification |
| `adjust_sel()` | selection bias |
| `adjust_uc()` | uncontrolled confounding |
| `adjust_emc_sel()` | exposure misclassification & selection bias |
| `adjust_emc_omc` | exposure misclassification & outcome misclassification |
| `adjust_omc_sel()` | outcome misclassification & selection bias |
| `adjust_uc_emc()` | uncontrolled confounding & exposure misclassificaiton |
| `adjust_uc_omc()` | uncontrolled confounding & outcome misclassification |
| `adjust_uc_sel()` | uncontrolled confounding & selection bias |
| `adjust_uc_emc_sel()` & `adjust_multinom_uc_emc_sel()` | uncontrolled confounding, exposure misclassification, & selection bias |
| `adjust_uc_omc_sel()` & `adjust_multinom_uc_omc_sel()` | uncontrolled confounding, outcome misclassification, & selection bias |

To use these functions without R programming, go to the [multibias Shiny app](https://pcbrendel.shinyapps.io/multibias/).

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

## Coming soon
* Final triple bias combinations
* Support for continuous exposure and outcome
