# multibias

<!-- badges: start -->
[![R-CMD-check](https://github.com/pcbrendel/multibias/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pcbrendel/multibias/actions/workflows/R-CMD-check.yaml)
[![cranlogs](https://cranlogs.r-pkg.org/badges/multibias)](https://cran.r-project.org/package=multibias)
<!-- badges: end -->

## Overview

The multibias package is used to adjust for multiple biases in causal inference. The underlying methods are explained in the [article](https://doi.org/10.1093/ije/dyad001):

Brendel PB, Torres AZ, Arah OA, Simultaneous adjustment of uncontrolled confounding, selection bias and misclassification in multiple-bias modelling, *International Journal of Epidemiology*, Volume 52, Issue 4, Pages 1220–1230

The functions provide odds ratio estimates adjusted for any combination of: uncontrolled confounding (**uc**), exposure misclassification (**emc**), outcome misclassification (**omc**), and selection bias (**sel**):

Single biases:
  - `adjust_emc()` adjusts for exposure misclassification
  - `adjust_omc()` adjusts for outcome misclassification
  - `adjust_sel()` adjusts for selection bias
  - `adjust_uc()` adjusts for uncontrolled confounding

Double biases:
  - `adjust_emc_sel()` adjusts for exposure misclassification and selection bias
  - `adjust_uc_emc()` & `adjust_multinom_uc_emc()` adjusts for uncontrolled confounding and exposure misclassificaiton
  - `adjust_uc_omc()` & `adjust_multinom_uc_omc()` adjusts for uncontrolled confounding and outcome misclassification
  - `adjust_uc_sel()` adjusts for uncontrolled confounding and selection bias

Triple biases:
  - `adjust_uc_emc_sel()` & `adjust_multinom_uc_emc_sel()` adjusts for uncontrolled confounding, exposure misclassification, and selection bias

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
* Bias adjustments with outcome misclassification
* Support for continuous exposure and outcome
