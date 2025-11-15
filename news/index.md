# Changelog

## multibias 1.7.2

CRAN release: 2025-06-15

- Added
  [`multibias_plot()`](http://www.paulbrendel.com/multibias/reference/multibias_plot.md)
  to visualize sensitivity analysis results
- When using validation data in
  [`multibias_adjust()`](http://www.paulbrendel.com/multibias/reference/multibias_adjust.md)
  the function now incorporates uncertainty of the effect estimates from
  the validation data by sampling from each estimate’s mean and SE. Now,
  when using validation data, the confidence intervals from multibias
  bootstrapped results will represent two sources of uncertainty: random
  error and systematic error.
- Added FAQ documentation

## multibias 1.7.1

CRAN release: 2025-05-10

- Updated code with dynamic formula construction so that there is no
  limit to the number of known confounders one can include when using
  `bias_params` as an input for
  [`multibias_adjust()`](http://www.paulbrendel.com/multibias/reference/multibias_adjust.md)
- [`multibias_adjust()`](http://www.paulbrendel.com/multibias/reference/multibias_adjust.md)
  now has built in bootstrapping
- Added [`summary()`](https://rdrr.io/r/base/summary.html) method to
  `data_observed`

## multibias 1.7

CRAN release: 2025-04-08

- Created `bias_params` class to handle bias parameter inputs to
  [`multibias_adjust()`](http://www.paulbrendel.com/multibias/reference/multibias_adjust.md)
- Replaced the various `adjust()` functions with a single
  [`multibias_adjust()`](http://www.paulbrendel.com/multibias/reference/multibias_adjust.md)
  function. Users now specify the biases they want to adjust for in the
  `data_observed` object. Bias adjustment formulas are now found in the
  `bias_params` documentation.
- The user now specifies biases for adjustment in the `bias` input of
  `data_observed`
- Removed `evans` data; now only used in vignette

## multibias 1.6.3

CRAN release: 2025-02-23

- Created a `pkgdown` web page: www.paulbrendel.com/multibias
- Refined the vignette, including a new NHANES analysis

## multibias 1.6.2

CRAN release: 2025-01-08

- The following functions now accept `data_validation` as an input for
  bias adjustment:
  - `adjust_om_sel.R`
  - `adjust_uc_sel.R`
  - `adjust_uc_em.R`
  - `adjust_uc_om.R`
  - `adjust_uc_em_sel.R`
  - `adjust_uc_om_sel.R`

## multibias 1.6.1

CRAN release: 2024-12-02

- The following functions now accept `data_validation` as an input for
  bias adjustment:
  - `adjust_em_om.R`
  - `adjust_em_sel.R`
- Bug fixes for validation data input in `adjust_em.R` and `adjust_om.R`
- Bug fixes for data and printing in `data_observed` and
  `data_validation`

## multibias 1.6

CRAN release: 2024-10-26

- Created new class `data_observed` to represent observed causal data
- All `adjust` functions now take `data_observed` as input
- Created new class `data_validation` to represent causal data that can
  be used as validaiton data for bias adjustment
- The following functions now accept `data_validation` as an input for
  bias adjustment:
  - `adjust_uc.R`
  - `adjust_em.R`
  - `adjust_om.R`
  - `adjust_sel.R`

## multibias 1.5.3

CRAN release: 2024-09-22

- All exposure misclassificaiton naming changed from `emc` changed to
  `em`
- All outcome misclassificaiton naming changed from `omc` changed to
  `om`
- Added lifecycle badges for above function renames
- Merged `adjust_multinom_uc_em_sel` into `adjust_uc_em_sel`
- Merged `adjust_multinom_uc_om_sel` into `adjust_uc_om_sel`
- The following functions now support more flexible combinations of
  continuous and binary exposure-outcome variables:
  - `adjust_uc_em_sel.R`
  - `adjust_uc_om_sel.R`

## multibias 1.5.2

CRAN release: 2024-08-21

- Merged `adjust_multinom_emc_omc` into `adjust_emc_omc`
- Merged `adjust_multinom_uc_emc` into `adjust_uc_emc`
- Merged `adjust_multinom_uc_omc` into `adjust_uc_omc`
- The following functions now support more flexible combinations of
  continuous and binary exposure-outcome variables:
  - `adjust_emc_sel` (exposure must be binary)
  - `adjust_omc_sel` (outcome must be binary)
  - `adjust_uc_emc` (exposure must be binary)
  - `adjust_uc_omc` (outcome must be binary)
  - `adjust_multinom_uc_emc` (exposure must be binary)
  - `adjust_multinom_uc_omc` (outcome must be binary)
- Expanded the number of known confounders in dataframes:
  - `df_omc_sel`
  - `df_omc_sel_source`

## multibias 1.5.1

CRAN release: 2024-06-20

- The following functions now support more flexible combinations of
  continuous and binary exposure-outcome variables:
  - `adjust_uc`
  - `adjust_emc` (exposure must be binary)
  - `adjust_omc` (outcome must be binary)
  - `adjust_sel`
  - `adjust_uc_sel`
- Expanded the number of known confounders in dataframes:
  - `df_uc_omc`
  - `df_uc_omc_source`
  - `df_uc_emc`
  - `df_uc_emc_source`
- Dataframes `df_uc` and `df_uc_source` now both have continuous and
  binary exposures and outcomes.

## multibias 1.5.0

CRAN release: 2024-05-05

### New features

- Added two functions for simultaneous adjustment of uncontrolled
  confounding, outcome misclassification, and selection bias:
  `adjust_uc_omc_sel` & `adjust_multinom_uc_omc_sel`.
- Added dataframes with uncontrolled confounding, outcome
  misclassification, and selection bias: `df_uc_omc_sel` and
  `df_uc_omc_sel_source`.
- Expanded the number of known confounders in dataframes:
  - `df_uc_sel`
  - `df_uc_sel_source`

## multibias 1.4.0

CRAN release: 2024-01-27

### New features

- Added two functions for simultaneous adjustment of exposure
  misclassification and outcome misclassification: `adjust_emc_omc` &
  `adjust_multinom_emc_omc`.
- Added dataframes with exposure misclassification and outcome
  misclassification: `df_emc_omc` and `df_emc_omc_source`.
- Expanded the number of known confounders in dataframes:
  - `df_emc_sel`
  - `df_emc_sel_source`

### Bug fixes

- Improved some of the documentation of equations.

## multibias 1.3.0

CRAN release: 2023-12-11

### New features

- Added a function for simultaneous adjustment of outcome
  misclassification and selection bias: `adjust_omc_sel`.
- Added dataframes with outcome misclassification and selection bias:
  `df_omc_sel` and `df_omc_sel_source`.
- Expanded the number of known confounders in dataframes:
  - `df_uc`
  - `df_uc_source`
  - `df_emc`
  - `df_emc_source`
  - `df_omc`
  - `df_omc_source`
  - `df_sel`
  - `df_sel_source`

### Bug fixes

- Fixed bug in `adjust_omc` that appears when using three confounders

## multibias 1.2.1

CRAN release: 2023-10-21

- Moved examples from README to vignette.

## multibias 1.2.0

CRAN release: 2023-10-16

### New features

- Added two functions for simultaneous adjustment of uncontrolled
  confounding and outcome misclassification: `adjust_uc_omc` and
  `adjust_multinom_uc_omc`.
- Added dataframes with uncontrolled confounding and outcome
  misclassification: `df_uc_omc` and `df_uc_omc_source`.

### Bug fixes

- None

## multibias 1.1.0

CRAN release: 2023-10-06

### New features

- Created new function to adjust for outcome misclassification:
  `adjust_omc`.
- Added dataframes for all single bias scenarios:
  - `df_emc`
  - `df_emc_source`
  - `df_omc`
  - `df_omc_source`
  - `df_sel`
  - `df_sel_source`
  - `df_uc`
  - `df_uc_source`

### Bug fixes

- `adjust_sel` had been weighing with the probability of selection
  instead of the **inverse** probability of selection.

## multibias 1.0.0

CRAN release: 2023-09-21

- Initial CRAN submission.
