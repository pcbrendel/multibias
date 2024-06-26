# multibias 1.5.1

* The following functions now support more flexible combinations of continuous
  and binary exposure-outcome variables:
  * `adjust_uc`
  * `adjust_emc` (exposure must be binary)
  * `adjust_omc` (outcome must be binary)
  * `adjust_sel`
  * `adjust_uc_sel`
* Expanded the number of known confounders in dataframes:
  * `df_uc_omc`
  * `df_uc_omc_source`
  * `df_uc_emc`
  * `df_uc_emc_source`
* Dataframes `df_uc` and `df_uc_source` now both have continuous and
  binary exposures and outcomes.

# multibias 1.5.0

## New features

* Added two functions for simultaneous adjustment of uncontrolled confounding, 
  outcome misclassification, and selection bias: `adjust_uc_omc_sel` & 
  `adjust_multinom_uc_omc_sel`.
* Added dataframes with uncontrolled confounding, outcome misclassification, 
  and selection bias: `df_uc_omc_sel` and `df_uc_omc_sel_source`.
* Expanded the number of known confounders in dataframes:
  * `df_uc_sel`
  * `df_uc_sel_source`

# multibias 1.4.0

## New features

* Added two functions for simultaneous adjustment of exposure misclassification
  and outcome misclassification: `adjust_emc_omc` & `adjust_multinom_emc_omc`.
* Added dataframes with exposure misclassification and outcome
  misclassification: `df_emc_omc` and `df_emc_omc_source`.
* Expanded the number of known confounders in dataframes:
  * `df_emc_sel`
  * `df_emc_sel_source`

## Bug fixes

* Improved some of the documentation of equations.

# multibias 1.3.0

## New features

* Added a function for simultaneous adjustment of outcome misclassification
  and selection bias: `adjust_omc_sel`.
* Added dataframes with outcome misclassification and selection bias:
  `df_omc_sel` and `df_omc_sel_source`.
* Expanded the number of known confounders in dataframes:
  * `df_uc`
  * `df_uc_source`
  * `df_emc`
  * `df_emc_source`
  * `df_omc`
  * `df_omc_source`
  * `df_sel`
  * `df_sel_source`

## Bug fixes

* Fixed bug in `adjust_omc` that appears when using three confounders

# multibias 1.2.1

* Moved examples from README to vignette.

# multibias 1.2.0

## New features

* Added two functions for simultaneous adjustment of uncontrolled confounding
  and outcome misclassification: `adjust_uc_omc` and `adjust_multinom_uc_omc`.
* Added dataframes with uncontrolled confounding and outcome misclassification:
  `df_uc_omc` and `df_uc_omc_source`.

## Bug fixes

* None

# multibias 1.1.0

## New features

* Created new function to adjust for outcome misclassification: `adjust_omc`.
* Added dataframes for all single bias scenarios:
  * `df_emc`
  * `df_emc_source`
  * `df_omc`
  * `df_omc_source`
  * `df_sel`
  * `df_sel_source`
  * `df_uc`
  * `df_uc_source`

## Bug fixes

* `adjust_sel` had been weighing with the probability of selection
  instead of the **inverse** probability of selection.

# multibias 1.0.0

* Initial CRAN submission.