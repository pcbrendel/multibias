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