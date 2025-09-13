# Test Confounders

# the following functions feed into test_confounders():
# 1. test_assoc_iso()
# 2. test_assoc_stepwise()
# 3. test_cie_iso()
# 4. test_cie_stepwise()

# 1. Associational criteria (isolation)
test_assoc_iso <- function(data_observed) {
  df <- data_observed$data
  exposure <- data_observed$exposure
  outcome <- data_observed$outcome
  confounders <- data_observed$confounders

  results <- list()

  for (i in seq_along(confounders)) {
    var <- confounders[i]
    df_c <- df[[var]]

    # 1. confounder-exposure
    confounder_exposure_association <- NA

    if (is_binary(df_c)) {
      if (is_continuous(df[[exposure]])) {
        # t-test
        c_x_ttest <- t.test(
          df[[exposure]][df_c == 1], #?
          df[[exposure]][df_c == 0] #?
        )
        c_x_cor <- list(
          test = "Student's t-Test (two-sided)",
          statistic = unname(c_x_ttest$statistic),
          p.val = c_x_ttest$p.value
        )
      } else if (is_binary(df[[exposure]])) {
        # chi-sq
        c_x_chisq
        c_x_cor
      }
    } else if (is_continuous(df_c)) {
      if (is_continuous(df[[exposure]])) {
        # corr coef
        c_x_cor <- cor.test(
          var_data, df[[exposure]],
          method = "pearson"
        )
        c_x_cor <- list(
          test = "Pearson Correlation Coefficient (two-sided)",
          statistic = unname(c_x_cor$statistic),
          p.val = c_x_cor$p.value
        )
      } else if (is_binary(df[[exposure]])) {
        # t-test
        c_x_ttest <- t.test(
          df_c[df[[exposure]] == 1],
          df_c[df[[exposure]] == 0]
        )
        c_x_cor <- list(
          test = "Student's t-Test (two-sided)",
          statistic = unname(c_x_ttest$statistic),
          p.val = c_x_ttest$p.value
        )
      }
    }

    # 2. confounder-outcome
    confounder_outcome_association <- NA

    formula_adj <- as.formula(paste(outcome, "~", exposure, "+", var))
    model_adj <- lm(formula_adj, data = data)

    for (i in seq_along(confounders)) {
      var <- confounders[i]
      df_c <- df[[var]]

    }

    # need this twice
    results[[var]] <- list(
      variable = var,
      confounder_exposure_association = c_x_cor,
      confounder_outcome_association = NA
    )
  }
  return(results)
}

# 2. Associational criteria (stepwise)
test_assoc_stepwise <- function(data_observed) {
  df <- data_observed$data
  exposure <- data_observed$exposure
  outcome <- data_observed$outcome
  confounders <- data_observed$confounders

  for (i in seq_along(confounders)) {
    var <- confounders[i]
    df_c <- df[[var]]

  }

}

# 3. Change-in-estimate criteria (isolation)
test_cie_iso <- function(data_observed) {
  df <- data_observed$data
  exposure <- data_observed$exposure
  outcome <- data_observed$outcome
  confounders <- data_observed$confounders

  formula_crude <- as.formula(paste(outcome, "~", exposure))
  model_crude <- lm(formula_crude, data = df)
  beta_x_crude <- coef(model_crude)[exposure]

  for (i in seq_along(confounders)) {
    var <- confounders[i]
    df_c <- df[[var]]
  }

}

# 4. Change-in-estimate critera (stepwise)
test_cie_stepwise <- function(data_observed) {
  df <- data_observed$data
  exposure <- data_observed$exposure
  outcome <- data_observed$outcome
  confounders <- data_observed$confounders

  formula_crude <- as.formula(paste(outcome, "~", exposure))
  model_crude <- lm(formula_crude, data = df)
  beta_x_crude <- coef(model_crude)[exposure]

  for (i in seq_along(confounders)) {
    var <- confounders[i]
    df_c <- df[[var]]
  }

}

# 5. Final return
#' Test for confounders in the observed data

#' @param data_observed Object of class `data_observed` corresponding to the
#' data to perform bias analysis on.
#'
#' @export

test_confounders <- function(data_observed, stepwise = FALSE) {

  if (stepwise) {
    x <- test_assoc_stepwise()
    y <- test_cie_stepwise()
  } else {
    x <- test_assoc_iso()
    y <- test_cie_iso()
  }

}