test_confounders <- function(data_observed) {
  df <- data_observed$data
  exposure <- data_observed$exposure
  outcome <- data_observed$outcome
  confounders <- data_observed$confounders

  results <- list()

  for (i in seq_along(confounders)) {
    var <- confounders[i]
    var_data <- df[[var]]
    confounder_exposure_association <- NA
    confounder_outcome_association <- NA

    if (is_binary(var_data)) {
      c_x_ttest <- t.test(
        var_data[df[[exposure]] == 1],
        var_data[df[[exposure]] == 0]
      )
      c_x_cor <- list(
        test = "Student's t-Test (two-sided)",
        statistic = unname(c_x_ttest$statistic),
        p.val = c_x_ttest$p.value
      )
    } else if (is_continuous(var_data)) {
      c_x_cor <- cor.test(
        var_data, df[[exposure]],
        method = "pearson"
      )
      c_x_cor <- list(
        test = "Pearson Correlation Coefficient (two-sided)",
        statistic = unname(c_x_cor$statistic),
        p.val = c_x_cor$p.value
      )
    }

    results[[var]] <- list(
      variable = var,
      confounder_exposure_association = c_x_cor,
      confounder_outcome_association = NA
    )
  }
  return(results)
}
