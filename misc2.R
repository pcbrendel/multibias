# causal data class
create_causal_data <- function(
  data,
  exposure,
  outcome,
  confounders = NULL,
  misclassified_exposure = NULL,
  misclassified_outcome = NULL,
  selection = NULL,
  remove = NULL
) {

  # Validate input
  stopifnot(
    is.data.frame(data),
    is.character(exposure) & length(exposure) == 1,
    is.character(outcome) & length(outcome) == 1,
    is.character(confounders),
    (is.character(misclassified_exposure) & length(misclassified_exposure) == 1) | is.null(misclassified_exposure),
    (is.character(misclassified_outcome) & length(misclassified_outcome) == 1) | is.null(misclassified_outcome),
    (is.character(selection) & length(selection) == 1) | is.null(selection)
  )

  # Create the object
  obj <- list(
    data = data[, !colnames(data) %in% remove],
    exposure = exposure,
    outcome = outcome,
    confounders = confounders,
    misclassified_exposure = misclassified_exposure,
    misclassified_outcome = misclassified_outcome,
    selection = selection
  )

  # Assign the class
  class(obj) <- "causal_data"
  return(obj)
}

print.causal_data <- function(x, ...) {
  cat("Causal Data Object\n")
  cat("------------------\n")
  cat("Exposure:", x$exposure, "\n")
  cat("Outcome:", x$outcome, "\n")
  if (!is.null(x$confounders)) {
    cat("Confounders:", paste(x$confounders, collapse = ", "), "\n")
  }
  if (!is.null(x$misclassified_exposure)) {
    cat("Misclassified exposure:", x$misclassified_exposure, "\n")
  }
  if (!is.null(x$misclassified_outcome)) {
    cat("Misclassified outcome:", x$misclassified_outcome, "\n")
  }
  if (!is.null(x$selection)) {
    cat("Selection:", x$selection, "\n")
  }
  cat("Data head: \n")
  print(head(x$data))
}

# ----

adjust_uc <- function(
  observed_data,
  validation_data,
  level = 0.95
) {

  # check that col names align
  if (observed_data$exposure != validation_data$exposure) {
    stop("Exposure in observed data must match that in validation data.")
  }
  if (observed_data$outcome != validation_data$outcome) {
    stop("Outcome in observed data must match that in validation data.")
  }
  if (all(observed_data$confounders %in% validation_data$confounders) == FALSE) {
    stop("All confounders in observed data must be present in validation data.")
  }

  # check that df_validation has one more confounder than df_observed
  if (length(validation_data$confounders) -
        length(observed_data$confounders) != 1) {
    stop("This function adjusts for unobserved confounding from one confounder:\nthe validation data should have one more confounder than the observed data.")
  }

  n <- nrow(observed_data$data)

  df <- data.frame(
    X = observed_data$data[, observed_data$exposure],
    Y = observed_data$data[, observed_data$outcome]
  )
  df <- bind_cols(df, observed_data$data[, observed_data$confounders])

  if (sum(df$Y %in% c(0, 1)) == n) {
    y_binary <- TRUE
  } else {
    y_binary <- FALSE
  }

  df_val <- data.frame(
    X = validation_data$data[, validation_data$exposure],
    Y = validation_data$data[, validation_data$outcome]
  )

  uc <- setdiff(validation_data$confounders, observed_data$confounders)
  df_val$U <- validation_data$data[, uc]
  df_val <- bind_cols(df_val, validation_data$data[, observed_data$confounders])

  if (sum(df_val$U %in% c(0, 1)) != n) {
    stop("Uncontrolled confounder from the validation data must be a binary integer.")
  }

  u_mod <- glm(U ~ X + Y + .,
               family = binomial(link = "logit"),
               data = df_val)

  # sequence along coefs to predict U in observed data
  u_mod_coefs <- coef(u_mod)
  u_pred <- u_mod_coefs[1]

  for (i in 2:length(u_mod_coefs)) {
    var_name <- names(u_mod_coefs)[i]
    u_pred <- u_pred + df[[var_name]] * u_mod_coefs[i]
  }

  df$Upred <- rbinom(n, 1, plogis(u_pred))

  if (y_binary) {
    final <- glm(
      Y ~ X + Upred + .,
      family = binomial(link = "logit"),
      data = df
    )
  } else {
    final <- lm(
      Y ~ X + Upred + .,
      data = df
    )
  }

  est <- summary(final)$coef[2, 1]
  se <- summary(final)$coef[2, 2]
  alpha <- 1 - level

  if (y_binary) {
    estimate <- exp(est)
    ci <- c(exp(est + se * qnorm(alpha / 2)),
            exp(est + se * qnorm(1 - alpha / 2)))
  } else {
    estimate <- est
    ci <- c(est + se * qnorm(alpha / 2),
            est + se * qnorm(1 - alpha / 2))
  }

  return(list(estimate = estimate, ci = ci))

}

df_observed <- create_causal_data(
  df_uc,
  exposure = "X_bi",
  outcome = "Y_bi",
  confounders = c("C1", "C2", "C3"),
  remove = c("X_cont", "Y_cont")
)

print(df_observed)

df_validation <- create_causal_data(
  df_uc_source,
  exposure = "X_bi",
  outcome = "Y_bi",
  confounders = c("C1", "C2", "C3", "U"),
  remove = c("X_cont", "Y_cont")
)

print(df_validation)

adjust_uc(df_observed, df_validation)

# ----

adjust_emc <- function(
  observed_data,
  validation_data,
  level = 0.95
) {

  # check that col names align
  if (observed_data$misclassified_exposure != validation_data$misclassified_exposure) {
    stop("Misclassified exposure in observed data must match that in validation data.")
  }
  if (observed_data$outcome != validation_data$outcome) {
    stop("Outcome in observed data must match that in validation data.")
  }
  if (all(observed_data$confounders %in% validation_data$confounders) == FALSE) {
    stop("All confounders in observed data must be present in validation data.")
  }

  n <- nrow(observed_data$data)

  df <- data.frame(
    Xstar = observed_data$data[, observed_data$misclassified_exposure],
    Y = observed_data$data[, observed_data$outcome]
  )
  df <- bind_cols(df, observed_data$data[, observed_data$confounders])

  if (sum(df$Xstar %in% c(0, 1)) != n) {
    stop("Observed, misclassified exposure must be a binary integer.")
  }

  if (sum(df$Y %in% c(0, 1)) == n) {
    y_binary <- TRUE
  } else {
    y_binary <- FALSE
  }

  df_val <- data.frame(
    X = validation_data$data[, validation_data$exposure],
    Y = validation_data$data[, validation_data$outcome],
    Xstar = validation_data$data[, validation_data$exposure]
  )

  if (sum(df_val$X %in% c(0, 1)) != n) {
    stop("Validation, misclassified exposure must be a binary integer.")
  }
  if (sum(df_val$X %in% c(0, 1)) != n) {
    stop("Validation, true exposure must be a binary integer.")
  }

  x_mod <- glm(X ~ Xstar + Y + .,
               family = binomial(link = "logit"),
               data = df_val)



    df$Xpred <- rbinom(n, 1, plogis(x1_0 + x1_xstar * df$Xstar +
                                      x1_y * df$Y + x1_c1 * df$C1))

    if (y_binary) {
      final <- glm(
        Y ~ Xpred + C1,
        family = binomial(link = "logit"),
        data = df
      )
    } else {
      final <- lm(
        Y ~ Xpred + C1,
        data = df
      )
    }



  est <- summary(final)$coef[2, 1]
  se <- summary(final)$coef[2, 2]
  alpha <- 1 - level

  if (y_binary) {
    estimate <- exp(est)
    ci <- c(exp(est + se * qnorm(alpha / 2)),
            exp(est + se * qnorm(1 - alpha / 2)))
  } else {
    estimate <- est
    ci <- c(est + se * qnorm(alpha / 2),
            est + se * qnorm(1 - alpha / 2))
  }

  return(list(estimate = estimate, ci = ci))

}

df_observed <- create_causal_data(
  df_emc,
  exposure = "Xstar",
  outcome = "Y",
  confounders = c("C1", "C2", "C3"),
)

print(df_observed)

df_validation <- create_causal_data(
  df_emc_source,
  exposure = "X",
  outcome = "Y",
  confounders = c("C1", "C2", "C3"),
  misclassified_exposure = "Xstar"
)

print(df_validation)

adjust_uc(df_observed, df_validation)