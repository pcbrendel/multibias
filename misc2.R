# causal data class
causal_data <- function(
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

  if (all(observed_data$confounders %in% validation_data$confounders) == FALSE) {
    stop("All confounders in observed data must be present in validation data.")
  }

  if (length(validation_data$confounders) -
        length(observed_data$confounders) != 1) {
    stop("This function is adjusting for unobserved confounding from one confounder.\nValidation data must have one more confounder than the observed data.")
  }

  n <- nrow(observed_data$data)

  df <- data.frame(
    X = observed_data$data[, observed_data$exposure],
    Y = observed_data$data[, observed_data$outcome]
  )
  df <- bind_cols(df, observed_data$data[, observed_data$confounders])

  if (all(df$Y %in% 0:1)) {
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

  if (all(df$X %in% 0:1)) {
    if (!all(df_val$X %in% 0:1)) {
      stop("Exposures from both datasets must all be binary or all be continuous.")
    }
  }

  if (all(df$Y %in% 0:1)) {
    if (!all(df_val$Y %in% 0:1)) {
      stop("Outcomes from both datasets must all be binary or all be continuous.")
    }
  }

  if (!all(df_val$U %in% 0:1)) {
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

df_observed <- causal_data(
  df_uc,
  exposure = "X_bi",
  outcome = "Y_bi",
  confounders = c("C1", "C2", "C3"),
  remove = c("X_cont", "Y_cont")
)

print(df_observed)

df_validation <- causal_data(
  df_uc_source,
  exposure = "X_bi",
  outcome = "Y_bi",
  confounders = c("C1", "C2", "C3", "U"),
  remove = c("X_cont", "Y_cont")
)

print(df_validation)

adjust_uc(df_observed, df_validation)

# ----

adjust_em <- function(
  observed_data,
  validation_data,
  level = 0.95
) {

  if (all(observed_data$confounders %in% validation_data$confounders) == FALSE) {
    stop("All confounders in observed data must be present in validation data.")
  }

  if (is.null(validation_data$misclassified_exposure)) {
    stop("This function is adjusting for a misclassified exposure.\nValidation data must have an exposure and misclassified exposure specified.")
  }

  n <- nrow(observed_data$data)

  df <- data.frame(
    Xstar = observed_data$data[, observed_data$exposure],
    Y = observed_data$data[, observed_data$outcome]
  )
  df <- bind_cols(df, observed_data$data[, observed_data$confounders])

  if (!all(df$Xstar %in% 0:1)) {
    stop("Observed, misclassified exposure must be a binary integer.")
  }

  if (all(df$Y %in% 0:1)) {
    y_binary <- TRUE
  } else {
    y_binary <- FALSE
  }

  df_val <- data.frame(
    X = validation_data$data[, validation_data$exposure],
    Y = validation_data$data[, validation_data$outcome],
    Xstar = validation_data$data[, validation_data$misclassified_exposure]
  )

  if (all(df$Xstar %in% 0:1)) {
    if (!all(df_val$Xstar %in% 0:1) || !all(df_val$X %in% 0:1)) {
      stop("Exposures from both datasets must all be binary or all be continuous.")
    }
  }

  if (all(df$Y %in% 0:1)) {
    if (!all(df_val$Y %in% 0:1)) {
      stop("Outcomes from both datasets must all be binary or all be continuous.")
    }
  }

  if (!all(df_val$Xstar %in% 0:1)) {
    stop("Validation, misclassified exposure must be a binary integer.")
  }
  if (!all(df_val$X %in% 0:1)) {
    stop("Validation, true exposure must be a binary integer.")
  }

  x_mod <- glm(X ~ Xstar + Y + .,
               family = binomial(link = "logit"),
               data = df_val)

  # sequence along coefs to predict X in observed data
  x_mod_coefs <- coef(x_mod)
  x_pred <- x_mod_coefs[1]

  for (i in 2:length(x_mod_coefs)) {
    var_name <- names(x_mod_coefs)[i]
    x_pred <- x_pred + df[[var_name]] * x_mod_coefs[i]
  }

  df$Xpred <- rbinom(n, 1, plogis(x_pred))

  if (y_binary) {
    final <- glm(
      Y ~ Xpred + .,
      family = binomial(link = "logit"),
      data = df
    )
  } else {
    final <- lm(
      Y ~ Xpred + .,
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

df_observed <- causal_data(
  df_emc,
  exposure = "Xstar",
  outcome = "Y",
  confounders = c("C1", "C2", "C3"),
)

print(df_observed)

df_validation <- causal_data(
  df_emc_source,
  exposure = "X",
  outcome = "Y",
  confounders = c("C1", "C2", "C3"),
  misclassified_exposure = "Xstar"
)

print(df_validation)

adjust_em(df_observed, df_validation)

# ----

adjust_om <- function(
  observed_data,
  validation_data,
  level = 0.95
) {

  if (all(observed_data$confounders %in% validation_data$confounders) == FALSE) {
    stop("All confounders in observed data must be present in validation data.")
  }

  if (is.null(validation_data$misclassified_outcome)) {
    stop("This function is adjusting for a misclassified outcome\nValidation data must have an outcome and misclassified outcome specified.")
  }

  n <- nrow(observed_data$data)

  df <- data.frame(
    X = observed_data$data[, observed_data$exposure],
    Ystar = observed_data$data[, observed_data$outcome]
  )
  df <- bind_cols(df, observed_data$data[, observed_data$confounders])

  if (!all(df$Ystar %in% 0:1)) {
    stop("Observed, misclassified outcome must be a binary integer.")
  }

  df_val <- data.frame(
    X = validation_data$data[, validation_data$exposure],
    Y = validation_data$data[, validation_data$outcome],
    Ystar = validation_data$data[, validation_data$misclassified_outcome]
  )

  if (all(df$X %in% 0:1)) {
    if (!all(df_val$X %in% 0:1)) {
      stop("Exposures from both datasets must all be binary or all be continuous.")
    }
  }

  if (all(df$Y %in% 0:1)) {
    if (!all(df_val$Y %in% 0:1) || !all(df_val$Ystar %in% 0:1)) {
      stop("Outcomes from both datasets must all be binary or all be continuous.")
    }
  }

  if (!all(df_val$Ystar %in% 0:1)) {
    stop("Validation, misclassified outcome must be a binary integer.")
  }
  if (!all(df_val$Y %in% 0:1)) {
    stop("Validation, true outcome must be a binary integer.")
  }

  y_mod <- glm(Y ~ X + Ystar + .,
               family = binomial(link = "logit"),
               data = df_val)

  # sequence along coefs to predict Y in observed data
  y_mod_coefs <- coef(y_mod)
  y_pred <- y_mod_coefs[1]

  for (i in 2:length(y_mod_coefs)) {
    var_name <- names(y_mod_coefs)[i]
    y_pred <- y_pred + df[[var_name]] * y_mod_coefs[i]
  }

  df$Ypred <- rbinom(n, 1, plogis(y_pred))

  final <- glm(
    Ypred ~ X + .,
    family = binomial(link = "logit"),
    data = df
  )

  est <- summary(final)$coef[2, 1]
  se <- summary(final)$coef[2, 2]
  alpha <- 1 - level


  estimate <- exp(est)
  ci <- c(exp(est + se * qnorm(alpha / 2)),
          exp(est + se * qnorm(1 - alpha / 2)))

  return(list(estimate = estimate, ci = ci))

}

df_observed <- causal_data(
  df_omc,
  exposure = "X",
  outcome = "Ystar",
  confounders = c("C1", "C2", "C3"),
)

print(df_observed)

df_validation <- causal_data(
  df_omc_source,
  exposure = "X",
  outcome = "Y",
  confounders = c("C1", "C2", "C3"),
  misclassified_outcome = "Ystar"
)

print(df_validation)

adjust_om(df_observed, df_validation)
