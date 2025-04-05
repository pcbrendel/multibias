# Adjust for outcome misclassification

# the following functions feed into adjust_om():
# adjust_om_val() (data_validation input),
# adjust_om_coef() (bias_params input)

adjust_om_val <- function(
    data_observed,
    data_validation) {
  if (!all(data_observed$confounders %in% data_validation$confounders)) {
    stop(
      "All confounders in observed data must be present in validation data.",
      call. = FALSE
    )
  }

  if (is.null(data_validation$misclassified_outcome)) {
    stop(
      paste0(
        "Attempting to adjust for a misclassified outcome.",
        "\n",
        "Validation data must have a true and misclassified outcome specified."
      ),
      call. = FALSE
    )
  }

  n <- nrow(data_observed$data)

  df <- data.frame(
    X = data_observed$data[, data_observed$exposure],
    Ystar = data_observed$data[, data_observed$outcome]
  )
  df <- bind_cols(
    df,
    data_observed$data %>%
      select(all_of(data_observed$confounders))
  )

  df_val <- data.frame(
    X = data_validation$data[, data_validation$true_exposure],
    Y = data_validation$data[, data_validation$true_outcome],
    Ystar = data_validation$data[, data_validation$misclassified_outcome]
  )
  df_val <- bind_cols(
    df_val,
    data_validation$data %>%
      select(all_of(data_validation$confounders))
  )

  force_match(
    df$X,
    df_val$X,
    "Exposures from both datasets must both be binary or both be continuous."
  )
  force_binary(
    df$Ystar,
    "Outcome in observed data must be a binary integer."
  )
  force_binary(
    df_val$Ystar,
    "Misclassified outcome in validation data must be a binary integer."
  )
  force_binary(
    df_val$Y,
    "True outcome in validation data must be a binary integer."
  )

  y_mod <- glm(Y ~ X + Ystar + .,
    family = binomial(link = "logit"),
    data = df_val
  )

  y_mod_coefs <- coef(y_mod)
  y_pred <- y_mod_coefs[1]

  for (i in 2:length(y_mod_coefs)) {
    var_name <- names(y_mod_coefs)[i]
    y_pred <- y_pred + df[[var_name]] * y_mod_coefs[i]
  }

  df$Ypred <- rbinom(n, 1, plogis(y_pred))

  final <- glm(
    Ypred ~ X + . - Ystar,
    family = binomial(link = "logit"),
    data = df
  )

  return(final)
}


adjust_om_coef <- function(
    data_observed,
    y_model_coefs) {
  data <- data_observed$data
  n <- nrow(data)
  confounders <- data_observed$confounders
  len_c <- length(confounders)
  len_y_coefs <- length(y_model_coefs)

  x <- data[, data_observed$exposure]
  ystar <- data[, data_observed$outcome]

  force_binary(ystar, "Outcome must be a binary integer.")
  force_len(
    len_y_coefs,
    3 + len_c,
    paste0(
      "Incorrect length of Y model coefficients. ",
      "Length should equal 3 + number of confounders."
    )
  )

  y1_0 <- y_model_coefs[1]
  y1_x <- y_model_coefs[2]
  y1_ystar <- y_model_coefs[3]

  if (is.null(confounders)) {
    df <- data.frame(X = x, Ystar = ystar)
    df$Ypred <- rbinom(n, 1, plogis(y1_0 + y1_x * df$X + y1_ystar * df$Ystar))

    final <- glm(
      Ypred ~ X,
      family = binomial(link = "logit"),
      data = df
    )
  } else if (len_c == 1) {
    c1 <- data[, confounders]
    df <- data.frame(X = x, Ystar = ystar, C1 = c1)

    y1_c1 <- y_model_coefs[4]

    df$Ypred <- rbinom(
      n, 1, plogis(
        y1_0 + y1_x * df$X + y1_ystar * df$Ystar +
          y1_c1 * df$C1
      )
    )

    final <- glm(
      Ypred ~ X + C1,
      family = binomial(link = "logit"),
      data = df
    )
  } else if (len_c == 2) {
    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]

    df <- data.frame(X = x, Ystar = ystar, C1 = c1, C2 = c2)

    y1_c1 <- y_model_coefs[4]
    y1_c2 <- y_model_coefs[5]

    df$Ypred <- rbinom(
      n, 1, plogis(
        y1_0 + y1_x * df$X + y1_ystar * df$Ystar +
          y1_c1 * df$C1 + y1_c2 * df$C2
      )
    )

    final <- glm(
      Ypred ~ X + C1 + C2,
      family = binomial(link = "logit"),
      data = df
    )
  } else if (len_c == 3) {
    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]
    c3 <- data[, confounders[3]]

    df <- data.frame(X = x, Ystar = ystar, C1 = c1, C2 = c2, C3 = c3)

    y1_c1 <- y_model_coefs[4]
    y1_c2 <- y_model_coefs[5]
    y1_c3 <- y_model_coefs[6]

    df$Ypred <- rbinom(
      n, 1, plogis(
        y1_0 + y1_x * df$X + y1_ystar * df$Ystar +
          y1_c1 * df$C1 + y1_c2 * df$C2 + y1_c3 * df$C3
      )
    )

    final <- glm(
      Ypred ~ X + C1 + C2 + C3,
      family = binomial(link = "logit"),
      data = df
    )
  } else if (len_c > 3) {
    stop(
      "This function is currently not compatible with >3 confounders.",
      call. = FALSE
    )
  }

  return(final)
}


adjust_om <- function(
    data_observed,
    data_validation = NULL,
    bias_params = NULL,
    level = 0.95) {
  if (
    (!is.null(data_validation) && !is.null(bias_params)) ||
      (is.null(data_validation) && is.null(bias_params))
  ) {
    stop(
      "One of data_validation or bias_params must be non-null.",
      call. = FALSE
    )
  }
  data <- data_observed$data

  if (!is.null(data_validation)) {
    final <- adjust_om_val(
      data_observed,
      data_validation
    )
  } else if (!is.null(bias_params)) {
    if (is.null(bias_params$coef_list$y)) {
      stop(
        "bias_params must specify parameters for y to adjust for outcome misclassification",
        call. = FALSE
      )
    }
    final <- adjust_om_coef(
      data_observed,
      y_model_coefs = bias_params$coef_list$y
    )
  }

  est <- summary(final)$coef[2, 1]
  se <- summary(final)$coef[2, 2]
  alpha <- 1 - level

  estimate <- exp(est)
  ci <- c(
    exp(est + se * qnorm(alpha / 2)),
    exp(est + se * qnorm(1 - alpha / 2))
  )

  return(list(estimate = estimate, ci = ci))
}
