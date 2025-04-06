# Adjust for exposure misclassification and selection bias

# the following functions feed into adjust_em_sel():
# adjust_em_sel_val() (data_validation input),
# adjust_em_sel_coef() (bias_params input)

adjust_em_sel_val <- function(
    data_observed,
    data_validation) {
  if (!all(data_observed$confounders %in% data_validation$confounders)) {
    stop(
      "All confounders in observed data must be present in validation data.",
      call. = FALSE
    )
  }

  if (
    is.null(data_validation$misclassified_exposure) ||
      is.null(data_validation$selection)
  ) {
    stop(
      paste0(
        "Attempting to adjust for a misclassified exposure and selection bias.",
        "\n",
        "Validation data must have: 1) a true and misclassified exposure specified, 2) a selection indicator column specified."
      ),
      call. = FALSE
    )
  }

  n <- nrow(data_observed$data)

  df <- data.frame(
    Xstar = data_observed$data[, data_observed$exposure],
    Y = data_observed$data[, data_observed$outcome]
  )
  df <- bind_cols(
    df,
    data_observed$data %>%
      select(all_of(data_observed$confounders))
  )

  if (all(df$Y %in% 0:1)) {
    y_binary <- TRUE
  } else {
    y_binary <- FALSE
  }

  df_val <- data.frame(
    X = data_validation$data[, data_validation$true_exposure],
    Y = data_validation$data[, data_validation$true_outcome],
    Xstar = data_validation$data[, data_validation$misclassified_exposure],
    S = data_validation$data[, data_validation$selection]
  )
  df_val <- bind_cols(
    df_val,
    data_validation$data %>%
      select(all_of(data_validation$confounders))
  )

  force_match(
    df$Y,
    df_val$Y,
    "Outcomes from both datasets must both be binary or both be continuous."
  )
  force_binary(
    df$Xstar,
    "Exposure in observed data must be a binary integer."
  )
  force_binary(
    df_val$Xstar,
    "Misclassified exposure in validation data must be a binary integer."
  )
  force_binary(
    df_val$X,
    "True exposure in validation data must be a binary integer."
  )
  force_binary(
    df_val$S,
    "Selection indicator in validation data must be a binary integer."
  )

  x_mod <- glm(X ~ Xstar + Y + . - S,
    family = binomial(link = "logit"),
    data = df_val
  )

  x_mod_coefs <- coef(x_mod)
  x_pred <- x_mod_coefs[1]

  for (i in 2:length(x_mod_coefs)) {
    var_name <- names(x_mod_coefs)[i]
    x_pred <- x_pred + df[[var_name]] * x_mod_coefs[i]
  }

  df$Xpred <- rbinom(n, 1, plogis(x_pred))

  s_mod <- glm(S ~ Xstar + Y + . - X,
    family = binomial(link = "logit"),
    data = df_val
  )

  s_mod_coefs <- coef(s_mod)
  s_pred <- s_mod_coefs[1]

  for (i in 2:length(s_mod_coefs)) {
    var_name <- names(s_mod_coefs)[i]
    s_pred <- s_pred + df[[var_name]] * s_mod_coefs[i]
  }

  df$Spred <- plogis(s_pred)

  if (y_binary) {
    suppressWarnings({
      final <- glm(
        Y ~ Xpred + . - Xstar - Spred,
        family = binomial(link = "logit"),
        weights = (1 / df$Spred),
        data = df
      )
    })
  } else {
    suppressWarnings({
      final <- lm(
        Y ~ Xpred + . - Xstar - Spred,
        weights = (1 / df$Spred),
        data = df
      )
    })
  }

  return(final)
}


adjust_em_sel_coef <- function(
    data_observed,
    x_model_coefs,
    s_model_coefs,
    level = 0.95) {
  data <- data_observed$data
  n <- nrow(data)
  confounders <- data_observed$confounders
  len_c <- length(confounders)
  len_x_coefs <- length(x_model_coefs)
  len_s_coefs <- length(s_model_coefs)

  xstar <- data[, data_observed$exposure]
  y <- data[, data_observed$outcome]

  force_binary(xstar, "Exposure must be a binary integer.")
  force_len(
    len_x_coefs,
    3 + len_c,
    paste0(
      "Incorrect length of X model coefficients. ",
      "Length should equal 3 + number of confounders."
    )
  )
  force_len(
    len_s_coefs,
    3 + len_c,
    paste0(
      "Incorrect length of S model coefficients. ",
      "Length should equal 3 + number of confounders."
    )
  )

  if (all(y %in% 0:1)) {
    y_binary <- TRUE
  } else {
    y_binary <- FALSE
  }

  s1_0 <- s_model_coefs[1]
  s1_xstar <- s_model_coefs[2]
  s1_y <- s_model_coefs[3]

  x1_0 <- x_model_coefs[1]
  x1_xstar <- x_model_coefs[2]
  x1_y <- x_model_coefs[3]

  if (is.null(confounders)) {
    df <- data.frame(Xstar = xstar, Y = y)

    x1_pred <- plogis(x1_0 + x1_xstar * xstar + x1_y * y)
    x1_pred <- rep(x1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(
        Xbar = rep(c(1, 0), each = n),
        pS = plogis(s1_0 + s1_xstar * .data$Xstar + s1_y * .data$Y),
        pX = case_when(
          Xbar == 1 ~ x1_pred,
          Xbar == 0 ~ 1 - x1_pred
        )
      )

    if (y_binary) {
      suppressWarnings({
        final <- glm(
          Y ~ Xbar,
          family = binomial(link = "logit"),
          weights = (combined$pX / combined$pS),
          data = combined
        )
      })
    } else {
      suppressWarnings({
        final <- lm(
          Y ~ Xbar,
          weights = (combined$pX / combined$pS),
          data = combined
        )
      })
    }
  } else if (len_c == 1) {
    c1 <- data[, confounders]
    df <- data.frame(Xstar = xstar, Y = y, C1 = c1)
    x1_c1 <- x_model_coefs[4]
    s1_c1 <- s_model_coefs[4]

    x1_pred <- plogis(x1_0 + x1_xstar * xstar + x1_y * y + x1_c1 * c1)
    x1_pred <- rep(x1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(
        Xbar = rep(c(1, 0), each = n),
        pS = plogis(
          s1_0 + s1_xstar * .data$Xstar + s1_y * .data$Y +
            s1_c1 * .data$C1
        ),
        pX = case_when(
          Xbar == 1 ~ x1_pred,
          Xbar == 0 ~ 1 - x1_pred
        )
      )

    if (y_binary) {
      suppressWarnings({
        final <- glm(
          Y ~ Xbar + C1,
          family = binomial(link = "logit"),
          weights = (combined$pX / combined$pS),
          data = combined
        )
      })
    } else {
      suppressWarnings({
        final <- lm(
          Y ~ Xbar + C1,
          weights = (combined$pX / combined$pS),
          data = combined
        )
      })
    }
  } else if (len_c == 2) {
    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]

    df <- data.frame(Xstar = xstar, Y = y, C1 = c1, C2 = c2)

    s1_c1 <- s_model_coefs[4]
    s1_c2 <- s_model_coefs[5]

    x1_c1 <- x_model_coefs[4]
    x1_c2 <- x_model_coefs[5]

    x1_pred <- plogis(
      x1_0 + x1_xstar * xstar +
        x1_y * y + x1_c1 * c1 + x1_c2 * c2
    )
    x1_pred <- rep(x1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(
        Xbar = rep(c(1, 0), each = n),
        pS = plogis(
          s1_0 + s1_xstar * .data$Xstar + s1_y * .data$Y +
            s1_c1 * .data$C1 + s1_c2 * .data$C2
        ),
        pX = case_when(
          Xbar == 1 ~ x1_pred,
          Xbar == 0 ~ 1 - x1_pred
        )
      )

    if (y_binary) {
      suppressWarnings({
        final <- glm(
          Y ~ Xbar + C1 + C2,
          family = binomial(link = "logit"),
          weights = (combined$pX / combined$pS),
          data = combined
        )
      })
    } else {
      suppressWarnings({
        final <- lm(
          Y ~ Xbar + C1 + C2,
          weights = (combined$pX / combined$pS),
          data = combined
        )
      })
    }
  } else if (len_c == 3) {
    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]
    c3 <- data[, confounders[3]]

    df <- data.frame(Xstar = xstar, Y = y, C1 = c1, C2 = c2, C3 = c3)

    s1_c1 <- s_model_coefs[4]
    s1_c2 <- s_model_coefs[5]
    s1_c3 <- s_model_coefs[6]

    x1_c1 <- x_model_coefs[4]
    x1_c2 <- x_model_coefs[5]
    x1_c3 <- x_model_coefs[6]

    x1_pred <- plogis(
      x1_0 + x1_xstar * xstar + x1_y * y + x1_c1 * c1 + x1_c2 * c2 + x1_c3 * c3
    )
    x1_pred <- rep(x1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(
        Xbar = rep(c(1, 0), each = n),
        pS = plogis(
          s1_0 + s1_xstar * .data$Xstar + s1_y * .data$Y +
            s1_c1 * .data$C1 + s1_c2 * .data$C2 + s1_c3 * .data$C3
        ),
        pX = case_when(
          Xbar == 1 ~ x1_pred,
          Xbar == 0 ~ 1 - x1_pred
        )
      )

    if (y_binary) {
      suppressWarnings({
        final <- glm(
          Y ~ Xbar + C1 + C2 + C3,
          family = binomial(link = "logit"),
          weights = (combined$pX / combined$pS),
          data = combined
        )
      })
    } else {
      suppressWarnings({
        final <- lm(
          Y ~ Xbar + C1 + C2 + C3,
          weights = (combined$pX / combined$pS),
          data = combined
        )
      })
    }
  } else if (len_c > 3) {
    stop(
      "This function is currently not compatible with >3 confounders.",
      call. = FALSE
    )
  }

  return(final)
}


adjust_em_sel <- function(
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
  xstar <- data[, data_observed$exposure]
  y <- data[, data_observed$outcome]

  if (all(y %in% 0:1)) {
    y_binary <- TRUE
  } else {
    y_binary <- FALSE
  }

  if (!is.null(data_validation)) {
    final <- adjust_em_sel_val(
      data_observed,
      data_validation
    )
  } else if (!is.null(bias_params)) {
    if (is.null(bias_params$coef_list$x) && is.null(bias_params$coef_list$s)) {
      stop(
        paste0(
          "bias_params must specify parameters for exposure ",
          "misclassification and selection bias"
        ),
        call. = FALSE
      )
    }
    final <- adjust_em_sel_coef(
      data_observed,
      bias_params$coef_list$x,
      bias_params$coef_list$s
    )
  }

  est <- summary(final)$coef[2, 1]
  se <- summary(final)$coef[2, 2]
  alpha <- 1 - level

  if (y_binary) {
    estimate <- exp(est)
    ci <- c(
      exp(est + se * qnorm(alpha / 2)),
      exp(est + se * qnorm(1 - alpha / 2))
    )
  } else {
    estimate <- est
    ci <- c(
      est + se * qnorm(alpha / 2),
      est + se * qnorm(1 - alpha / 2)
    )
  }

  return(list(estimate = estimate, ci = ci))
}
