# Adjust for outcome misclassification and selection bias

# the following functions feed into adjust_om_sel():
# adjust_om_sel_val() (data_validation input),
# adjust_om_sel_coef() (bias_params input)

adjust_om_sel_val <- function(
    data_observed,
    data_validation) {
  if (!all(data_observed$confounders %in% data_validation$confounders)) {
    stop(
      "All confounders in observed data must be present in validation data.",
      call. = FALSE
    )
  }

  if (
    is.null(data_validation$misclassified_outcome) ||
      is.null(data_validation$selection)
  ) {
    stop(
      paste0(
        "This function is adjusting for a misclassified outcome and selection bias.",
        "\n",
        "Validation data must have: 1) a true and misclassified outcome specified, 2) a selection indicator column specified."
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
    Ystar = data_validation$data[, data_validation$misclassified_outcome],
    S = data_validation$data[, data_validation$selection]
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
  force_binary(
    df_val$S,
    "Selection indicator in validation data must be a binary integer."
  )

  y_mod <- glm(Y ~ X + Ystar + . - S,
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

  s_mod <- glm(S ~ X + Ystar + . - Y,
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

  suppressWarnings({
    final <- glm(
      Ypred ~ X + . - Ystar - Spred,
      family = binomial(link = "logit"),
      weights = (1 / df$Spred),
      data = df
    )
  })

  return(final)
}


adjust_om_sel_coef <- function(
    data_observed,
    y_model_coefs,
    s_model_coefs) {
  data <- data_observed$data
  n <- nrow(data)
  confounders <- data_observed$confounders
  len_c <- length(confounders)
  len_y_coefs <- length(y_model_coefs)
  len_s_coefs <- length(s_model_coefs)

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
  force_len(
    len_s_coefs,
    3 + len_c,
    paste0(
      "Incorrect length of S model coefficients. ",
      "Length should equal 3 + number of confounders."
    )
  )

  s1_0 <- s_model_coefs[1]
  s1_x <- s_model_coefs[2]
  s1_ystar <- s_model_coefs[3]

  y1_0 <- y_model_coefs[1]
  y1_x <- y_model_coefs[2]
  y1_ystar <- y_model_coefs[3]

  if (is.null(confounders)) {
    df <- data.frame(X = x, Ystar = ystar)

    y1_pred <- plogis(y1_0 + y1_x * x + y1_ystar * ystar)
    y1_pred <- rep(y1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(
        Ybar = rep(c(1, 0), each = n),
        pS = plogis(s1_0 + s1_x * .data$X + s1_ystar * .data$Ystar),
        pY = case_when(
          Ybar == 1 ~ y1_pred,
          Ybar == 0 ~ 1 - y1_pred
        )
      )

    suppressWarnings({
      final <- glm(
        Ybar ~ X,
        family = binomial(link = "logit"),
        weights = (combined$pY / combined$pS),
        data = combined
      )
    })
  } else if (len_c == 1) {
    c1 <- data[, confounders]
    df <- data.frame(X = x, Ystar = ystar, C1 = c1)
    y1_c1 <- y_model_coefs[4]
    s1_c1 <- s_model_coefs[4]

    y1_pred <- plogis(y1_0 + y1_x * x + y1_ystar * ystar + y1_c1 * c1)
    y1_pred <- rep(y1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(
        Ybar = rep(c(1, 0), each = n),
        pS = plogis(
          s1_0 + s1_x * .data$X + s1_ystar * .data$Ystar +
            s1_c1 * .data$C1
        ),
        pY = case_when(
          Ybar == 1 ~ y1_pred,
          Ybar == 0 ~ 1 - y1_pred
        )
      )

    suppressWarnings({
      final <- glm(
        Ybar ~ X + C1,
        family = binomial(link = "logit"),
        weights = (combined$pY / combined$pS),
        data = combined
      )
    })
  } else if (len_c == 2) {
    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]

    df <- data.frame(X = x, Ystar = ystar, C1 = c1, C2 = c2)

    s1_c1 <- s_model_coefs[4]
    s1_c2 <- s_model_coefs[5]

    y1_c1 <- y_model_coefs[4]
    y1_c2 <- y_model_coefs[5]

    y1_pred <- plogis(
      y1_0 + y1_x * x +
        y1_ystar * ystar + y1_c1 * c1 + y1_c2 * c2
    )
    y1_pred <- rep(y1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(
        Ybar = rep(c(1, 0), each = n),
        pS = plogis(
          s1_0 + s1_x * .data$X + s1_ystar * .data$Ystar +
            s1_c1 * .data$C1 + s1_c2 * .data$C2
        ),
        pY = case_when(
          Ybar == 1 ~ y1_pred,
          Ybar == 0 ~ 1 - y1_pred
        )
      )

    suppressWarnings({
      final <- glm(
        Ybar ~ X + C1 + C2,
        family = binomial(link = "logit"),
        weights = (combined$pY / combined$pS),
        data = combined
      )
    })
  } else if (len_c == 3) {
    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]
    c3 <- data[, confounders[3]]

    df <- data.frame(X = x, Ystar = ystar, C1 = c1, C2 = c2, C3 = c3)

    s1_c1 <- s_model_coefs[4]
    s1_c2 <- s_model_coefs[5]
    s1_c3 <- s_model_coefs[6]

    y1_c1 <- y_model_coefs[4]
    y1_c2 <- y_model_coefs[5]
    y1_c3 <- y_model_coefs[6]

    y1_pred <- plogis(
      y1_0 + y1_x * x + y1_ystar * ystar + y1_c1 * c1 + y1_c2 * c2 + y1_c3 * c3
    )
    y1_pred <- rep(y1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(
        Ybar = rep(c(1, 0), each = n),
        pS = plogis(
          s1_0 + s1_x * .data$X + s1_ystar * .data$Ystar +
            s1_c1 * .data$C1 + s1_c2 * .data$C2 + s1_c3 * .data$C3
        ),
        pY = case_when(
          Ybar == 1 ~ y1_pred,
          Ybar == 0 ~ 1 - y1_pred
        )
      )

    suppressWarnings({
      final <- glm(
        Ybar ~ X + C1 + C2 + C3,
        family = binomial(link = "logit"),
        weights = (combined$pY / combined$pS),
        data = combined
      )
    })
  } else if (len_c > 3) {
    stop(
      "This function is currently not compatible with >3 confounders.",
      call. = FALSE
    )
  }

  return(final)
}


adjust_om_sel <- function(
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

  if (!is.null(data_validation)) {
    final <- adjust_om_sel_val(
      data_observed,
      data_validation
    )
  } else if (!is.null(bias_params)) {
    if (is.null(bias_params$coef_list$y) && is.null(bias_params$coef_list$s)) {
      stop(
        paste0(
          "bias_params must specify parameters for outcome ",
          "misclassification and selection bias"
        ),
        call. = FALSE
      )
    }
    final <- adjust_om_sel_coef(
      data_observed,
      bias_params$coef_list$y,
      bias_params$coef_list$s
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
