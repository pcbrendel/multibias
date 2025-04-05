# Adjust for uncontrolled confounding and exposure misclassification

# the following functions feed into adjust_uc_em():
# adjust_uc_em_val() (data_validation input),
# adjust_uc_em_coef() (bias_params input)

adjust_uc_em_val <- function(
    data_observed,
    data_validation) {
  if (!all(data_observed$confounders %in% data_validation$confounders)) {
    stop(
      "All confounders in observed data must be present in validation data.",
      call. = FALSE
    )
  }

  if (
    (length(data_validation$confounders) - length(data_observed$confounders) != 1) ||
      (is.null(data_validation$misclassified_exposure))
  ) {
    stop(
      paste0(
        "Attempting to adjust for unobserved confounding from one confounder and exposure misclassification.",
        "\n",
        "Validation data must have: 1) one more confounder than the observed data, 2) a true and misclassified exposure specified."
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
    Xstar = data_validation$data[, data_validation$misclassified_exposure]
  )

  uc <- setdiff(data_validation$confounders, data_observed$confounders)
  df_val$U <- data_validation$data[, uc]
  df_val <- bind_cols(
    df_val,
    data_validation$data %>%
      select(all_of(data_observed$confounders))
  )

  force_match(
    df$Y,
    df_val$Y,
    "Outcomes from both datasets must both be binary or both be continuous."
  )
  force_binary(
    df_val$U,
    "Uncontrolled confounder in validation data must be a binary integer."
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

  x_mod <- glm(X ~ Xstar + Y + . - U,
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

  u_mod <- glm(U ~ X + Y,
    family = binomial(link = "logit"),
    data = df_val
  )

  u_mod_coefs <- coef(u_mod)
  u_pred <- u_mod_coefs[1]

  for (i in 2:length(u_mod_coefs)) {
    var_name <- names(u_mod_coefs)[i]
    var_name <- gsub("X", "Xpred", var_name) # col X is not in df
    u_pred <- u_pred + df[[var_name]] * u_mod_coefs[i]
  }

  df$Upred <- rbinom(n, 1, plogis(u_pred))

  if (y_binary) {
    final <- glm(
      Y ~ Xpred + Upred + . - Xstar,
      family = binomial(link = "logit"),
      data = df
    )
  } else {
    final <- lm(
      Y ~ Xpred + Upred + . - Xstar,
      data = df
    )
  }

  return(final)
}


adjust_uc_em_coef_single <- function(
    data_observed,
    u_model_coefs,
    x_model_coefs) {
  data <- data_observed$data
  n <- nrow(data)
  confounders <- data_observed$confounders
  len_c <- length(confounders)
  len_u_coefs <- length(u_model_coefs)
  len_x_coefs <- length(x_model_coefs)

  xstar <- data[, data_observed$exposure]
  y <- data[, data_observed$outcome]

  force_binary(xstar, "Exposure must be a binary integer.")
  force_len(
    len_u_coefs,
    3,
    paste0(
      "Incorrect length of U model coefficients. ",
      "Length should equal 3."
    )
  )
  force_len(
    len_x_coefs,
    3 + len_c,
    paste0(
      "Incorrect length of X model coefficients. ",
      "Length should equal 3 + number of confounders."
    )
  )

  if (all(y %in% 0:1)) {
    y_binary <- TRUE
  } else {
    y_binary <- FALSE
  }

  u1_0 <- u_model_coefs[1]
  u1_x <- u_model_coefs[2]
  u1_y <- u_model_coefs[3]

  x1_0 <- x_model_coefs[1]
  x1_xstar <- x_model_coefs[2]
  x1_y <- x_model_coefs[3]

  if (is.null(confounders)) {
    df <- data.frame(Xstar = xstar, Y = y)
    df$Xpred <- rbinom(n, 1, plogis(x1_0 + x1_xstar * df$Xstar + x1_y * df$Y))
    df$Upred <- rbinom(n, 1, plogis(u1_0 + u1_x * df$Xpred + u1_y * df$Y))

    if (y_binary) {
      final <- glm(
        Y ~ Xpred + Upred,
        family = binomial(link = "logit"),
        data = df
      )
    } else {
      final <- lm(
        Y ~ Xpred + Upred,
        data = df
      )
    }
  } else if (len_c == 1) {
    c1 <- data[, confounders]
    df <- data.frame(Xstar = xstar, Y = y, C1 = c1)

    x1_c1 <- x_model_coefs[4]

    df$Xpred <- rbinom(
      n, 1, plogis(
        x1_0 + x1_xstar * df$Xstar +
          x1_y * df$Y + x1_c1 * df$C1
      )
    )
    df$Upred <- rbinom(
      n, 1, plogis(
        u1_0 + u1_x * df$Xpred + u1_y * df$Y
      )
    )

    if (y_binary) {
      final <- glm(
        Y ~ Xpred + C1 + Upred,
        family = binomial(link = "logit"),
        data = df
      )
    } else {
      final <- lm(
        Y ~ Xpred + C1 + Upred,
        data = df
      )
    }
  } else if (len_c == 2) {
    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]

    df <- data.frame(Xstar = xstar, Y = y, C1 = c1, C2 = c2)

    x1_c1 <- x_model_coefs[4]
    x1_c2 <- x_model_coefs[5]

    df$Xpred <- rbinom(
      n, 1, plogis(
        x1_0 + x1_xstar * df$Xstar + x1_y * df$Y +
          x1_c1 * df$C1 + x1_c2 * df$C2
      )
    )
    df$Upred <- rbinom(
      n, 1, plogis(u1_0 + u1_x * df$Xpred + u1_y * df$Y)
    )

    if (y_binary) {
      final <- glm(
        Y ~ Xpred + C1 + C2 + Upred,
        family = binomial(link = "logit"),
        data = df
      )
    } else {
      final <- lm(
        Y ~ Xpred + C1 + C2 + Upred,
        data = df
      )
    }
  } else if (len_c == 3) {
    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]
    c3 <- data[, confounders[3]]

    df <- data.frame(Xstar = xstar, Y = y, C1 = c1, C2 = c2, C3 = c3)

    x1_c1 <- x_model_coefs[4]
    x1_c2 <- x_model_coefs[5]
    x1_c3 <- x_model_coefs[6]

    df$Xpred <- rbinom(
      n, 1,
      plogis(
        x1_0 + x1_xstar * df$Xstar + x1_y * df$Y +
          x1_c1 * df$C1 + x1_c2 * df$C2 + x1_c3 * df$C3
      )
    )
    df$Upred <- rbinom(n, 1, plogis(u1_0 + u1_x * df$Xpred + u1_y * df$Y))

    if (y_binary) {
      final <- glm(
        Y ~ Xpred + C1 + C2 + C3 + Upred,
        family = binomial(link = "logit"),
        data = df
      )
    } else {
      final <- lm(
        Y ~ Xpred + C1 + C2 + C3 + Upred,
        data = df
      )
    }
  } else if (len_c > 3) {
    stop(
      "This function is currently not compatible with >3 confounders.",
      call. = FALSE
    )
  }

  return(final)
}


adjust_uc_em_coef_multinom <- function(
    data_observed,
    x1u0_model_coefs,
    x0u1_model_coefs,
    x1u1_model_coefs) {
  data <- data_observed$data
  n <- nrow(data)
  confounders <- data_observed$confounders
  len_c <- length(confounders)
  len_x1u0_coefs <- length(x1u0_model_coefs)
  len_x0u1_coefs <- length(x0u1_model_coefs)
  len_x1u1_coefs <- length(x1u1_model_coefs)

  xstar <- data[, data_observed$exposure]
  y <- data[, data_observed$outcome]

  force_binary(xstar, "Exposure must be a binary integer.")
  force_len(
    len_x1u0_coefs,
    3 + len_c,
    paste0(
      "Incorrect length of X1U0 model coefficients. ",
      "Length should equal 3 + number of confounders."
    )
  )
  force_len(
    len_x0u1_coefs,
    3 + len_c,
    paste0(
      "Incorrect length of X0U1 model coefficients. ",
      "Length should equal 3 + number of confounders."
    )
  )
  force_len(
    len_x1u1_coefs,
    3 + len_c,
    paste0(
      "Incorrect length of X1U1 model coefficients. ",
      "Length should equal 3 + number of confounders."
    )
  )

  if (all(y %in% 0:1)) {
    y_binary <- TRUE
  } else {
    y_binary <- FALSE
  }

  x1u0_0 <- x1u0_model_coefs[1]
  x1u0_xstar <- x1u0_model_coefs[2]
  x1u0_y <- x1u0_model_coefs[3]

  x0u1_0 <- x0u1_model_coefs[1]
  x0u1_xstar <- x0u1_model_coefs[2]
  x0u1_y <- x0u1_model_coefs[3]

  x1u1_0 <- x1u1_model_coefs[1]
  x1u1_xstar <- x1u1_model_coefs[2]
  x1u1_y <- x1u1_model_coefs[3]

  if (is.null(confounders)) {
    df <- data.frame(Xstar = xstar, Y = y)

    p_x1u0 <- exp(x1u0_0 + x1u0_xstar * df$Xstar + x1u0_y * df$Y)
    p_x0u1 <- exp(x0u1_0 + x0u1_xstar * df$Xstar + x0u1_y * df$Y)
    p_x1u1 <- exp(x1u1_0 + x1u1_xstar * df$Xstar + x1u1_y * df$Y)

    denom <- (1 + p_x1u0 + p_x0u1 + p_x1u1)

    x0u0_pred <- 1 / denom
    x1u0_pred <- p_x1u0 / denom
    x0u1_pred <- p_x0u1 / denom
    x1u1_pred <- p_x1u1 / denom

    df_xu_pred <- data.frame(
      X0U0 = x0u0_pred,
      X1U0 = x1u0_pred,
      X0U1 = x0u1_pred,
      X1U1 = x1u1_pred
    )
    df_xu_pred4 <- bind_rows(df_xu_pred, df_xu_pred, df_xu_pred, df_xu_pred)

    combined <- bind_rows(df, df, df, df) %>%
      bind_cols(df_xu_pred4) %>%
      mutate(
        Xbar = rep(c(1, 0, 1, 0), each = n),
        Ubar = rep(c(1, 1, 0, 0), each = n),
        pXU = case_when(
          Xbar == 0 & Ubar == 0 ~ X0U0,
          Xbar == 1 & Ubar == 0 ~ X1U0,
          Xbar == 0 & Ubar == 1 ~ X0U1,
          Xbar == 1 & Ubar == 1 ~ X1U1
        )
      )

    if (y_binary) {
      suppressWarnings({
        final <- glm(
          Y ~ Xbar + Ubar,
          family = binomial(link = "logit"),
          weights = combined$pXU,
          data = combined
        )
      })
    } else {
      suppressWarnings({
        final <- lm(
          Y ~ Xbar + Ubar,
          weights = combined$pXU,
          data = combined
        )
      })
    }
  } else if (len_c == 1) {
    c1 <- data[, confounders]

    df <- data.frame(Xstar = xstar, Y = y, C1 = c1)

    x1u0_c1 <- x1u0_model_coefs[4]
    x0u1_c1 <- x0u1_model_coefs[4]
    x1u1_c1 <- x1u1_model_coefs[4]

    p_x1u0 <- exp(
      x1u0_0 + x1u0_xstar * df$Xstar + x1u0_y * df$Y +
        x1u0_c1 * df$C1
    )
    p_x0u1 <- exp(
      x0u1_0 + x0u1_xstar * df$Xstar + x0u1_y * df$Y +
        x0u1_c1 * df$C1
    )
    p_x1u1 <- exp(
      x1u1_0 + x1u1_xstar * df$Xstar + x1u1_y * df$Y +
        x1u1_c1 * df$C1
    )

    denom <- (1 + p_x1u0 + p_x0u1 + p_x1u1)

    x0u0_pred <- 1 / denom
    x1u0_pred <- p_x1u0 / denom
    x0u1_pred <- p_x0u1 / denom
    x1u1_pred <- p_x1u1 / denom

    df_xu_pred <- data.frame(
      X0U0 = x0u0_pred,
      X1U0 = x1u0_pred,
      X0U1 = x0u1_pred,
      X1U1 = x1u1_pred
    )
    df_xu_pred4 <- bind_rows(df_xu_pred, df_xu_pred, df_xu_pred, df_xu_pred)

    combined <- bind_rows(df, df, df, df) %>%
      bind_cols(df_xu_pred4) %>%
      mutate(
        Xbar = rep(c(1, 0, 1, 0), each = n),
        Ubar = rep(c(1, 1, 0, 0), each = n),
        pXU = case_when(
          Xbar == 0 & Ubar == 0 ~ X0U0,
          Xbar == 1 & Ubar == 0 ~ X1U0,
          Xbar == 0 & Ubar == 1 ~ X0U1,
          Xbar == 1 & Ubar == 1 ~ X1U1
        )
      )

    if (y_binary) {
      suppressWarnings({
        final <- glm(
          Y ~ Xbar + C1 + Ubar,
          family = binomial(link = "logit"),
          weights = combined$pXU,
          data = combined
        )
      })
    } else {
      suppressWarnings({
        final <- lm(
          Y ~ Xbar + C1 + Ubar,
          weights = combined$pXU,
          data = combined
        )
      })
    }
  } else if (len_c == 2) {
    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]

    df <- data.frame(Xstar = xstar, Y = y, C1 = c1, C2 = c2)

    x1u0_c1 <- x1u0_model_coefs[4]
    x1u0_c2 <- x1u0_model_coefs[5]

    x0u1_c1 <- x0u1_model_coefs[4]
    x0u1_c2 <- x0u1_model_coefs[5]

    x1u1_c1 <- x1u1_model_coefs[4]
    x1u1_c2 <- x1u1_model_coefs[5]

    p_x1u0 <- exp(
      x1u0_0 + x1u0_xstar * df$Xstar + x1u0_y * df$Y +
        x1u0_c1 * df$C1 + x1u0_c2 * df$C2
    )
    p_x0u1 <- exp(
      x0u1_0 + x0u1_xstar * df$Xstar + x0u1_y * df$Y +
        x0u1_c1 * df$C1 + x0u1_c2 * df$C2
    )
    p_x1u1 <- exp(
      x1u1_0 + x1u1_xstar * df$Xstar + x1u1_y * df$Y +
        x1u1_c1 * df$C1 + x1u1_c2 * df$C2
    )

    denom <- (1 + p_x1u0 + p_x0u1 + p_x1u1)

    x0u0_pred <- 1 / denom
    x1u0_pred <- p_x1u0 / denom
    x0u1_pred <- p_x0u1 / denom
    x1u1_pred <- p_x1u1 / denom

    df_xu_pred <- data.frame(
      X0U0 = x0u0_pred,
      X1U0 = x1u0_pred,
      X0U1 = x0u1_pred,
      X1U1 = x1u1_pred
    )
    df_xu_pred4 <- bind_rows(df_xu_pred, df_xu_pred, df_xu_pred, df_xu_pred)

    combined <- bind_rows(df, df, df, df) %>%
      bind_cols(df_xu_pred4) %>%
      mutate(
        Xbar = rep(c(1, 0, 1, 0), each = n),
        Ubar = rep(c(1, 1, 0, 0), each = n),
        pXU = case_when(
          Xbar == 0 & Ubar == 0 ~ X0U0,
          Xbar == 1 & Ubar == 0 ~ X1U0,
          Xbar == 0 & Ubar == 1 ~ X0U1,
          Xbar == 1 & Ubar == 1 ~ X1U1
        )
      )

    if (y_binary) {
      suppressWarnings({
        final <- glm(
          Y ~ Xbar + C1 + C2 + Ubar,
          family = binomial(link = "logit"),
          weights = combined$pXU,
          data = combined
        )
      })
    } else {
      suppressWarnings({
        final <- lm(
          Y ~ Xbar + C1 + C2 + Ubar,
          weights = combined$pXU,
          data = combined
        )
      })
    }
  } else if (len_c == 3) {
    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]
    c3 <- data[, confounders[3]]

    df <- data.frame(Xstar = xstar, Y = y, C1 = c1, C2 = c2, C3 = c3)

    x1u0_c1 <- x1u0_model_coefs[4]
    x1u0_c2 <- x1u0_model_coefs[5]
    x1u0_c3 <- x1u0_model_coefs[6]

    x0u1_c1 <- x0u1_model_coefs[4]
    x0u1_c2 <- x0u1_model_coefs[5]
    x0u1_c3 <- x0u1_model_coefs[6]

    x1u1_c1 <- x1u1_model_coefs[4]
    x1u1_c2 <- x1u1_model_coefs[5]
    x1u1_c3 <- x1u1_model_coefs[6]

    p_x1u0 <- exp(
      x1u0_0 + x1u0_xstar * df$Xstar + x1u0_y * df$Y +
        x1u0_c1 * df$C1 + x1u0_c2 * df$C2 + x1u0_c3 * df$C3
    )
    p_x0u1 <- exp(
      x0u1_0 + x0u1_xstar * df$Xstar + x0u1_y * df$Y +
        x0u1_c1 * df$C1 + x0u1_c2 * df$C2 + x0u1_c3 * df$C3
    )
    p_x1u1 <- exp(
      x1u1_0 + x1u1_xstar * df$Xstar + x1u1_y * df$Y +
        x1u1_c1 * df$C1 + x1u1_c2 * df$C2 + x1u1_c3 * df$C3
    )

    denom <- (1 + p_x1u0 + p_x0u1 + p_x1u1)

    x0u0_pred <- 1 / denom
    x1u0_pred <- p_x1u0 / denom
    x0u1_pred <- p_x0u1 / denom
    x1u1_pred <- p_x1u1 / denom

    df_xu_pred <- data.frame(
      X0U0 = x0u0_pred,
      X1U0 = x1u0_pred,
      X0U1 = x0u1_pred,
      X1U1 = x1u1_pred
    )
    df_xu_pred4 <- bind_rows(df_xu_pred, df_xu_pred, df_xu_pred, df_xu_pred)

    combined <- bind_rows(df, df, df, df) %>%
      bind_cols(df_xu_pred4) %>%
      mutate(
        Xbar = rep(c(1, 0, 1, 0), each = n),
        Ubar = rep(c(1, 1, 0, 0), each = n),
        pXU = case_when(
          Xbar == 0 & Ubar == 0 ~ X0U0,
          Xbar == 1 & Ubar == 0 ~ X1U0,
          Xbar == 0 & Ubar == 1 ~ X0U1,
          Xbar == 1 & Ubar == 1 ~ X1U1
        )
      )

    if (y_binary) {
      suppressWarnings({
        final <- glm(
          Y ~ Xbar + C1 + C2 + C3 + Ubar,
          family = binomial(link = "logit"),
          weights = combined$pXU,
          data = combined
        )
      })
    } else {
      suppressWarnings({
        final <- lm(
          Y ~ Xbar + C1 + C2 + C3 + Ubar,
          weights = combined$pXU,
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


adjust_uc_em <- function(
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
    final <- adjust_uc_em_val(
      data_observed,
      data_validation
    )
  } else if (!is.null(bias_params)) {
    if (all(c("x", "u") %in% names(bias_params$coef_list))) {
      final <- adjust_uc_em_coef_single(
        data_observed,
        u_model_coefs = bias_params$coef_list$u,
        x_model_coefs = bias_params$coef_list$x
      )
    } else if (
      all(
        c("x1u0", "x0u1", "x1u1") %in%
          names(bias_params$coef_list)
      )
    ) {
      final <- adjust_uc_em_coef_multinom(
        data_observed,
        x1u0_model_coefs = bias_params$coef_list$x1u0,
        x0u1_model_coefs = bias_params$coef_list$x0u1,
        x1u1_model_coefs = bias_params$coef_list$x1u1
      )
    } else {
      (
        stop(
          paste0(
            "bias_params must specify parameters for ",
            "exposure misclassification and uncontrolled confounding"
          ),
          call. = FALSE
        )
      )
    }
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
