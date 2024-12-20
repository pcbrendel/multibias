adjust_uc_em_val <- function(
    data_observed,
    data_validation) {
  if (!all(data_observed$confounders %in% data_validation$confounders)) {
    stop("All confounders in observed data must be present in validation data.")
  }

  if (
    length(data_validation$confounders) - length(data_observed$confounders) != 1
  ) {
    stop(
      paste0(
        "This function adjusts for unobserved confounding from one confounder.",
        "\n",
        "Validation data must have one more confounder than the observed data."
      )
    )
  }
  if (is.null(data_validation$misclassified_exposure)) {
    stop(
      paste0(
        "This function is adjusting for a misclassified exposure.",
        "\n",
        "Validation data must have a true and misclassified exposure specified."
      )
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
    stop("This function is currently not compatible with >3 confounders.")
  }

  return(final)
}


# bias adjust with multinomial coefs

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
    stop("This function is currently not compatible with >3 confounders.")
  }

  return(final)
}


#' Adust for uncontrolled confounding and exposure misclassification.
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `adjust_uc_emc()` was renamed to `adjust_uc_em()`
#' @keywords internal
#'
#' @export
adjust_uc_emc <- function(
    data_observed,
    u_model_coefs = NULL,
    x_model_coefs = NULL,
    x1u0_model_coefs = NULL,
    x0u1_model_coefs = NULL,
    x1u1_model_coefs = NULL,
    level = 0.95) {
  lifecycle::deprecate_warn(
    "1.5.3", "adjust_uc_emc()", "adjust_uc_em()"
  )
  adjust_uc_em(
    data_observed,
    u_model_coefs,
    x_model_coefs,
    x1u0_model_coefs,
    x0u1_model_coefs,
    x1u1_model_coefs,
    level
  )
}


#' Adust for uncontrolled confounding and exposure misclassification.
#'
#' `adjust_uc_em` returns the exposure-outcome odds ratio and confidence
#' interval, adjusted for uncontrolled confounding and exposure
#' misclassificaiton. Two different options for the bias parameters are
#' available here: 1) parameters from separate models of *U* and *X*
#' (`u_model_coefs` and `x_model_coefs`) or 2) parameters from a
#' joint model of *U* and *X* (`x1u0_model_coefs`,
#' `x0u1_model_coefs`, and `x1u1_model_coefs`).
#'
#' Values for the regression coefficients can be applied as
#' fixed values or as single draws from a probability
#' distribution (ex: `rnorm(1, mean = 2, sd = 1)`). The latter has
#' the advantage of allowing the researcher to capture the uncertainty
#' in the bias parameter estimates. To incorporate this uncertainty in the
#' estimate and confidence interval, this function should be run in loop across
#' bootstrap samples of the dataframe for analysis. The estimate and
#' confidence interval would then be obtained from the median and quantiles
#' of the distribution of odds ratio estimates.
#'
#' @param data_observed Object of class `data_observed` corresponding to the
#' data to perform bias analysis on.
#' @param data_validation Object of class `data_validation` corresponding to
#' the validation data used to adjust for bias in the observed data. Here, the
#' validation data should have data for the same variables as in the observed
#' data, plus data for the true and misclassified exposure corresponding to the
#' observed exposure in `data_observed`.
#' There should also be data for the confounder missing in `data_observed`.
#' @param u_model_coefs The regression coefficients corresponding to the model:
#' \ifelse{html}{\out{logit(P(U=1)) = &alpha;<sub>0</sub> + &alpha;<sub>1</sub>X + &alpha;<sub>2</sub>Y, }}{\eqn{logit(P(U=1)) = \alpha_0 + \alpha_1 X + \alpha_2 Y, }}
#' where *U* is the binary unmeasured confounder, *X* is the
#' binary true exposure, and *Y* is the outcome.
#' The number of parameters therefore equals 3.
#' @param x_model_coefs The regression coefficients corresponding to the model:
#' \ifelse{html}{\out{logit(P(X=1)) = &delta;<sub>0</sub> + &delta;<sub>1</sub>X* + &delta;<sub>2</sub>Y + &delta;<sub>2+j</sub>C<sub>j</sub>, }}{\eqn{logit(P(X=1)) = \delta_0 + \delta_1 X^* + \delta_2 Y + \delta_{2+j} C_j, }}
#' where *X* represents the binary true exposure,
#' *X** is the binary misclassified exposure, *Y* is the
#' outcome, and *C* represents the vector of measured confounders
#' (if any), and *j* corresponds to the number of measured confounders.
#' The number of parameters therefore equals 3 + *j*.
#' @param x1u0_model_coefs The regression coefficients corresponding to the
#' model:
#' \ifelse{html}{\out{log(P(X=1,U=0)/P(X=0,U=0)) = &gamma;<sub>1,0</sub> + &gamma;<sub>1,1</sub>X* + &gamma;<sub>1,2</sub>Y + &gamma;<sub>1,2+j</sub>C<sub>j</sub>, }}{\eqn{log(P(X=1,U=0)/P(X=0,U=0)) = \gamma_{1,0} + \gamma_{1,1} X^* + \gamma_{1,2} Y + \gamma_{1,2+j} C_j, }}
#' where *X* is the binary true exposure, *U* is the binary unmeasured
#' confounder, *X** is the binary misclassified exposure, *Y* is the
#' outcome, *C* represents the vector of measured confounders (if any),
#' and *j* corresponds to the number of measured confounders.
#' @param x0u1_model_coefs The regression coefficients corresponding to the
#' model:
#' \ifelse{html}{\out{log(P(X=0,U=1)/P(X=0,U=0)) = &gamma;<sub>2,0</sub> + &gamma;<sub>2,1</sub>X* + &gamma;<sub>2,2</sub>Y + &gamma;<sub>2,2+j</sub>C<sub>j</sub>, }}{\eqn{log(P(X=0,U=1)/P(X=0,U=0)) = \gamma_{2,0} + \gamma_{2,1} X^* + \gamma_{2,2} Y + \gamma_{2,2+j} C_j, }}
#' where *X* is the binary true exposure, *U* is the binary unmeasured
#' confounder, *X** is the binary misclassified exposure, *Y* is the
#' outcome, *C* represents the vector of measured confounders (if any),
#' and *j* corresponds to the number of measured confounders.
#' @param x1u1_model_coefs The regression coefficients corresponding to the
#' model:
#' \ifelse{html}{\out{log(P(X=1,U=1)/P(X=0,U=0)) = &gamma;<sub>3,0</sub> + &gamma;<sub>3,1</sub>X* + &gamma;<sub>3,2</sub>Y + &gamma;<sub>3,2+j</sub>C<sub>j</sub>, }}{\eqn{log(P(X=1,U=1)/P(X=0,U=0)) = \gamma_{3,0} + \gamma_{3,1} X^* + \gamma_{3,2} Y + \gamma_{3,2+j} C_j, }}
#' where *X* is the binary true exposure, *U* is the binary unmeasured
#' confounder, *X** is the binary misclassified exposure, *Y* is the
#' outcome, *C* represents the vector of measured confounders (if any),
#' and *j* corresponds to the number of measured confounders.
#' @param level Value from 0-1 representing the full range of the confidence
#' interval. Default is 0.95.
#'
#' @return A list where the first item is the odds ratio estimate of the
#' effect of the exposure on the outcome and the second item is the
#' confidence interval as the vector: (lower bound, upper bound).
#'
#' @examples
#' df_observed <- data_observed(
#'   data = df_uc_em,
#'   exposure = "Xstar",
#'   outcome = "Y",
#'   confounders = "C1"
#' )
#'
#' # Using validation data -----------------------------------------------------
#' df_validation <- data_validation(
#'   data = df_uc_em_source,
#'   true_exposure = "X",
#'   true_outcome = "Y",
#'   confounders = c("C1", "U"),
#'   misclassified_exposure = "Xstar",
#' )
#'
#' adjust_uc_em(
#'   data_observed = df_observed,
#'   data_validation = df_validation
#' )
#'
#' # Using u_model_coefs and x_model_coefs -------------------------------------
#' adjust_uc_em(
#'   data_observed = df_observed,
#'   u_model_coefs = c(-0.23, 0.63, 0.66),
#'   x_model_coefs = c(-2.47, 1.62, 0.73, 0.32)
#' )
#'
#' # Using x1u0_model_coefs, x0u1_model_coefs, x1u1_model_coefs ----------------
#' adjust_uc_em(
#'   data_observed = df_observed,
#'   x1u0_model_coefs = c(-2.82, 1.62, 0.68, -0.06),
#'   x0u1_model_coefs = c(-0.20, 0.00, 0.68, -0.05),
#'   x1u1_model_coefs = c(-2.36, 1.62, 1.29, 0.27)
#' )
#'
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom stats binomial
#' @importFrom stats glm
#' @importFrom stats qnorm
#' @importFrom stats rbinom
#' @importFrom stats plogis
#' @importFrom rlang .data
#'
#' @export

adjust_uc_em <- function(
    data_observed,
    data_validation = NULL,
    u_model_coefs = NULL,
    x_model_coefs = NULL,
    x1u0_model_coefs = NULL,
    x0u1_model_coefs = NULL,
    x1u1_model_coefs = NULL,
    level = 0.95) {
  if (!is.null(data_validation)) {
    if (!all(is.null(u_model_coefs), is.null(x_model_coefs), is.null(x1u0_model_coefs), is.null(x0u1_model_coefs), is.null(x1u1_model_coefs))) {
      stop("No bias parameters should be specified when 'data_validation' is used.")
    }
  } else if (!is.null(u_model_coefs) && !is.null(x_model_coefs)) {
    if (!all(is.null(data_validation), is.null(x1u0_model_coefs), is.null(x0u1_model_coefs), is.null(x1u1_model_coefs))) {
      stop("No other bias-adjusting inputs should be specified when 'u_model_coefs' and 'x_model_coefs' are used.")
    }
  } else if (!is.null(x1u0_model_coefs) && !is.null(x0u1_model_coefs) && !is.null(x1u1_model_coefs)) {
    if (!all(is.null(data_validation), is.null(u_model_coefs), is.null(x_model_coefs))) {
      stop("No other bias-adjusting inputs should be specified when 'x1u0_model_coefs', 'x0u1_model_coefs', and 'x1u1_model_coefs' are used.")
    }
  } else {
    stop(
      paste(
        "One of:",
        "1. data_validation",
        "2. (u_model_coefs & x_model_coefs)",
        "3. (x1u0_model_coefs, x0u1_model_coefs, x1u1_model_coefs)",
        "must be non-null.",
        sep = "\n"
      )
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
  } else if (!is.null(x_model_coefs)) {
    final <- adjust_uc_em_coef_single(
      data_observed = data_observed,
      u_model_coefs = u_model_coefs,
      x_model_coefs = x_model_coefs
    )
  } else if (!is.null(x1u0_model_coefs)) {
    final <- adjust_uc_em_coef_multinom(
      data_observed = data_observed,
      x1u0_model_coefs = x1u0_model_coefs,
      x0u1_model_coefs = x0u1_model_coefs,
      x1u1_model_coefs = x1u1_model_coefs
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
