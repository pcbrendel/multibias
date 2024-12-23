adjust_uc_om_val <- function(
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
  if (is.null(data_validation$misclassified_outcome)) {
    stop(
      paste0(
        "This function is adjusting for a misclassified outcome.",
        "\n",
        "Validation data must have a true and misclassified outcome specified."
      )
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

  uc <- setdiff(data_validation$confounders, data_observed$confounders)
  df_val$U <- data_validation$data[, uc]
  df_val <- bind_cols(
    df_val,
    data_validation$data %>%
      select(all_of(data_observed$confounders))
  )

  force_match(
    df$X,
    df_val$X,
    "Outcomes from both datasets must both be binary or both be continuous."
  )
  force_binary(
    df_val$U,
    "Uncontrolled confounder in validation data must be a binary integer."
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

  y_mod <- glm(Y ~ X + Ystar + . - U,
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

  u_mod <- glm(U ~ X + Y,
    family = binomial(link = "logit"),
    data = df_val
  )

  u_mod_coefs <- coef(u_mod)
  u_pred <- u_mod_coefs[1]

  for (i in 2:length(u_mod_coefs)) {
    var_name <- names(u_mod_coefs)[i]
    var_name <- gsub("Y", "Ypred", var_name) # col Y is not in df
    u_pred <- u_pred + df[[var_name]] * u_mod_coefs[i]
  }

  df$Upred <- rbinom(n, 1, plogis(u_pred))

  final <- glm(
    Ypred ~ X + Upred + . - Ystar,
    family = binomial(link = "logit"),
    data = df
  )

  return(final)
}


uc_om_single <- function(
    data_observed,
    u_model_coefs,
    y_model_coefs) {
  data <- data_observed$data
  n <- nrow(data)
  confounders <- data_observed$confounders
  len_c <- length(confounders)
  len_u_coefs <- length(u_model_coefs)
  len_y_coefs <- length(y_model_coefs)

  x <- data[, data_observed$exposure]
  ystar <- data[, data_observed$outcome]

  force_binary(ystar, "Outcome must be a binary integer.")
  force_len(
    len_u_coefs,
    3,
    paste0(
      "Incorrect length of U model coefficients. ",
      "Length should equal 3."
    )
  )
  force_len(
    len_y_coefs,
    3 + len_c,
    paste0(
      "Incorrect length of Y model coefficients. ",
      "Length should equal 3 + number of confounders."
    )
  )

  u1_0 <- u_model_coefs[1]
  u1_x <- u_model_coefs[2]
  u1_y <- u_model_coefs[3]

  y1_0 <- y_model_coefs[1]
  y1_x <- y_model_coefs[2]
  y1_ystar <- y_model_coefs[3]

  if (is.null(confounders)) {
    df <- data.frame(X = x, Ystar = ystar)
    df$Ypred <- rbinom(n, 1, plogis(y1_0 + y1_x * df$X + y1_ystar * df$Ystar))
    df$Upred <- rbinom(n, 1, plogis(u1_0 + u1_x * df$X + u1_y * df$Ypred))

    final <- glm(
      Ypred ~ X + Upred,
      family = binomial(link = "logit"),
      data = df
    )
  } else if (len_c == 1) {
    c1 <- data[, confounders]
    df <- data.frame(X = x, Ystar = ystar, C1 = c1)

    y1_c1 <- y_model_coefs[4]

    df$Ypred <- rbinom(
      n, 1, plogis(
        y1_0 + y1_x * df$X + y1_ystar * df$Ystar + y1_c1 * df$C1
      )
    )
    df$Upred <- rbinom(n, 1, plogis(u1_0 + u1_x * df$X + u1_y * df$Ypred))

    final <- glm(
      Ypred ~ X + C1 + Upred,
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
    df$Upred <- rbinom(n, 1, plogis(u1_0 + u1_x * df$X + u1_y * df$Ypred))

    final <- glm(
      Ypred ~ X + C1 + C2 + Upred,
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
      n, 1,
      plogis(
        y1_0 + y1_x * df$X + y1_ystar * df$Ystar +
          y1_c1 * df$C1 + y1_c2 * df$C2 + y1_c3 * df$C3
      )
    )
    df$Upred <- rbinom(n, 1, plogis(u1_0 + u1_x * df$X + u1_y * df$Ypred))

    final <- glm(
      Ypred ~ X + C1 + C2 + C3 + Upred,
      family = binomial(link = "logit"),
      data = df
    )
  } else if (len_c > 3) {
    stop("This function is currently not compatible with >3 confounders.")
  }

  return(final)
}


uc_om_multinom <- function(
    data_observed,
    u1y0_model_coefs,
    u0y1_model_coefs,
    u1y1_model_coefs) {
  data <- data_observed$data
  n <- nrow(data)
  confounders <- data_observed$confounders
  len_c <- length(confounders)
  len_u1y0_coefs <- length(u1y0_model_coefs)
  len_u0y1_coefs <- length(u0y1_model_coefs)
  len_u1y1_coefs <- length(u1y1_model_coefs)

  x <- data[, data_observed$exposure]
  ystar <- data[, data_observed$outcome]

  force_binary(ystar, "Outcome must be a binary integer.")
  force_len(
    len_u1y0_coefs,
    3 + len_c,
    paste0(
      "Incorrect length of U1Y0 model coefficients. ",
      "Length should equal 3 + number of confounders."
    )
  )
  force_len(
    len_u0y1_coefs,
    3 + len_c,
    paste0(
      "Incorrect length of U0Y1 model coefficients. ",
      "Length should equal 3 + number of confounders."
    )
  )
  force_len(
    len_u1y1_coefs,
    3 + len_c,
    paste0(
      "Incorrect length of U1Y1 model coefficients. ",
      "Length should equal 3 + number of confounders."
    )
  )

  u1y0_0 <- u1y0_model_coefs[1]
  u1y0_x <- u1y0_model_coefs[2]
  u1y0_ystar <- u1y0_model_coefs[3]

  u0y1_0 <- u0y1_model_coefs[1]
  u0y1_x <- u0y1_model_coefs[2]
  u0y1_ystar <- u0y1_model_coefs[3]

  u1y1_0 <- u1y1_model_coefs[1]
  u1y1_x <- u1y1_model_coefs[2]
  u1y1_ystar <- u1y1_model_coefs[3]

  if (is.null(confounders)) {
    df <- data.frame(X = x, Ystar = ystar)

    p_u1y0 <- exp(u1y0_0 + u1y0_x * df$X + u1y0_ystar * df$Ystar)
    p_u0y1 <- exp(u0y1_0 + u0y1_x * df$X + u0y1_ystar * df$Ystar)
    p_u1y1 <- exp(u1y1_0 + u1y1_x * df$X + u1y1_ystar * df$Ystar)

    denom <- (1 + p_u1y0 + p_u0y1 + p_u1y1)

    u0y0_pred <- 1 / denom
    u1y0_pred <- p_u1y0 / denom
    u0y1_pred <- p_u0y1 / denom
    u1y1_pred <- p_u1y1 / denom

    df_uy_pred <- data.frame(
      U0Y0 = u0y0_pred,
      U1Y0 = u1y0_pred,
      U0Y1 = u0y1_pred,
      U1Y1 = u1y1_pred
    )
    df_uy_pred4 <- bind_rows(df_uy_pred, df_uy_pred, df_uy_pred, df_uy_pred)

    combined <- bind_rows(df, df, df, df) %>%
      bind_cols(df_uy_pred4) %>%
      mutate(
        Ubar = rep(c(1, 0, 1, 0), each = n),
        Ybar = rep(c(1, 1, 0, 0), each = n),
        pUY = case_when(
          Ubar == 0 & Ybar == 0 ~ U0Y0,
          Ubar == 1 & Ybar == 0 ~ U1Y0,
          Ubar == 0 & Ybar == 1 ~ U0Y1,
          Ubar == 1 & Ybar == 1 ~ U1Y1
        )
      )
    suppressWarnings({
      final <- glm(
        Ybar ~ X + Ubar,
        family = binomial(link = "logit"),
        weights = combined$pUY,
        data = combined
      )
    })
  } else if (len_c == 1) {
    c1 <- data[, confounders]

    df <- data.frame(X = x, Ystar = ystar, C1 = c1)

    u1y0_c1 <- u1y0_model_coefs[4]
    u0y1_c1 <- u0y1_model_coefs[4]
    u1y1_c1 <- u1y1_model_coefs[4]

    p_u1y0 <- exp(
      u1y0_0 + u1y0_x * df$X + u1y0_ystar * df$Ystar +
        u1y0_c1 * df$C1
    )
    p_u0y1 <- exp(
      u0y1_0 + u0y1_x * df$X + u0y1_ystar * df$Ystar +
        u0y1_c1 * df$C1
    )
    p_u1y1 <- exp(
      u1y1_0 + u1y1_x * df$X + u1y1_ystar * df$Ystar +
        u1y1_c1 * df$C1
    )

    denom <- (1 + p_u1y0 + p_u0y1 + p_u1y1)

    u0y0_pred <- 1 / denom
    u1y0_pred <- p_u1y0 / denom
    u0y1_pred <- p_u0y1 / denom
    u1y1_pred <- p_u1y1 / denom

    df_uy_pred <- data.frame(
      U0Y0 = u0y0_pred,
      U1Y0 = u1y0_pred,
      U0Y1 = u0y1_pred,
      U1Y1 = u1y1_pred
    )
    df_uy_pred4 <- bind_rows(df_uy_pred, df_uy_pred, df_uy_pred, df_uy_pred)

    combined <- bind_rows(df, df, df, df) %>%
      bind_cols(df_uy_pred4) %>%
      mutate(
        Ubar = rep(c(1, 0, 1, 0), each = n),
        Ybar = rep(c(1, 1, 0, 0), each = n),
        pUY = case_when(
          Ubar == 0 & Ybar == 0 ~ U0Y0,
          Ubar == 1 & Ybar == 0 ~ U1Y0,
          Ubar == 0 & Ybar == 1 ~ U0Y1,
          Ubar == 1 & Ybar == 1 ~ U1Y1
        )
      )

    suppressWarnings({
      final <- glm(
        Ybar ~ X + C1 + Ubar,
        family = binomial(link = "logit"),
        weights = combined$pUY,
        data = combined
      )
    })
  } else if (len_c == 2) {
    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]

    df <- data.frame(X = x, Ystar = ystar, C1 = c1, C2 = c2)

    u1y0_c1 <- u1y0_model_coefs[4]
    u1y0_c2 <- u1y0_model_coefs[5]

    u0y1_c1 <- u0y1_model_coefs[4]
    u0y1_c2 <- u0y1_model_coefs[5]

    u1y1_c1 <- u1y1_model_coefs[4]
    u1y1_c2 <- u1y1_model_coefs[5]

    p_u1y0 <- exp(
      u1y0_0 + u1y0_x * df$X + u1y0_ystar * df$Ystar +
        u1y0_c1 * df$C1 + u1y0_c2 * df$C2
    )
    p_u0y1 <- exp(
      u0y1_0 + u0y1_x * df$X + u0y1_ystar * df$Ystar +
        u0y1_c1 * df$C1 + u0y1_c2 * df$C2
    )
    p_u1y1 <- exp(
      u1y1_0 + u1y1_x * df$X + u1y1_ystar * df$Ystar +
        u1y1_c1 * df$C1 + u1y1_c2 * df$C2
    )

    denom <- (1 + p_u1y0 + p_u0y1 + p_u1y1)

    u0y0_pred <- 1 / denom
    u1y0_pred <- p_u1y0 / denom
    u0y1_pred <- p_u0y1 / denom
    u1y1_pred <- p_u1y1 / denom

    df_uy_pred <- data.frame(
      U0Y0 = u0y0_pred,
      U1Y0 = u1y0_pred,
      U0Y1 = u0y1_pred,
      U1Y1 = u1y1_pred
    )
    df_uy_pred4 <- bind_rows(df_uy_pred, df_uy_pred, df_uy_pred, df_uy_pred)

    combined <- bind_rows(df, df, df, df) %>%
      bind_cols(df_uy_pred4) %>%
      mutate(
        Ubar = rep(c(1, 0, 1, 0), each = n),
        Ybar = rep(c(1, 1, 0, 0), each = n),
        pUY = case_when(
          Ubar == 0 & Ybar == 0 ~ U0Y0,
          Ubar == 1 & Ybar == 0 ~ U1Y0,
          Ubar == 0 & Ybar == 1 ~ U0Y1,
          Ubar == 1 & Ybar == 1 ~ U1Y1
        )
      )

    suppressWarnings({
      final <- glm(
        Ybar ~ X + C1 + C2 + Ubar,
        family = binomial(link = "logit"),
        weights = combined$pUY,
        data = combined
      )
    })
  } else if (len_c == 3) {
    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]
    c3 <- data[, confounders[3]]

    df <- data.frame(X = x, Ystar = ystar, C1 = c1, C2 = c2, C3 = c3)

    u1y0_c1 <- u1y0_model_coefs[4]
    u1y0_c2 <- u1y0_model_coefs[5]
    u1y0_c3 <- u1y0_model_coefs[6]

    u0y1_c1 <- u0y1_model_coefs[4]
    u0y1_c2 <- u0y1_model_coefs[5]
    u0y1_c3 <- u0y1_model_coefs[6]

    u1y1_c1 <- u1y1_model_coefs[4]
    u1y1_c2 <- u1y1_model_coefs[5]
    u1y1_c3 <- u1y1_model_coefs[6]

    p_u1y0 <- exp(
      u1y0_0 + u1y0_x * df$X + u1y0_ystar * df$Ystar +
        u1y0_c1 * df$C1 + u1y0_c2 * df$C2 + u1y0_c3 * df$C3
    )
    p_u0y1 <- exp(
      u0y1_0 + u0y1_x * df$X + u0y1_ystar * df$Ystar +
        u0y1_c1 * df$C1 + u0y1_c2 * df$C2 + u0y1_c3 * df$C3
    )
    p_u1y1 <- exp(
      u1y1_0 + u1y1_x * df$X + u1y1_ystar * df$Ystar +
        u1y1_c1 * df$C1 + u1y1_c2 * df$C2 + u1y1_c3 * df$C3
    )

    denom <- (1 + p_u1y0 + p_u0y1 + p_u1y1)

    u0y0_pred <- 1 / denom
    u1y0_pred <- p_u1y0 / denom
    u0y1_pred <- p_u0y1 / denom
    u1y1_pred <- p_u1y1 / denom

    df_uy_pred <- data.frame(
      U0Y0 = u0y0_pred,
      U1Y0 = u1y0_pred,
      U0Y1 = u0y1_pred,
      U1Y1 = u1y1_pred
    )
    df_uy_pred4 <- bind_rows(df_uy_pred, df_uy_pred, df_uy_pred, df_uy_pred)

    combined <- bind_rows(df, df, df, df) %>%
      bind_cols(df_uy_pred4) %>%
      mutate(
        Ubar = rep(c(1, 0, 1, 0), each = n),
        Ybar = rep(c(1, 1, 0, 0), each = n),
        pUY = case_when(
          Ubar == 0 & Ybar == 0 ~ U0Y0,
          Ubar == 1 & Ybar == 0 ~ U1Y0,
          Ubar == 0 & Ybar == 1 ~ U0Y1,
          Ubar == 1 & Ybar == 1 ~ U1Y1
        )
      )

    suppressWarnings({
      final <- glm(
        Ybar ~ X + C1 + C2 + C3 + Ubar,
        family = binomial(link = "logit"),
        weights = combined$pUY,
        data = combined
      )
    })
  } else if (len_c > 3) {
    stop("This function is currently not compatible with >3 confounders.")
  }

  return(final)
}


#' Adust for uncontrolled confounding and outcome misclassification.
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `adjust_uc_omc_sel()` was renamed to `adjust_uc_om_sel()`
#' @keywords internal
#'
#' @export
adjust_uc_omc <- function(
    data_observed,
    u_model_coefs = NULL,
    y_model_coefs = NULL,
    u0y1_model_coefs = NULL,
    u1y0_model_coefs = NULL,
    u1y1_model_coefs = NULL,
    level = 0.95) {
  lifecycle::deprecate_warn(
    "1.5.3", "adjust_uc_omc()", "adjust_uc_om()"
  )
  adjust_uc_om(
    data_observed,
    u_model_coefs,
    y_model_coefs,
    u0y1_model_coefs,
    u1y0_model_coefs,
    u1y1_model_coefs,
    level
  )
}


#' Adust for uncontrolled confounding and outcome misclassification.
#'
#' `adjust_uc_om` returns the exposure-outcome odds ratio and confidence
#' interval, adjusted for uncontrolled confounding and outcome
#' misclassificaiton. Two different options for the bias parameters are
#' available here: 1) parameters from separate models of *U* and *Y*
#' (`u_model_coefs` and `y_model_coefs`) or 2) parameters from
#' a joint model of *U* and *Y* (`u1y0_model_coefs`,
#' `u0y1_model_coefs`, and `u1y1_model_coefs`).
#'
#' Bias adjustment can be performed by inputting either a validation dataset or
#' the necessary bias parameters.
#' Values for the bias parameters can be applied as
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
#' data, plus data for the true and misclassified outcome corresponding to the
#' observed exposure in `data_observed`.
#' There should also be data for the confounder missing in `data_observed`.
#' @param u_model_coefs The regression coefficients corresponding to the model:
#' \ifelse{html}{\out{logit(P(U=1)) = &alpha;<sub>0</sub> + &alpha;<sub>1</sub>X + &alpha;<sub>2</sub>Y, }}{\eqn{logit(P(U=1)) = \alpha_0 + \alpha_1 X + \alpha_2 Y, }}
#' where *U* is the binary unmeasured confounder, *X* is the
#' exposure, *Y* is the binary true outcome. The number of parameters
#' therefore equals 3.
#' @param y_model_coefs The regression coefficients corresponding to the model:
#' \ifelse{html}{\out{logit(P(Y=1)) = &delta;<sub>0</sub> + &delta;<sub>1</sub>X + &delta;<sub>2</sub>Y* + &delta;<sub>2+j</sub>C<sub>j</sub>, }}{\eqn{logit(P(Y=1)) = \delta_0 + \delta_1 X + \delta_2 Y^* + \delta_{2+j} C_j, }}
#' where *Y* represents binary true outcome, *X* is the exposure,
#' *Y** is the binary misclassified outcome,
#' *C* represents the vector of measured confounders (if any),
#' and *j* corresponds to the number of measured confounders.
#' The number of parameters therefore equals 3 + *j*.
#' @param u1y0_model_coefs The regression coefficients corresponding to the
#' model:
#' \ifelse{html}{\out{log(P(U=1,Y=0)/P(U=0,Y=0)) = &gamma;<sub>1,0</sub> + &gamma;<sub>1,1</sub>X + &gamma;<sub>1,2</sub>Y* + &gamma;<sub>1,2+j</sub>C<sub>j</sub>, }}{\eqn{log(P(U=1,Y=0)/P(U=0,Y=0)) = \gamma_{1,0} + \gamma_{1,1} X + \gamma_{1,2} Y^* + \gamma_{1,2+j} C_j, }}
#' where *U* is the binary unmeasured confounder, *Y* is the
#' binary true outcome, *X* is the exposure, *Y** is the binary
#' misclassified outcome, *C* represents the vector of measured
#' confounders (if any), and *j* corresponds to the number of
#' measured confounders.
#' @param u0y1_model_coefs The regression coefficients corresponding to the
#' model:
#' \ifelse{html}{\out{log(P(U=0,Y=1)/P(U=0,Y=0)) = &gamma;<sub>2,0</sub> + &gamma;<sub>2,1</sub>X + &gamma;<sub>2,2</sub>Y* + &gamma;<sub>2,2+j</sub>C<sub>j</sub>, }}{\eqn{log(P(U=0,Y=1)/P(U=0,Y=0)) = \gamma_{2,0} + \gamma_{2,1} X + \gamma_{2,2} Y^* + \gamma_{2,2+j} C_j,}}
#' where *U* is the binary unmeasured confounder, *Y* is the
#' binary true outcome, *X* is the exposure, *Y** is the binary
#' misclassified outcome, *C* represents the vector of measured
#' confounders (if any), and *j* corresponds to the number of
#' measured confounders.
#' @param u1y1_model_coefs The regression coefficients corresponding to the
#' model:
#' \ifelse{html}{\out{log(P(U=1,Y=1)/P(U=0,Y=0)) = &gamma;<sub>3,0</sub> + &gamma;<sub>3,1</sub>X + &gamma;<sub>3,2</sub>Y* + &gamma;<sub>3,2+j</sub>C<sub>j</sub>, }}{\eqn{log(P(U=1,Y=1)/P(U=0,Y=0)) = \gamma_{3,0} + \gamma_{3,1} X + \gamma_{3,2} Y^* + \gamma_{3,2+j} C_j,}}
#' where *U* is the binary unmeasured confounder, *Y* is the
#' binary true outcome, *X* is the exposure, *Y** is the binary
#' misclassified outcome, *C* represents the vector of measured
#' confounders (if any), and *j* corresponds to the number of
#' measured confounders.
#' @param level Value from 0-1 representing the full range of the confidence
#' interval. Default is 0.95.
#'
#' @return A list where the first item is the odds ratio estimate of the
#' effect of the exposure on the outcome and the second item is the
#' confidence interval as the vector: (lower bound, upper bound).
#'
#' @examples
#' df_observed <- data_observed(
#'   data = df_uc_om,
#'   exposure = "X",
#'   outcome = "Ystar",
#'   confounders = "C1"
#' )
#'
#' # Using validation data -----------------------------------------------------
#' df_validation <- data_validation(
#'   data = df_uc_om_source,
#'   true_exposure = "X",
#'   true_outcome = "Y",
#'   confounders = c("C1", "U"),
#'   misclassified_outcome = "Ystar"
#' )
#'
#' adjust_uc_om(
#'   data_observed = df_observed,
#'   data_validation = df_validation
#' )
#'
#' # Using u_model_coefs and y_model_coefs -------------------------------------
#' adjust_uc_om(
#'   data_observed = df_observed,
#'   u_model_coefs = c(-0.22, 0.61, 0.70),
#'   y_model_coefs = c(-2.85, 0.73, 1.60, 0.38)
#' )
#'
#' # Using u1y0_model_coefs, u0y1_model_coefs, u1y1_model_coefs ----------------
#' adjust_uc_om(
#'   data_observed = df_observed,
#'   u1y0_model_coefs = c(-0.19, 0.61, 0.00, -0.07),
#'   u0y1_model_coefs = c(-3.21, 0.60, 1.60, 0.36),
#'   u1y1_model_coefs = c(-2.72, 1.24, 1.59, 0.34)
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

adjust_uc_om <- function(
    data_observed,
    data_validation = NULL,
    u_model_coefs = NULL,
    y_model_coefs = NULL,
    u1y0_model_coefs = NULL,
    u0y1_model_coefs = NULL,
    u1y1_model_coefs = NULL,
    level = 0.95) {
  if (!is.null(data_validation)) {
    if (!all(is.null(u_model_coefs), is.null(y_model_coefs), is.null(u1y0_model_coefs), is.null(u0y1_model_coefs), is.null(u1y1_model_coefs))) {
      stop("No bias parameters should be specified when 'data_validation' is used.")
    }
  } else if (!is.null(u_model_coefs) && !is.null(y_model_coefs)) {
    if (!all(is.null(data_validation), is.null(u1y0_model_coefs), is.null(u0y1_model_coefs), is.null(u1y1_model_coefs))) {
      stop("No other bias-adjusting inputs should be specified when 'u_model_coefs' and 'y_model_coefs' are used.")
    }
  } else if (!is.null(u1y0_model_coefs) && !is.null(u0y1_model_coefs) && !is.null(u1y1_model_coefs)) {
    if (!all(is.null(data_validation), is.null(u_model_coefs), is.null(y_model_coefs))) {
      stop("No other bias-adjusting inputs should be specified when 'u1y0_model_coefs', 'u0y1_model_coefs', and 'u1y1_model_coefs' are used.")
    }
  } else {
    stop(
      paste(
        "One of:",
        "1. data_validation",
        "2. (u_model_coefs & y_model_coefs)",
        "3. (u1y0_model_coefs, u0y1_model_coefs, u1y1_model_coefs)",
        "must be non-null.",
        sep = "\n"
      )
    )
  }

  data <- data_observed$data

  if (!is.null(data_validation)) {
    final <- adjust_uc_om_val(
      data_observed,
      data_validation
    )
  } else if (!is.null(u_model_coefs)) {
    final <- uc_om_single(
      data_observed = data_observed,
      u_model_coefs = u_model_coefs,
      y_model_coefs = y_model_coefs
    )
  } else if (!is.null(u1y0_model_coefs)) {
    final <- uc_om_multinom(
      data_observed = data_observed,
      u1y0_model_coefs = u1y0_model_coefs,
      u0y1_model_coefs = u0y1_model_coefs,
      u1y1_model_coefs = u1y1_model_coefs
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
