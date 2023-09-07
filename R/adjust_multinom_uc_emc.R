#' Adust for uncontrolled confounding and exposure misclassification.
#'
#' \code{adjust_multinom_uc_emc} returns the exposure-outcome odds ratio and
#' confidence interval, adjusted for uncontrolled confounding and exposure
#' misclassificaiton.
#'
#' This function uses one bias model, a multinomial logistic regression model,
#' to predict the uncontrolled confounder (U) and exposure (X). If separate bias
#' models for X and U are desired, use \code{adjust_uc_emc_2}.
#'
#' @inheritParams adjust_emc_sel
#' @param x1u0_model_coefs The regression coefficients corresponding to the
#'  model: \ifelse{html}{\out{log(P(X=1,U=0)/P(X=0,U=0)) =
#'  &gamma;<sub>1,0</sub> + &gamma;<sub>1,1</sub>X* +
#'  &gamma;<sub>1,2</sub>Y + &gamma;<sub>1,2+j</sub>C<sub>j</sub>, }
#'  where X is the true (binary) exposure, U is the (binary) unmeasured
#'  confounder, X* is the (binary) misclassified exposure, Y is the (binary)
#'  outcome, C represents the vector of (binary) measured confounders (if any),
#'  and j corresponds to the number of measured
#'  confounders.}{\eqn{log(P(X=1,U=0)/P(X=0,U=0)) =}}
#' @param x0u1_model_coefs The regression coefficients corresponding to the
#'  model: \ifelse{html}{\out{log(P(X=0,U=1)/P(X=0,U=0)) =
#'  &gamma;<sub>2,0</sub> + &gamma;<sub>2,1</sub>X* +
#'  &gamma;<sub>2,2</sub>Y + &gamma;<sub>2,2+j</sub>C<sub>j</sub>, }
#'  where X is the true (binary) exposure, U is the (binary) unmeasured
#'  confounder, X* is the (binary) misclassified exposure, Y is the (binary)
#'  outcome, C represents the vector of (binary) measured confounders (if any),
#'  and j corresponds to the number of measured
#'  confounders.}{\eqn{log(P(X=0,U=1)/P(X=0,U=0)) =}}
#' @param x1u1_model_coefs The regression coefficients corresponding to the
#'  model: \ifelse{html}{\out{log(P(X=1,U=1)/P(X=0,U=0)) =
#'  &gamma;<sub>3,0</sub> + &gamma;<sub>3,1</sub>X* +
#'  &gamma;<sub>3,2</sub>Y + &gamma;<sub>3,2+j</sub>C<sub>j</sub>, }
#'  where X is the true (binary) exposure, U is the (binary) unmeasured
#'  confounder, X* is the (binary) misclassified exposure, Y is the (binary)
#'  outcome, C represents the vector of (binary) measured confounders (if any),
#'  and j corresponds to the number of measured
#'  confounders.}{\eqn{log(P(X=1,U=1)/P(X=0,U=0)) =}}
#'
#' @examples
#' adjust_multinom_uc_emc(
#'  df_uc_emc,
#'  exposure = "Xstar",
#'  outcome = "Y",
#'  confounders = "C1",
#'  x1u0_model_coefs = c(-2.82, 1.62, 0.68, -0.06),
#'  x0u1_model_coefs = c(-0.20, 0.00, 0.68, -0.05),
#'  x1u1_model_coefs = c(-2.36, 1.62, 1.29, 0.27)
#' )
#'
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom stats binomial
#' @importFrom stats glm
#' @importFrom stats qnorm
#'
#' @export

adjust_multinom_uc_emc <- function(
  data,
  exposure,
  outcome,
  confounders = NULL,
  x1u0_model_coefs,
  x0u1_model_coefs,
  x1u1_model_coefs,
  level = 0.95
) {

  n <- nrow(data)
  len_c <- length(confounders)
  len_x1u0_coefs <- length(x1u0_model_coefs)
  len_x0u1_coefs <- length(x0u1_model_coefs)
  len_x1u1_coefs <- length(x1u1_model_coefs)

  xstar <- data[, exposure]
  y <- data[, outcome]

  if (sum(xstar %in% c(0, 1)) != n) {
    stop("Exposure must be a binary integer.")
  }
  if (sum(y %in% c(0, 1)) != n) {
    stop("Outcome must be a binary integer.")
  }
  if (len_x1u0_coefs != 3 + len_c) {
    stop(
      paste0(
        "Incorrect length of X1U0 model coefficients. ",
        "Length should equal 3 + number of confounders."
      )
    )
  }
  if (len_x0u1_coefs != 3 + len_c) {
    stop(
      paste0(
        "Incorrect length of X0U1 model coefficients. ",
        "Length should equal 3 + number of confounders."
      )
    )
  }
  if (len_x1u1_coefs != 3 + len_c) {
    stop(
      paste0(
        "Incorrect length of X1U1 model coefficients. ",
        "Length should equal 3 + number of confounders."
      )
    )
  }

  x1u0_0     <- x1u0_model_coefs[1]
  x1u0_xstar <- x1u0_model_coefs[2]
  x1u0_y     <- x1u0_model_coefs[3]

  x0u1_0     <- x0u1_model_coefs[1]
  x0u1_xstar <- x0u1_model_coefs[2]
  x0u1_y     <- x0u1_model_coefs[3]

  x1u1_0     <- x1u1_model_coefs[1]
  x1u1_xstar <- x1u1_model_coefs[2]
  x1u1_y     <- x1u1_model_coefs[3]

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
      mutate(Xbar = rep(c(1, 0, 1, 0), each = n),
             Ubar = rep(c(1, 1, 0, 0), each = n),
             pXU = case_when(Xbar == 0 & Ubar == 0 ~ X0U0,
                             Xbar == 1 & Ubar == 0 ~ X1U0,
                             Xbar == 0 & Ubar == 1 ~ X0U1,
                             Xbar == 1 & Ubar == 1 ~ X1U1))

    final <- glm(
      Y ~ Xbar + Ubar,
      family = binomial(link = "logit"),
      weights = combined$pXU,
      data = combined
    )

    est <- summary(final)$coef[2, 1]
    se <- summary(final)$coef[2, 2]
    alpha <- 1 - level

    return(
      list(
        exp(est),
        c(
          exp(est + se * qnorm(alpha / 2)),
          exp(est + se * qnorm(1 - alpha / 2)))
      )
    )

  } else if (len_c == 1) {

    c1 <- data[, confounders]

    df <- data.frame(Xstar = xstar, Y = y, C1 = c1)

    x1u0_c1 <- x1u0_model_coefs[4]
    x0u1_c1 <- x0u1_model_coefs[4]
    x1u1_c1 <- x1u1_model_coefs[4]

    p_x1u0 <- exp(x1u0_0 + x1u0_xstar * df$Xstar + x1u0_y * df$Y +
                    x1u0_c1 * df$C1)
    p_x0u1 <- exp(x0u1_0 + x0u1_xstar * df$Xstar + x0u1_y * df$Y +
                    x0u1_c1 * df$C1)
    p_x1u1 <- exp(x1u1_0 + x1u1_xstar * df$Xstar + x1u1_y * df$Y +
                    x1u1_c1 * df$C1)

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
      mutate(Xbar = rep(c(1, 0, 1, 0), each = n),
             Ubar = rep(c(1, 1, 0, 0), each = n),
             pXU = case_when(Xbar == 0 & Ubar == 0 ~ X0U0,
                             Xbar == 1 & Ubar == 0 ~ X1U0,
                             Xbar == 0 & Ubar == 1 ~ X0U1,
                             Xbar == 1 & Ubar == 1 ~ X1U1))

    final <- glm(
      Y ~ Xbar + C1 + Ubar,
      family = binomial(link = "logit"),
      weights = combined$pXU,
      data = combined
    )

    est <- summary(final)$coef[2, 1]
    se <- summary(final)$coef[2, 2]
    alpha <- 1 - level

    return(
      list(
        exp(est),
        c(
          exp(est + se * qnorm(alpha / 2)),
          exp(est + se * qnorm(1 - alpha / 2))
        )
      )
    )

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

    p_x1u0 <- exp(x1u0_0 + x1u0_xstar * df$Xstar + x1u0_y * df$Y +
                    x1u0_c1 * df$C1 + x1u0_c2 * df$C2)
    p_x0u1 <- exp(x0u1_0 + x0u1_xstar * df$Xstar + x0u1_y * df$Y +
                    x0u1_c1 * df$C1 + x0u1_c2 * df$C2)
    p_x1u1 <- exp(x1u1_0 + x1u1_xstar * df$Xstar + x1u1_y * df$Y +
                    x1u1_c1 * df$C1 + x1u1_c2 * df$C2)

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
      mutate(Xbar = rep(c(1, 0, 1, 0), each = n),
             Ubar = rep(c(1, 1, 0, 0), each = n),
             pXU = case_when(Xbar == 0 & Ubar == 0 ~ X0U0,
                             Xbar == 1 & Ubar == 0 ~ X1U0,
                             Xbar == 0 & Ubar == 1 ~ X0U1,
                             Xbar == 1 & Ubar == 1 ~ X1U1))

    final <- glm(
      Y ~ Xbar + C1 + C2 + Ubar,
      family = binomial(link = "logit"),
      weights = combined$pXU,
      data = combined
    )

    est <- summary(final)$coef[2, 1]
    se <- summary(final)$coef[2, 2]
    alpha <- 1 - level

    return(
      list(
        exp(est),
        c(
          exp(est + se * qnorm(alpha / 2)),
          exp(est + se * qnorm(1 - alpha / 2))
        )
      )
    )

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

    p_x1u0 <- exp(x1u0_0 + x1u0_xstar * df$Xstar + x1u0_y * df$Y +
                    x1u0_c1 * df$C1 + x1u0_c2 * df$C2 + x1u0_c3 * df$C3)
    p_x0u1 <- exp(x0u1_0 + x0u1_xstar * df$Xstar + x0u1_y * df$Y +
                    x0u1_c1 * df$C1 + x0u1_c2 * df$C2 + x0u1_c3 * df$C3)
    p_x1u1 <- exp(x1u1_0 + x1u1_xstar * df$Xstar + x1u1_y * df$Y +
                    x1u1_c1 * df$C1 + x1u1_c2 * df$C2 + x1u1_c3 * df$C3)

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
      mutate(Xbar = rep(c(1, 0, 1, 0), each = n),
             Ubar = rep(c(1, 1, 0, 0), each = n),
             pXU = case_when(Xbar == 0 & Ubar == 0 ~ X0U0,
                             Xbar == 1 & Ubar == 0 ~ X1U0,
                             Xbar == 0 & Ubar == 1 ~ X0U1,
                             Xbar == 1 & Ubar == 1 ~ X1U1))

    final <- glm(
      Y ~ Xbar + C1 + C2 + C3 + Ubar,
      family = binomial(link = "logit"),
      weights = combined$pXU,
      data = combined
    )

    est <- summary(final)$coef[2, 1]
    se <- summary(final)$coef[2, 2]
    alpha <- 1 - level

    return(
      list(
        exp(est),
        c(
          exp(est + se * qnorm(alpha / 2)),
          exp(est + se * qnorm(1 - alpha / 2))
        )
      )
    )

  } else if (len_c > 3) {
    stop("This function is currently not compatible with >3 confounders.")
  }
}
