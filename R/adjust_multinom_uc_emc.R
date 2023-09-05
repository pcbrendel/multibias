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
#' @param data The data set.
#' @param exposure The name of the exposure variable in the data.
#' @param outcome The name of the outcome variable in the data.
#' @param confounders The variable name(s) of the confounder(s) in the data.
#' A maximum of three confounders are allowed.
#' @param px1_u0_parameters The regression coefficients corresponding to the
#'  model: \ifelse{html}{\out{log(P(X=1,U=0)/P(X=0,U=0)) =
#'  &gamma;<sub>1,0</sub> + &gamma;<sub>1,1</sub>X* +
#'  &gamma;<sub>1,2</sub>Y + &gamma;<sub>1,2+j</sub>C<sub>j</sub>, }
#'  where X is the true (binary) exposure, U is the (binary) unmeasured
#'  confounder, X* is the (binary) misclassified exposure, Y is the (binary)
#'  outcome, C represents the vector of (binary) measured confounders (if any),
#'  and j corresponds to the number of measured
#'  confounders.}{\eqn{log(P(X=1,U=0)/P(X=0,U=0)) =}}
#' @param px0_u1_parameters The regression coefficients corresponding to the
#'  model: \ifelse{html}{\out{log(P(X=0,U=1)/P(X=0,U=0)) =
#'  &gamma;<sub>2,0</sub> + &gamma;<sub>2,1</sub>X* +
#'  &gamma;<sub>2,2</sub>Y + &gamma;<sub>2,2+j</sub>C<sub>j</sub>, }
#'  where X is the true (binary) exposure, U is the (binary) unmeasured
#'  confounder, X* is the (binary) misclassified exposure, Y is the (binary)
#'  outcome, C represents the vector of (binary) measured confounders (if any),
#'  and j corresponds to the number of measured
#'  confounders.}{\eqn{log(P(X=0,U=1)/P(X=0,U=0)) =}}
#' @param px1_u1_parameters The regression coefficients corresponding to the
#'  model: \ifelse{html}{\out{log(P(X=1,U=1)/P(X=0,U=0)) =
#'  &gamma;<sub>3,0</sub> + &gamma;<sub>3,1</sub>X* +
#'  &gamma;<sub>3,2</sub>Y + &gamma;<sub>3,2+j</sub>C<sub>j</sub>, }
#'  where X is the true (binary) exposure, U is the (binary) unmeasured
#'  confounder, X* is the (binary) misclassified exposure, Y is the (binary)
#'  outcome, C represents the vector of (binary) measured confounders (if any),
#'  and j corresponds to the number of measured
#'  confounders.}{\eqn{log(P(X=1,U=1)/P(X=0,U=0)) =}}
#' @param level Numeric value from 0-1 representing the range of the confidence
#' interval. The default value is 0.95.
#'
#' @examples
#' adjust_multinom_uc_emc(
#'  df_uc_mc,
#'  exposure = "Xstar",
#'  outcome = "Y",
#'  confounders = c("C1", "C2", "C3"),
#'  px1_u0_parameters = c(-1.37, 1.64, .71, -.43, -.43, .43),
#'  px0_u1_parameters = c(-.45, -.05, .50, .12, .10, -.11),
#'  px1_u1_parameters = c(-1.45, 1.69, 1.20, -.30, -.31, .29)
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
  px1_u0_parameters,
  px0_u1_parameters,
  px1_u1_parameters,
  level = 0.95
) {

  n <- nrow(data)
  len_c <- length(confounders)
  len_p1 <- length(px1_u0_parameters)
  len_p2 <- length(px0_u1_parameters)
  len_p3 <- length(px1_u1_parameters)

  xstar <- data[, exposure]
  y <- data[, outcome]

  if (sum(xstar %in% c(0, 1)) != n) {
    stop("Exposure must be binary.")
  }
  if (sum(y %in% c(0, 1)) != n) {
    stop("Outcome must be binary.")
  }
  if (len_p1 != len_c + 3) {
    stop("Incorrect X1_U0 parameter length.")
  }
  if (len_p2 != len_c + 3) {
    stop("Incorrect X0_U1 parameter length.")
  }
  if (len_p3 != len_c + 3) {
    stop("Incorrect X1_U1 parameter length.")
  }

  x1u0_0 <- px1_u0_parameters[1]
  x1u0_xstar <- px1_u0_parameters[2]
  x1u0_y <- px1_u0_parameters[3]

  x0u1_0 <- px0_u1_parameters[1]
  x0u1_xstar <- px0_u1_parameters[2]
  x0u1_y <- px0_u1_parameters[3]

  x1u1_0 <- px1_u1_parameters[1]
  x1u1_xstar <- px1_u1_parameters[2]
  x1u1_y <- px1_u1_parameters[3]

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

    x1u0_c1 <- px1_u0_parameters[4]
    x0u1_c1 <- px0_u1_parameters[4]
    x1u1_c1 <- px1_u1_parameters[4]

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

    x1u0_c1 <- px1_u0_parameters[4]
    x1u0_c2 <- px1_u0_parameters[5]

    x0u1_c1 <- px0_u1_parameters[4]
    x0u1_c2 <- px0_u1_parameters[5]

    x1u1_c1 <- px1_u1_parameters[4]
    x1u1_c2 <- px1_u1_parameters[5]

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

    x1u0_c1 <- px1_u0_parameters[4]
    x1u0_c2 <- px1_u0_parameters[5]
    x1u0_c3 <- px1_u0_parameters[6]

    x0u1_c1 <- px0_u1_parameters[4]
    x0u1_c2 <- px0_u1_parameters[5]
    x0u1_c3 <- px0_u1_parameters[6]

    x1u1_c1 <- px1_u1_parameters[4]
    x1u1_c2 <- px1_u1_parameters[5]
    x1u1_c3 <- px1_u1_parameters[6]

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
