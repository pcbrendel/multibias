#' Adust for uncontrolled confounding, exposure misclassification, and selection
#' bias.
#'
#' \code{adjust_multinom_uc_emc_sel} returns the exposure-outcome odds ratio and
#' confidence interval, adjusted for uncontrolled confounding, exposure
#' misclassificaiton, and selection bias.
#'
#' This function uses one bias model, a multinomial logistic regression model,
#' to predict the uncontrolled confounder (U) and exposure (X). If separate bias
#' models for X and U are desired, use \code{adjust_uc_mc_2}.
#'
#' @param data The data set.
#' @param exposure The name of the exposure variable in the data.
#' @param outcome The name of the outcome variable in the data.
#' @param confounders The variable name(s) of the confounder(s) in the data.
#' A maximum of three confounders are allowed.
#' @param px1_u0_parameters The regression coefficients corresponding to the
#'  model: \ifelse{html}{\out{log(P(X=1,U=0)/P(X=0,U=0)) =
#'  &gamma;<sub>1,0</sub> + &gamma;<sub>1,1</sub>X* +
#'  &gamma;<sub>1,2</sub>Y + &gamma;<sub>1,2+j</sub>C<sub>j</sub>, }where X
#'  is the true (binary) exposure, U is the (binary) unmeasured confounder,
#'  X* is the (binary) misclassified exposure, Y is the (binary) outcome, C
#'  represents the vector of (binary) measured confounders (if any), and
#'  j corresponds to the number of measured
#'  confounders.}{\eqn{log(P(X=1,U=0)/P(X=0,U=0)) =}}
#' @param px0_u1_parameters The regression coefficients corresponding to the
#'  model: \ifelse{html}{\out{log(P(X=0,U=1)/P(X=0,U=0)) =
#'  &gamma;<sub>2,0</sub> + &gamma;<sub>2,1</sub>X* +
#'  &gamma;<sub>2,2</sub>Y + &gamma;<sub>2,2+j</sub>C<sub>j</sub>, } where X
#'  is the true (binary) exposure, U is the (binary) unmeasured confounder,
#'  X* is the (binary) misclassified exposure, Y is the (binary) outcome,
#'  C represents the vector of (binary) measured confounders (if any), and
#'  j corresponds to the number of measured
#'  confounders.}{\eqn{log(P(X=0,U=1)/P(X=0,U=0)) =}}
#' @param px1_u1_parameters The regression coefficients corresponding to the
#'  model: \ifelse{html}{\out{log(P(X=1,U=1)/P(X=0,U=0)) =
#'  &gamma;<sub>3,0</sub> + &gamma;<sub>3,1</sub>X* +
#'  &gamma;<sub>3,2</sub>Y + &gamma;<sub>3,2+j</sub>C<sub>j</sub>, } where X
#'  is the true (binary) exposure, U is the (binary) unmeasured confounder,
#'  X* is the (binary) misclassified exposure, Y is the (binary) outcome,
#'  C represents the vector of (binary) measured confounders (if any),
#'  and j corresponds to the number of measured
#'  confounders.}{\eqn{log(P(X=1,U=1)/P(X=0,U=0)) =}}
#' @param ps1_parameters The regression coefficients corresponding to the model:
#'  \ifelse{html}{\out{logit(P(S=1)) = &beta;<sub>0</sub> +
#'  &beta;<sub>1</sub>X* + &beta;<sub>2</sub>Y +
#'  &beta;<sub>2+j</sub>C<sub>j</sub>, } where S represents (binary) selection,
#'  X* is the (binary) misclassified exposure, Y is the (binary) outcome,
#'  C represents the vector of (binary) measured confounders (if any), and
#'  j corresponds to the number of measured confounders.}{\eqn{logit(P(S=1)) =}}
#' @param level Number from 0-1 representing the range of the confidence
#'  interval. Default is 0.95.
#'
#' @examples
#' adjust_multinom_uc_emc_sel(
#'  df_uc_mc_sel,
#'  exposure = "Xstar",
#'  outcome = "Y",
#'  confounders = c("C1", "C2", "C3"),
#'  px1_u0_parameters = c(-1.79, 2.71, .57, -.43, -.43, .42),
#'  px0_u1_parameters = c(-.47, 0, .50, .12, .10, -.11),
#'  px1_u1_parameters = c(-1.84, 2.70, 1.07, -.30, -.31, .29),
#'  ps1_parameters = c(-.39, .40, .75, -.04, -.04, .05)
#' )
#'
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom stats binomial
#' @importFrom stats glm
#' @importFrom stats qnorm
#' @importFrom stats plogis
#' @importFrom rlang .data
#'
#' @export

adjust_multinom_uc_emc_sel <- function(
  data,
  exposure,
  outcome,
  confounders = NULL,
  px1_u0_parameters,
  px0_u1_parameters,
  px1_u1_parameters,
  ps1_parameters,
  level = 0.95
) {

  n <- nrow(data)
  len_c <- length(confounders)
  len_p1 <- length(px1_u0_parameters)
  len_p2 <- length(px0_u1_parameters)
  len_p3 <- length(px1_u1_parameters)
  len_p4 <- length(ps1_parameters)

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
  if (len_p4 != len_c + 3) {
    stop("Incorrect S1 parameter length.")
  }

  s1_0     <- ps1_parameters[1]
  s1_xstar <- ps1_parameters[2]
  s1_y     <- ps1_parameters[3]

  x1u0_0     <- px1_u0_parameters[1]
  x1u0_xstar <- px1_u0_parameters[2]
  x1u0_y     <- px1_u0_parameters[3]

  x0u1_0     <- px0_u1_parameters[1]
  x0u1_xstar <- px0_u1_parameters[2]
  x0u1_y     <- px0_u1_parameters[3]

  x1u1_0     <- px1_u1_parameters[1]
  x1u1_xstar <- px1_u1_parameters[2]
  x1u1_y     <- px1_u1_parameters[3]

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
                             Xbar == 1 & Ubar == 1 ~ X1U1),
             pS = plogis(s1_0 + s1_xstar * .data$Xstar + s1_y * .data$Y))

    final <- glm(
      Y ~ Xbar + Ubar,
      family = binomial(link = "logit"),
      weights = (combined$pXU / combined$pS),
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

  }else if (len_c == 1) {

    c1 <- data[, confounders]
    df <- data.frame(Xstar = xstar, Y = y, C1 = c1)

    s1_c1   <- ps1_parameters[4]
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
                             Xbar == 1 & Ubar == 1 ~ X1U1),
             pS = plogis(s1_0 + s1_xstar * .data$Xstar + s1_y * .data$Y +
                           s1_c1 * .data$C1))

    final <- glm(
      Y ~ Xbar + C1 + Ubar,
      family = binomial(link = "logit"),
      weights = (combined$pXU / combined$pS),
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

    s1_c1 <- ps1_parameters[4]
    s1_c2 <- ps1_parameters[5]

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
                             Xbar == 1 & Ubar == 1 ~ X1U1),
             pS = plogis(s1_0 + s1_xstar * .data$Xstar + s1_y * .data$Y +
                           s1_c1 * .data$C1 + s1_c2 * .data$C2))

    final <- glm(
      Y ~ Xbar + C1 + C2 + Ubar,
      family = binomial(link = "logit"),
      weights = (combined$pXU / combined$pS),
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

    s1_c1 <- ps1_parameters[4]
    s1_c2 <- ps1_parameters[5]
    s1_c3 <- ps1_parameters[6]

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
                             Xbar == 1 & Ubar == 1 ~ X1U1),
             pS = plogis(s1_0 + s1_xstar * .data$Xstar + s1_y * .data$Y +
                           s1_c1 * .data$C1 + s1_c2 * .data$C2 +
                           s1_c3 * .data$C3))

    final <- glm(
      Y ~ Xbar + C1 + C2 + C3 + Ubar,
      family = binomial(link = "logit"),
      weights = (combined$pXU / combined$pS),
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

  } else if (len_c > 3) {
    stop("This function is currently not compatible with >3 confounders.")
  }
}