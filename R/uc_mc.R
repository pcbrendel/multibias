#' Adust for uncontrolled confounding and exposure misclassification.
#'
#' \code{adjust_uc_mc} returns the exposure-outcome odds ratio and confidence
#' interval, adjusted for uncontrolled confounding and exposure
#' misclassificaiton.
#'
#' Details
#'
#' @param data The data set.
#' @param exposure The variable corresponding to the exposure in the data.
#' @param outcome The variable corresponding to the outcome in the data.
#' @param confounders The variable(s) corresponding to the confounder(s) in the data.
#' A maximum of three confounders are allowed.
#' @param px1_u0_parameters The regression coefficients corresponding to the model: .
#' @param px0_u1_parameters The regression coefficients corresponding to the model: .
#' @param px1_u1_parameters The regression coefficients corresponding to the model: .
#'
#' @author Paul Brendel
#'
#' @import dplyr
#' @importFrom magrittr %>%
#'
#' @return The exposure-outcome odds ratio, adjusted for uncontrolled
#'   confounding and exposure misclassification, and its confidence interval.
#'
#' @export

adjust_uc_mc <- function (data, exposure, outcome, confounders = NULL, px1_u0_parameters,
                          px0_u1_parameters, px1_u1_parameters, level = .95) {

  n <- nrow(data)
  c <- length(confounders)
  p1 <- length(px1_u0_parameters)
  p2 <- length(px0_u1_parameters)
  p3 <- length(px1_u1_parameters)

  Xstar <- data[,exposure]
  Y <- data[,outcome]

  if (sum(Xstar %in% c(0, 1)) != n) {stop('Exposure must be binary')}
  if (sum(Y %in% c(0, 1)) != n) {stop('Outcome must be binary')}
  if (p1 != c + 3) {stop('Incorrect X1_U0 parameter length')}
  if (p2 != c + 3) {stop('Incorrect X0_U1 parameter length')}
  if (p3 != c + 3) {stop('Incorrect X1_U1 parameter length')}

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

    df <- data.frame(Xstar, Y)

    A <- exp(x1u0_0 + x1u0_xstar * df$Xstar + x1u0_y * df$Y)
    B <- exp(x0u1_0 + x0u1_xstar * df$Xstar + x0u1_y * df$Y)
    C <- exp(x1u1_0 + x1u1_xstar * df$Xstar + x1u1_y * df$Y)

    denom <- (1 + A + B + C)

    x0u0_pred <- 1 / denom
    x1u0_pred <- A / denom
    x0u1_pred <- B / denom
    x1u1_pred <- C / denom
    XUpred <- data.frame(x0u0_pred, x1u0_pred, x0u1_pred, x1u1_pred)
    XUpred4 <- bind_rows(XUpred, XUpred, XUpred, XUpred)

    combined <- bind_rows(df, df, df, df) %>%
      mutate(Xbar = rep(c(1, 0, 1, 0), each = n),
             Ubar = rep(c(1, 1, 0, 0), each = n),
             pXU = case_when(Xbar == 0 & Ubar == 0 ~ XUpred4[, 1],
                             Xbar == 1 & Ubar == 0 ~ XUpred4[, 2],
                             Xbar == 0 & Ubar == 1 ~ XUpred4[, 3],
                             Xbar == 1 & Ubar == 1 ~ XUpred4[, 4]))

    Final <- glm(Y ~ Xbar + Ubar, family = binomial(link = "logit"), weights = pXU, data = combined)
    est <- summary(Final)$coef[2, 1]
    se <- summary(Final)$coef[2, 2]
    alpha <- 1 - level
    return(list(exp(est), c(exp(est + se * qnorm(alpha / 2)), exp(est + se * qnorm(1 - alpha / 2)))))

  }

  if (c == 1) {

    C <- data[,confounders]

    df <- data.frame(Xstar, Y, C)

    x1u0_c <- px1_u0_parameters[4]
    x0u1_c <- px0_u1_parameters[4]
    x1u1_c <- px1_u1_parameters[4]

    A <- exp(x1u0_0 + x1u0_xstar * df$Xstar + x1u0_y * df$Y + x1u0_c * df$C)
    B <- exp(x0u1_0 + x0u1_xstar * df$Xstar + x0u1_y * df$Y + x0u1_c * df$C)
    C <- exp(x1u1_0 + x1u1_xstar * df$Xstar + x1u1_y * df$Y + x1u1_c * df$C)

    denom <- (1 + A + B + C)

    x0u0_pred <- 1 / denom
    x1u0_pred <- A / denom
    x0u1_pred <- B / denom
    x1u1_pred <- C / denom
    XUpred <- data.frame(x0u0_pred, x1u0_pred, x0u1_pred, x1u1_pred)
    XUpred4 <- bind_rows(XUpred, XUpred, XUpred, XUpred)

    combined <- bind_rows(df, df, df, df) %>%
      mutate(Xbar = rep(c(1, 0, 1, 0), each = n),
             Ubar = rep(c(1, 1, 0, 0), each = n),
             pXU = case_when(Xbar == 0 & Ubar == 0 ~ XUpred4[, 1],
                             Xbar == 1 & Ubar == 0 ~ XUpred4[, 2],
                             Xbar == 0 & Ubar == 1 ~ XUpred4[, 3],
                             Xbar == 1 & Ubar == 1 ~ XUpred4[, 4]))

    Final <- glm(Y ~ Xbar + C + Ubar, family = binomial(link = "logit"), weights = pXU, data = combined)
    est <- summary(Final)$coef[2, 1]
    se <- summary(Final)$coef[2, 2]
    alpha <- 1 - level
    return(list(exp(est), c(exp(est + se * qnorm(alpha / 2)), exp(est + se * qnorm(1 - alpha / 2)))))

  }

  else if (c == 2) {

    C1 <- data[,confounders[1]]
    C2 <- data[,confounders[2]]

    df <- data.frame(Xstar, Y, C1, C2)

    x1u0_c1 <- px1_u0_parameters[4]
    x1u0_c2 <- px1_u0_parameters[5]

    x0u1_c1 <- px0_u1_parameters[4]
    x0u1_c2 <- px0_u1_parameters[5]

    x1u1_c1 <- px1_u1_parameters[4]
    x1u1_c2 <- px1_u1_parameters[5]

    A <- exp(x1u0_0 + x1u0_xstar * df$Xstar + x1u0_y * df$Y + x1u0_c1 * df$C1 + x1u0_c2 * df$C2)
    B <- exp(x0u1_0 + x0u1_xstar * df$Xstar + x0u1_y * df$Y + x0u1_c1 * df$C1 + x0u1_c2 * df$C2)
    C <- exp(x1u1_0 + x1u1_xstar * df$Xstar + x1u1_y * df$Y + x1u1_c1 * df$C1 + x1u1_c2 * df$C2)

    denom <- (1 + A + B + C)

    x0u0_pred <- 1 / denom
    x1u0_pred <- A / denom
    x0u1_pred <- B / denom
    x1u1_pred <- C / denom
    XUpred <- data.frame(x0u0_pred, x1u0_pred, x0u1_pred, x1u1_pred)
    XUpred4 <- bind_rows(XUpred, XUpred, XUpred, XUpred)

    combined <- bind_rows(df, df, df, df) %>%
      mutate(Xbar = rep(c(1, 0, 1, 0), each = n),
             Ubar = rep(c(1, 1, 0, 0), each = n),
             pXU = case_when(Xbar == 0 & Ubar == 0 ~ XUpred4[, 1],
                             Xbar == 1 & Ubar == 0 ~ XUpred4[, 2],
                             Xbar == 0 & Ubar == 1 ~ XUpred4[, 3],
                             Xbar == 1 & Ubar == 1 ~ XUpred4[, 4]))

    Final <- glm(Y ~ Xbar + C1 + C2 + Ubar, family = binomial(link = "logit"), weights = pXU, data = combined)
    est <- summary(Final)$coef[2, 1]
    se <- summary(Final)$coef[2, 2]
    alpha <- 1 - level
    return(list(exp(est), c(exp(est + se * qnorm(alpha / 2)), exp(est + se * qnorm(1 - alpha / 2)))))

  }
  else if (c == 3) {

    C1 <- data[,confounders[1]]
    C2 <- data[,confounders[2]]
    C3 <- data[,confounders[3]]

    df <- data.frame(Xstar, Y, C1, C2, C3)

    x1u0_c1 <- px1_u0_parameters[4]
    x1u0_c2 <- px1_u0_parameters[5]
    x1u0_c3 <- px1_u0_parameters[6]

    x0u1_c1 <- px0_u1_parameters[4]
    x0u1_c2 <- px0_u1_parameters[5]
    x0u1_c3 <- px0_u1_parameters[6]

    x1u1_c1 <- px1_u1_parameters[4]
    x1u1_c2 <- px1_u1_parameters[5]
    x1u1_c3 <- px1_u1_parameters[6]

    A <- exp(x1u0_0 + x1u0_xstar * df$Xstar + x1u0_y * df$Y + x1u0_c1 * df$C1 + x1u0_c2 * df$C2 + x1u0_c3 * df$C3)
    B <- exp(x0u1_0 + x0u1_xstar * df$Xstar + x0u1_y * df$Y + x0u1_c1 * df$C1 + x0u1_c2 * df$C2 + x0u1_c3 * df$C3)
    C <- exp(x1u1_0 + x1u1_xstar * df$Xstar + x1u1_y * df$Y + x1u1_c1 * df$C1 + x1u1_c2 * df$C2 + x1u1_c3 * df$C3)

    denom <- (1 + A + B + C)

    x0u0_pred <- 1 / denom
    x1u0_pred <- A / denom
    x0u1_pred <- B / denom
    x1u1_pred <- C / denom
    XUpred <- data.frame(x0u0_pred, x1u0_pred, x0u1_pred, x1u1_pred)
    XUpred4 <- bind_rows(XUpred, XUpred, XUpred, XUpred)

    combined <- bind_rows(df, df, df, df) %>%
      mutate(Xbar = rep(c(1, 0, 1, 0), each = n),
             Ubar = rep(c(1, 1, 0, 0), each = n),
             pXU = case_when(Xbar == 0 & Ubar == 0 ~ XUpred4[, 1],
                             Xbar == 1 & Ubar == 0 ~ XUpred4[, 2],
                             Xbar == 0 & Ubar == 1 ~ XUpred4[, 3],
                             Xbar == 1 & Ubar == 1 ~ XUpred4[, 4]))

    Final <- glm(Y ~ Xbar + C1 + C2 + C3 + Ubar, family = binomial(link = "logit"), weights = pXU, data = combined)
    est <- summary(Final)$coef[2, 1]
    se <- summary(Final)$coef[2, 2]
    alpha <- 1 - level
    return(list(exp(est), c(exp(est + se * qnorm(alpha / 2)), exp(est + se * qnorm(1 - alpha / 2)))))

  }

  if (c > 3) {

    stop('This function is currently not compatible with >3 confounders')

  }

}


