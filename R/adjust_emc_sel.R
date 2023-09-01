#' Adust for exposure misclassification and selection bias.
#'
#' \code{adjust_emc_sel} returns the exposure-outcome odds ratio and confidence
#' interval, adjusted for exposure misclassification and selection bias.
#'
#' Details
#'
#' @param data Dataframe for analysis.
#' @param exposure The variable corresponding to the exposure in the data.
#' @param outcome The variable corresponding to the outcome in the data.
#' @param confounders The variable(s) corresponding to the confounder(s) in the
#'  data. A maximum of three confounders are allowed.
#' @param px1_parameters The regression coefficients corresponding to the model:
#'  \ifelse{html}{\out{logit(P(X=1)) = &delta;<sub>0</sub> +
#'  &delta;<sub>1</sub>X* + &delta;<sub>2</sub>Y +
#'  &delta;<sub>2+j</sub>C<sub>j</sub>, }
#'  where X represents (binary) true exposure, X* is the (binary) misclassified
#'  exposure, Y is the (binary) outcome, C represents the vector of (binary)
#'  measured confounders (if any), and j corresponds to the number of measured
#'  confounders. The number of parameters is therefore
#'  3 + j.}{\eqn{logit(P(X=1)) =}}
#' @param ps1_parameters The regression coefficients corresponding to the model:
#'  \ifelse{html}{\out{logit(P(S=1)) = &beta;<sub>0</sub> +
#'  &beta;<sub>1</sub>X* + &beta;<sub>2</sub>Y +
#'  &beta;<sub>2+j</sub>C<sub>j</sub>, } where S represents (binary) selection,
#'  X* is the (binary) misclassified exposure, Y is the (binary) outcome,
#'  C represents the vector of (binary) measured confounders (if any), and j
#'  corresponds to the number of measured confounders. The number of
#'  parameters is therefore 3 + j.}{\eqn{logit(P(S=1)) =}}
#' @param level Value from 0-1 representing the full range of the confidence
#'  interval. Default is 0.95.
#'
#' @examples
#' adjust_emc_sel(
#'  df_mc_sel,
#'  exposure = "Xstar",
#'  outcome = "Y",
#'  confounders = c("C1", "C2", "C3"),
#'  px1_parameters = c(-1.35, 1.64, 0.70, -0.42, -0.42, 0.40),
#'  ps1_parameters = c(-0.08, 0.55, 2.02, -0.16, -0.14, 0.15)
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

adjust_emc_sel <- function(
  data,
  exposure,
  outcome,
  confounders = NULL,
  px1_parameters,
  ps1_parameters,
  level = 0.95
) {

  n <- nrow(data)
  len_c <- length(confounders)
  len_px1 <- length(px1_parameters)
  len_ps1 <- length(ps1_parameters)

  xstar <- data[, exposure]
  y <- data[, outcome]

  if (sum(xstar %in% c(0, 1)) != n) {
    stop("Exposure must be binary")
  }
  if (sum(y %in% c(0, 1)) != n) {
    stop("Outcome must be binary")
  }
  if (len_px1 != len_c + 3) {
    stop("Incorrect X1 parameter length")
  }
  if (len_ps1 != len_c + 3) {
    stop("Incorrect S1 parameter length")
  }

  s1_0     <- ps1_parameters[1]
  s1_xstar <- ps1_parameters[2]
  s1_y     <- ps1_parameters[3]

  x1_0     <- px1_parameters[1]
  x1_xstar <- px1_parameters[2]
  x1_y     <- px1_parameters[3]

  if (is.null(confounders)) {

    df <- data.frame(Xstar = xstar, Y = y)

    x1_pred <- plogis(x1_0 + x1_xstar * xstar + x1_y * y)
    x1_pred <- rep(x1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(Xbar = rep(c(1, 0), each = n),
             pS = plogis(s1_0 + s1_xstar * .data$Xstar + s1_y * .data$Y),
             pX = case_when(Xbar == 1 ~ x1_pred,
                            Xbar == 0 ~ 1 - x1_pred))

    final <- glm(
      Y ~ Xbar,
      family = binomial(link = "logit"),
      weights = (combined$pX / combined$pS),
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

  } else if (len_c == 1) {

    c1 <- data[, confounders]
    df <- data.frame(Xstar = xstar, Y = y, C1 = c1)
    x1_c1 <- px1_parameters[4]
    s1_c1    <- ps1_parameters[4]

    x1_pred <- plogis(x1_0 + x1_xstar * xstar + x1_y * y + x1_c1 * c1)
    x1_pred <- rep(x1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(
        Xbar = rep(c(1, 0), each = n),
        pS = plogis(s1_0 + s1_xstar * .data$Xstar + s1_y * .data$Y +
                      s1_c1 * .data$C1),
        pX = case_when(Xbar == 1 ~ x1_pred,
                       Xbar == 0 ~ 1 - x1_pred)
      )

    final <- glm(
      Y ~ Xbar + C1,
      family = binomial(link = "logit"),
      weights = (combined$pX / combined$pS),
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

    s1_c1    <- ps1_parameters[4]
    s1_c2    <- ps1_parameters[5]

    x1_c1    <- px1_parameters[4]
    x1_c2    <- px1_parameters[5]

    x1_pred <- plogis(x1_0 + x1_xstar * xstar +
                        x1_y * y + x1_c1 * c1 + x1_c2 * c2)
    x1_pred <- rep(x1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(
        Xbar = rep(c(1, 0), each = n),
        pS = plogis(s1_0 + s1_xstar * .data$Xstar + s1_y * .data$Y +
                      s1_c1 * .data$C1 + s1_c2 * .data$C2),
        pX = case_when(Xbar == 1 ~ x1_pred,
                       Xbar == 0 ~ 1 - x1_pred)
      )

    final <- glm(
      Y ~ Xbar + C1 + C2,
      family = binomial(link = "logit"),
      weights = (combined$pX / combined$pS),
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

    s1_c1    <- ps1_parameters[4]
    s1_c2    <- ps1_parameters[5]
    s1_c3    <- ps1_parameters[6]

    x1_c1    <- px1_parameters[4]
    x1_c2    <- px1_parameters[5]
    x1_c3    <- px1_parameters[6]

    x1_pred <- plogis(
      x1_0 + x1_xstar * xstar + x1_y * y + x1_c1 * c1 + x1_c2 * c2 + x1_c3 * c3
    )
    x1_pred <- rep(x1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(
        Xbar = rep(c(1, 0), each = n),
        pS = plogis(s1_0 + s1_xstar * .data$Xstar + s1_y * .data$Y +
                      s1_c1 * .data$C1 + s1_c2 * .data$C2 + s1_c3 * .data$C3),
        pX = case_when(Xbar == 1 ~ x1_pred,
                       Xbar == 0 ~ 1 - x1_pred)
      )

    final <- glm(
      Y ~ Xbar + C1 + C2 + C3,
      family = binomial(link = "logit"),
      weights = (combined$pX / combined$pS),
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