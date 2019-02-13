#' Adust for exposure misclassification and selection bias.
#'
#' \code{adjust_mc_sel} returns the exposure-outcome odds ratio and confidence
#' interval, adjusted for exposure misclassification and selection bias.
#'
#' @param data The data set.
#' @param exposure The variable corresponding to the exposure in the data.
#' @param outcome The variable corresponding to the outcome in the data.
#' @param confounders The variable(s) corresponding to the confounder(s) in the data.
#' A maximum of three confounders are allowed.
#' @param px1_parameters The regression coefficients corresponding to the model:
#' \ifelse{html}{\out{logit(P(X=1)) = &delta;<sub>0</sub> + &delta;<sub>1</sub>X* + &delta;<sub>2</sub>Y +
#' &delta;<sub>2+j</sub>C<sub>2+j</sub>, }where X represents (binary) true exposure, X* is the
#' (binary) misclassified exposure, Y is the (binary) outcome, C represents the vector of (binary)
#' measured confounders (if any), and j corresponds to the number of measured 
#' confounders. The number of parameters is therefore 3 + j.}{\eqn{logit(P(X=1)) =}}
#' @param ps1_parameters The regression coefficients corresponding to the model:
#' \ifelse{html}{\out{logit(P(S=1)) = &beta;<sub>0</sub> + &beta;<sub>1</sub>X* + &beta;<sub>2</sub>Y +
#' &beta;<sub>2+j</sub>C<sub>2+j</sub>, }where S represents (binary) selection, X* is the (binary) 
#' misclassified exposure, Y is the (binary) outcome, C represents the vector of (binary) 
#' measured confounders (if any), and j corresponds to the number of measured 
#' confounders. The number of parameters is therefore 3 + j.}{\eqn{logit(P(S=1)) =}}
#' @param level Number from 0-1 representing the range of the confidence interval. Default is .95.
#'
#' @examples
#' adjust_mc_sel(df_mc_sel, exposure = "Xstar", outcome = "Y",
#' confounders = c("C1", "C2", "C3"),px1_parameters =
#' c(-1.35, 1.64, 0.70, -0.42, -0.42, 0.40), ps1_parameters =
#' c(-0.08, 0.55, 2.02, -0.16, -0.14, 0.15))
#'
#' @import dplyr
#' @importFrom magrittr %>%
#'
#' @export

adjust_mc_sel <- function (data, exposure, outcome, confounders = NULL,
                           px1_parameters, ps1_parameters, level = .95) {

  n <- nrow(data)
  c <- length(confounders)
  p1 <- length(px1_parameters)
  p2 <- length(ps1_parameters)

  Xstar <- data[,exposure]
  Y <- data[,outcome]

  if (sum(Xstar %in% c(0, 1)) != n) {stop('Exposure must be binary')}
  if (sum(Y %in% c(0, 1)) != n) {stop('Outcome must be binary')}
  if (p1 != c + 3) {stop('Incorrect X1 parameter length')}
  if (p2 != c + 3) {stop('Incorrect S1 parameter length')}

  s1_0     <- ps1_parameters[1]
  s1_xstar <- ps1_parameters[2]
  s1_y     <- ps1_parameters[3]

  x1_0     <- px1_parameters[1]
  x1_xstar <- px1_parameters[2]
  x1_y     <- px1_parameters[3]

  if (is.null(confounders)) {

    df <- data.frame(Xstar, Y)

    x1_pred <- expit(x1_0 + x1_xstar * Xstar + x1_y * Y)
    x1_pred <- rep(x1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(Xbar = rep(c(1, 0), each = n),
             pS = expit(s1_0 + s1_xstar * Xstar + s1_y * Y),
             pX = case_when(Xbar == 1 ~ x1_pred,
                            Xbar == 0 ~ 1 - x1_pred))

    final <- glm(Y ~ Xbar, family = binomial(link = "logit"), weights = (pX / pS), data = combined)
    est <- summary(final)$coef[2, 1]
    se <- summary(final)$coef[2, 2]
    alpha <- 1 - level
    return(list(exp(est), c(exp(est + se * qnorm(alpha / 2)), exp(est + se * qnorm(1 - alpha / 2)))))

  }

  if (c == 1) {

    C <- data[,confounders]
    df <- data.frame(Xstar, Y, C)
    x1_c <- px1_parameters[4]

    x1_pred <- expit(x1_0 + x1_xstar * Xstar + x1_y * Y + x1_c * C)
    x1_pred <- rep(x1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(Xbar = rep(c(1, 0), each = n),
             pS = expit(s1_0 + s1_xstar * Xstar + s1_y * Y),
             pX = case_when(Xbar == 1 ~ x1_pred,
                            Xbar == 0 ~ 1 - x1_pred))

    final <- glm(Y ~ Xbar + C, family = binomial(link = "logit"), weights = (pX / pS), data = combined)
    est <- summary(final)$coef[2, 1]
    se <- summary(final)$coef[2, 2]
    alpha <- 1 - level
    return(list(exp(est), c(exp(est + se * qnorm(alpha / 2)), exp(est + se * qnorm(1 - alpha / 2)))))

  }

  if (c == 2) {

    C1 <- data[,confounders[1]]
    C2 <- data[,confounders[2]]

    df <- data.frame(Xstar, Y, C1, C2)

    x1_c1    <- px1_parameters[4]
    x1_c2    <- px1_parameters[5]

    x1_pred <- expit(x1_0 + x1_xstar * Xstar + x1_y * Y + x1_c1 * C1 + x1_c2 * C2)
    x1_pred <- rep(x1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(Xbar = rep(c(1, 0), each = n),
             pS = expit(s1_0 + s1_xstar * Xstar + s1_y * Y),
             pX = case_when(Xbar == 1 ~ x1_pred,
                            Xbar == 0 ~ 1 - x1_pred))

    final <- glm(Y ~ Xbar + C1 + C2, family = binomial(link = "logit"), weights = (pX / pS), data = combined)
    est <- summary(final)$coef[2, 1]
    se <- summary(final)$coef[2, 2]
    alpha <- 1 - level
    return(list(exp(est), c(exp(est + se * qnorm(alpha / 2)), exp(est + se * qnorm(1 - alpha / 2)))))

  }

  if (c == 3) {

    C1 <- data[,confounders[1]]
    C2 <- data[,confounders[2]]
    C3 <- data[,confounders[3]]

    df <- data.frame(Xstar, Y, C1, C2, C3)

    s1_c1    <- ps1_parameters[4]
    s1_c2    <- ps1_parameters[5]
    s1_c3    <- ps1_parameters[6]

    x1_c1    <- px1_parameters[4]
    x1_c2    <- px1_parameters[5]
    x1_c3    <- px1_parameters[6]

    x1_pred <- expit(x1_0 + x1_xstar * Xstar + x1_y * Y + x1_c1 * C1 + x1_c2 * C2 + x1_c3 * C3)
    x1_pred <- rep(x1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(Xbar = rep(c(1, 0), each = n),
             pS = expit(s1_0 + s1_xstar * Xstar + s1_y * Y + s1_c1 * C1 + s1_c2 * C2 + s1_c3 * C3),
             pX = case_when(Xbar == 1 ~ x1_pred,
                            Xbar == 0 ~ 1 - x1_pred))

    final <- glm(Y ~ Xbar + C1 + C2 + C3, family = binomial(link = "logit"), weights = (pX / pS), data = combined)
    est <- summary(final)$coef[2, 1]
    se <- summary(final)$coef[2, 2]
    alpha <- 1 - level
    return(list(exp(est), c(exp(est + se * qnorm(alpha / 2)), exp(est + se * qnorm(1 - alpha / 2)))))

  }

  if (c > 3) {

    stop('This function is currently not compatible with >3 confounders')

  }

}


