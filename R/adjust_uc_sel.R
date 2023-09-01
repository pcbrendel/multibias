#' Adust for uncontrolled confounding and selection bias.
#'
#'\code{adjust_uc_sel} returns the exposure-outcome odds ratio and confidence
#' interval, adjusted for uncontrolled confounding and exposure
#' misclassificaiton.
#'
#' Details
#'
#' @param data The data set.
#' @param exposure The name of the exposure variable in the data.
#' @param outcome The name of the outcome variable in the data.
#' @param confounders The variable name(s) of the confounder(s) in the data.
#'  A maximum of three confounders are allowed.
#' @param pu1_parameters The regression coefficients corresponding to the model:
#'  \ifelse{html}{\out{logit(P(U=1)) = &alpha;<sub>0</sub> +
#'  &alpha;<sub>1</sub>X + &alpha;<sub>2</sub>Y +
#'  &alpha;<sub>2+j</sub>C<sub>j</sub>, } where U is the (binary) unmeasured
#'  confounder, X is the (binary) exposure, Y is the (binary) outcome, C
#'  represents the vector of (binary) measured confounders (if any), and j
#'  corresponds to the number of measured confounders. The number of parameters
#'  therefore equals 3 + j.}{\eqn{logit(P(U=1)) =}}
#' @param ps1_parameters The regression coefficients corresponding to the model:
#'  \ifelse{html}{\out{logit(P(S=1)) = &beta;<sub>0</sub> +
#'  &beta;<sub>1</sub>X + &beta;<sub>2</sub>Y,} where S represents (binary)
#'  selection, X is the (binary) exposure, and Y is the (binary) outcome.
#'  The number of parameters therefore equals 3.}{\eqn{logit(P(S=1)) =}}
#' @param level Number from 0-1 representing the range of the confidence
#'  interval. Default is 0.95.
#'
#' @examples
#' adjust_uc_sel(
#'  df_uc_sel,
#'  exposure = "X",
#'  outcome = "Y",
#'  confounders = c("C1", "C2", "C3"),
#'  pu1_parameters = c(-.48, .41, .52, .12, .12, -.12),
#'  ps1_parameters = c(-.52, 2.01, 2.00)
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

adjust_uc_sel <- function(
  data,
  exposure,
  outcome,
  confounders = NULL,
  pu1_parameters,
  ps1_parameters,
  level = 0.95
) {

  n  <- nrow(data)
  len_c  <- length(confounders)
  len_pu1 <- length(pu1_parameters)
  len_ps1 <- length(ps1_parameters)

  x <- data[, exposure]
  y <- data[, outcome]

  if (sum(x %in% c(0, 1)) != n) {
    stop("Exposure must be binary.")
  }
  if (sum(y %in% c(0, 1)) != n) {
    stop("Outcome must be binary.")
  }
  if (len_pu1 != 3 + len_c) {
    stop(
      paste0(
        "Incorrect U1 parameter length.",
        "U1 parameter length should equal 3 + length(confounders)."
      )
    )
  }
  if (len_ps1 != 3) {
    stop("Incorrect S1 parameter length. S1 parameter length should equal 3.")
  }

  s1_0 <- ps1_parameters[1]
  s1_x <- ps1_parameters[2]
  s1_y <- ps1_parameters[3]

  u1_0 <- pu1_parameters[1]
  u1_x <- pu1_parameters[2]
  u1_y <- pu1_parameters[3]

  if (is.null(confounders)) {

    df <- data.frame(X = x, Y = y)

    u1_pred <- plogis(u1_0 + u1_x * x + u1_y * y)
    u1_pred <- rep(u1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(Ubar = rep(c(1, 0), each = n),
             pS = plogis(s1_0 + s1_x * .data$X + s1_y * .data$Y),
             pU = case_when(Ubar == 1 ~ u1_pred,
                            Ubar == 0 ~ 1 - u1_pred))

    final <- glm(
      Y ~ X + Ubar,
      family = binomial(link = "logit"),
      weights = (combined$pU / combined$pS),
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
    df <- data.frame(X = x, Y = y, C1 = c1)
    u1_c1 <- pu1_parameters[4]

    u1_pred <- plogis(u1_0 + u1_x * x + u1_y * y + u1_c1 * c1)
    u1_pred <- rep(u1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(Ubar = rep(c(1, 0), each = n),
             pS = plogis(s1_0 + s1_x * .data$X + s1_y * .data$Y),
             pU = case_when(Ubar == 1 ~ u1_pred,
                            Ubar == 0 ~ 1 - u1_pred))

    final <- glm(
      Y ~ X + C1 + Ubar,
      family = binomial(link = "logit"),
      weights = (combined$pU / combined$pS),
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

    df <- data.frame(X = x, Y = y, C1 = c1, C2 = c2)

    u1_c1 <- pu1_parameters[4]
    u1_c2 <- pu1_parameters[5]

    u1_pred <- plogis(u1_0 + u1_x * x + u1_y * y + u1_c1 * c1 + u1_c2 * c2)
    u1_pred <- rep(u1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(Ubar = rep(c(1, 0), each = n),
             pS = plogis(s1_0 + s1_x * .data$X + s1_y * .data$Y),
             pU = case_when(Ubar == 1 ~ u1_pred,
                            Ubar == 0 ~ 1 - u1_pred))

    final <- glm(
      Y ~ X + C1 + C2 + Ubar,
      family = binomial(link = "logit"),
      weights = (combined$pU / combined$pS),
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

    df <- data.frame(X = x, Y = y, C1 = c1, C2 = c2, C3 = c3)

    u1_c1 <- pu1_parameters[4]
    u1_c2 <- pu1_parameters[5]
    u1_c3 <- pu1_parameters[6]

    u1_pred <- plogis(u1_0 + u1_x * x + u1_y * y +
                        u1_c1 * c1 + u1_c2 * c2 + u1_c3 * c3)
    u1_pred <- rep(u1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(Ubar = rep(c(1, 0), each = n),
             pS = plogis(s1_0 + s1_x * .data$X + s1_y * .data$Y),
             pU = case_when(Ubar == 1 ~ u1_pred,
                            Ubar == 0 ~ 1 - u1_pred))

    final <- glm(
      Y ~ X + C1 + C2 + C3 + Ubar,
      family = binomial(link = "logit"),
      weights = (combined$pU / combined$pS),
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