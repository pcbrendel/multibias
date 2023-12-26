#' Adust for exposure misclassification and outcome misclassification
#'
#' \code{adjust_multinom_emc_omc} returns the exposure-outcome odds ratio and
#' confidence interval, adjusted for exposure misclassificaiton and outcome
#' misclassification.
#'
#' This function uses one bias model, a multinomial logistic regression model,
#' to predict the exposure (X) and outcome (Y). If separate bias
#' models for X and Y are desired, use \code{adjust_emc_omc}.
#'
#' Values for the regression coefficients can be applied as
#' fixed values or as single draws from a probability
#' distribution (ex: \code{rnorm(1, mean = 2, sd = 1)}). The latter has
#' the advantage of allowing the researcher to capture the uncertainty
#' in the bias parameter estimates. To incorporate this uncertainty in the
#' estimate and confidence interval, this function should be run in loop across
#' bootstrap samples of the dataframe for analysis. The estimate and
#' confidence interval would then be obtained from the median and quantiles
#' of the distribution of odds ratio estimates.
#'
#' @inheritParams adjust_emc_sel
#' @param x1y0_model_coefs The regression coefficients corresponding to the
#'  model:
#'  \ifelse{html}{\out{log(P(X=1,Y=0) / P(X=0,Y=0)) = &gamma;<sub>1,0</sub> + &gamma;<sub>1,1</sub>X* + &gamma;<sub>1,2</sub>Y* + &gamma;<sub>1,2+j</sub>C<sub>j</sub>, }}{\eqn{log(P(X=1,Y=0) / P(X=0,Y=0)) = $\gamma_{1,0}$ + $\gamma_{1,1}$X* + $\gamma_{1,2}$Y* + $\gamma_{1,2+j}$C_j, }}
#'  where X is the binary true exposure, Y is the binary true outcome,
#'  X* is the binary misclassified exposure, Y* is the binary misclassified
#'  outcome, C represents the vector of binary measured confounders (if any),
#'  and j corresponds to the number of measured confounders.
#' @param x0y1_model_coefs The regression coefficients corresponding to the
#'  model:
#'  \ifelse{html}{\out{log(P(X=0,Y=1) / P(X=0,Y=0)) = &gamma;<sub>2,0</sub> + &gamma;<sub>2,1</sub>X* + &gamma;<sub>2,2</sub>Y* + &gamma;<sub>2,2+j</sub>C<sub>j</sub>, }}{\eqn{log(P(X=0,U=1) / P(X=0,U=0)) = $\gamma_{2,0}$ + $\gamma_{2,1}$X* + $\gamma_{2,2}$Y* + $\gamma_{2,2+j}$C_j, }}
#'  where X is the binary true exposure, Y is the binary true outcome,
#'  X* is the binary misclassified exposure, Y* is the binary misclassified
#'  outcome, C represents the vector of binary measured confounders (if any),
#'  and j corresponds to the number of measured confounders.
#' @param x1y1_model_coefs The regression coefficients corresponding to the
#'  model:
#'  \ifelse{html}{\out{log(P(X=1,Y=1) / P(X=0,Y=0)) = &gamma;<sub>3,0</sub> + &gamma;<sub>3,1</sub>X* + &gamma;<sub>3,2</sub>Y* + &gamma;<sub>3,2+j</sub>C<sub>j</sub>, }}{\eqn{log(P(X=1,Y=1) / P(X=0,Y=0)) = $\gamma_{3,0}$ + $\gamma_{3,1}$X* + $\gamma_{3,2}$Y* + $\gamma_{3,2+j}$C_j, }}
#'  where X is the binary true exposure, Y is the binary true outcome,
#'  X* is the binary misclassified exposure, Y* is the binary misclassified
#'  outcome, C represents the vector of binary measured confounders (if any),
#'  and j corresponds to the number of measured confounders.
#' @return A list where the first item is the odds ratio estimate of the
#'  effect of the exposure on the outcome and the second item is the
#'  confidence interval as the vector: (lower bound, upper bound).
#'
#' @examples
#' adjust_multinom_emc_omc(
#'   df_emc_omc,
#'   exposure = "Xstar",
#'   outcome = "Ystar",
#'   confounders = c("C1", "C2", "C3"),
#'   x1y0_model_coefs = c(-2.86, 1.63, 0.23, 0.37, -0.22, 0.87),
#'   x0y1_model_coefs = c(-3.26, 0.22, 1.60, 0.41, -0.93, 0.28),
#'   x1y1_model_coefs = c(-5.62, 1.83, 1.83, 0.74, -1.15, 1.19)
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

adjust_multinom_emc_omc <- function(
  data,
  exposure,
  outcome,
  confounders = NULL,
  x1y0_model_coefs,
  x0y1_model_coefs,
  x1y1_model_coefs,
  level = 0.95
) {

  n <- nrow(data)
  len_c <- length(confounders)
  len_x1y0_coefs <- length(x1y0_model_coefs)
  len_x0y1_coefs <- length(x0y1_model_coefs)
  len_x1y1_coefs <- length(x1y1_model_coefs)

  xstar <- data[, exposure]
  ystar <- data[, outcome]

  if (sum(xstar %in% c(0, 1)) != n) {
    stop("Exposure must be a binary integer.")
  }
  if (sum(ystar %in% c(0, 1)) != n) {
    stop("Outcome must be a binary integer.")
  }
  if (len_x1y0_coefs != 3 + len_c) {
    stop(
      paste0(
        "Incorrect length of X1Y0 model coefficients. ",
        "Length should equal 3 + number of confounders."
      )
    )
  }
  if (len_x0y1_coefs != 3 + len_c) {
    stop(
      paste0(
        "Incorrect length of X0Y1 model coefficients. ",
        "Length should equal 3 + number of confounders."
      )
    )
  }
  if (len_x1y1_coefs != 3 + len_c) {
    stop(
      paste0(
        "Incorrect length of X1Y1 model coefficients. ",
        "Length should equal 3 + number of confounders."
      )
    )
  }

  x1y0_0     <- x1y0_model_coefs[1]
  x1y0_xstar <- x1y0_model_coefs[2]
  x1y0_ystar <- x1y0_model_coefs[3]

  x0y1_0     <- x0y1_model_coefs[1]
  x0y1_xstar <- x0y1_model_coefs[2]
  x0y1_ystar <- x0y1_model_coefs[3]

  x1y1_0     <- x1y1_model_coefs[1]
  x1y1_xstar <- x1y1_model_coefs[2]
  x1y1_ystar <- x1y1_model_coefs[3]

  if (is.null(confounders)) {

    df <- data.frame(Xstar = xstar, Ystar = ystar)

    p_x1y0 <- exp(x1y0_0 + x1y0_xstar * df$Xstar + x1y0_ystar * df$Ystar)
    p_x0y1 <- exp(x0y1_0 + x0y1_xstar * df$Xstar + x0y1_ystar * df$Ystar)
    p_x1y1 <- exp(x1y1_0 + x1y1_xstar * df$Xstar + x1y1_ystar * df$Ystar)

    denom <- (1 + p_x1y0 + p_x0y1 + p_x1y1)

    x0y0_pred <- 1 / denom
    x1y0_pred <- p_x1y0 / denom
    x0y1_pred <- p_x0y1 / denom
    x1y1_pred <- p_x1y1 / denom

    df_xy_pred <- data.frame(
      X0Y0 = x0y0_pred,
      X1Y0 = x1y0_pred,
      X0Y1 = x0y1_pred,
      X1Y1 = x1y1_pred
    )
    df_xy_pred4 <- bind_rows(df_xy_pred, df_xy_pred, df_xy_pred, df_xy_pred)

    combined <- bind_rows(df, df, df, df) %>%
      bind_cols(df_xy_pred4) %>%
      mutate(Xbar = rep(c(1, 0, 1, 0), each = n),
             Ybar = rep(c(1, 1, 0, 0), each = n),
             pXY = case_when(Xbar == 0 & Ybar == 0 ~ X0Y0,
                             Xbar == 1 & Ybar == 0 ~ X1Y0,
                             Xbar == 0 & Ybar == 1 ~ X0Y1,
                             Xbar == 1 & Ybar == 1 ~ X1Y1))
    suppressWarnings({
      final <- glm(
        Ybar ~ Xbar,
        family = binomial(link = "logit"),
        weights = combined$pXY,
        data = combined
      )
    })

  } else if (len_c == 1) {

    c1 <- data[, confounders]

    df <- data.frame(Xstar = xstar, Ystar = ystar, C1 = c1)

    x1y0_c1 <- x1y0_model_coefs[4]
    x0y1_c1 <- x0y1_model_coefs[4]
    x1y1_c1 <- x1y1_model_coefs[4]

    p_x1y0 <- exp(x1y0_0 + x1y0_xstar * df$Xstar + x1y0_ystar * df$Ystar +
                    x1y0_c1 * df$C1)
    p_x0y1 <- exp(x0y1_0 + x0y1_xstar * df$Xstar + x0y1_ystar * df$Ystar +
                    x0y1_c1 * df$C1)
    p_x1y1 <- exp(x1y1_0 + x1y1_xstar * df$Xstar + x1y1_ystar * df$Ystar +
                    x1y1_c1 * df$C1)

    denom <- (1 + p_x1y0 + p_x0y1 + p_x1y1)

    x0y0_pred <- 1 / denom
    x1y0_pred <- p_x1y0 / denom
    x0y1_pred <- p_x0y1 / denom
    x1y1_pred <- p_x1y1 / denom

    df_xy_pred <- data.frame(
      X0Y0 = x0y0_pred,
      X1Y0 = x1y0_pred,
      X0Y1 = x0y1_pred,
      X1Y1 = x1y1_pred
    )
    df_xy_pred4 <- bind_rows(df_xy_pred, df_xy_pred, df_xy_pred, df_xy_pred)

    combined <- bind_rows(df, df, df, df) %>%
      bind_cols(df_xy_pred4) %>%
      mutate(Xbar = rep(c(1, 0, 1, 0), each = n),
             Ybar = rep(c(1, 1, 0, 0), each = n),
             pXY = case_when(Xbar == 0 & Ybar == 0 ~ X0Y0,
                             Xbar == 1 & Ybar == 0 ~ X1Y0,
                             Xbar == 0 & Ybar == 1 ~ X0Y1,
                             Xbar == 1 & Ybar == 1 ~ X1Y1))

    suppressWarnings({
      final <- glm(
        Ybar ~ Xbar + C1,
        family = binomial(link = "logit"),
        weights = combined$pXY,
        data = combined
      )
    })

  } else if (len_c == 2) {

    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]

    df <- data.frame(Xstar = xstar, Ystar = ystar, C1 = c1, C2 = c2)

    x1y0_c1 <- x1y0_model_coefs[4]
    x1y0_c2 <- x1y0_model_coefs[5]

    x0y1_c1 <- x0y1_model_coefs[4]
    x0y1_c2 <- x0y1_model_coefs[5]

    x1y1_c1 <- x1y1_model_coefs[4]
    x1y1_c2 <- x1y1_model_coefs[5]

    p_x1y0 <- exp(x1y0_0 + x1y0_xstar * df$Xstar + x1y0_ystar * df$Ystar +
                    x1y0_c1 * df$C1 + x1y0_c2 * df$C2)
    p_x0y1 <- exp(x0y1_0 + x0y1_xstar * df$Xstar + x0y1_ystar * df$Ystar +
                    x0y1_c1 * df$C1 + x0y1_c2 * df$C2)
    p_x1y1 <- exp(x1y1_0 + x1y1_xstar * df$Xstar + x1y1_ystar * df$Ystar +
                    x1y1_c1 * df$C1 + x1y1_c2 * df$C2)

    denom <- (1 + p_x1y0 + p_x0y1 + p_x1y1)

    x0y0_pred <- 1 / denom
    x1y0_pred <- p_x1y0 / denom
    x0y1_pred <- p_x0y1 / denom
    x1y1_pred <- p_x1y1 / denom

    df_xy_pred <- data.frame(
      X0Y0 = x0y0_pred,
      X1Y0 = x1y0_pred,
      X0Y1 = x0y1_pred,
      X1Y1 = x1y1_pred
    )
    df_xy_pred4 <- bind_rows(df_xy_pred, df_xy_pred, df_xy_pred, df_xy_pred)

    combined <- bind_rows(df, df, df, df) %>%
      bind_cols(df_xy_pred4) %>%
      mutate(Xbar = rep(c(1, 0, 1, 0), each = n),
             Ybar = rep(c(1, 1, 0, 0), each = n),
             pXY = case_when(Xbar == 0 & Ybar == 0 ~ X0Y0,
                             Xbar == 1 & Ybar == 0 ~ X1Y0,
                             Xbar == 0 & Ybar == 1 ~ X0Y1,
                             Xbar == 1 & Ybar == 1 ~ X1Y1))

    suppressWarnings({
      final <- glm(
        Ybar ~ Xbar + C1 + C2,
        family = binomial(link = "logit"),
        weights = combined$pXY,
        data = combined
      )
    })

  } else if (len_c == 3) {

    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]
    c3 <- data[, confounders[3]]

    df <- data.frame(Xstar = xstar, Ystar = ystar, C1 = c1, C2 = c2, C3 = c3)

    x1y0_c1 <- x1y0_model_coefs[4]
    x1y0_c2 <- x1y0_model_coefs[5]
    x1y0_c3 <- x1y0_model_coefs[6]

    x0y1_c1 <- x0y1_model_coefs[4]
    x0y1_c2 <- x0y1_model_coefs[5]
    x0y1_c3 <- x0y1_model_coefs[6]

    x1y1_c1 <- x1y1_model_coefs[4]
    x1y1_c2 <- x1y1_model_coefs[5]
    x1y1_c3 <- x1y1_model_coefs[6]

    p_x1y0 <- exp(x1y0_0 + x1y0_xstar * df$Xstar + x1y0_ystar * df$Ystar +
                    x1y0_c1 * df$C1 + x1y0_c2 * df$C2 + x1y0_c3 * df$C3)
    p_x0y1 <- exp(x0y1_0 + x0y1_xstar * df$Xstar + x0y1_ystar * df$Ystar +
                    x0y1_c1 * df$C1 + x0y1_c2 * df$C2 + x0y1_c3 * df$C3)
    p_x1y1 <- exp(x1y1_0 + x1y1_xstar * df$Xstar + x1y1_ystar * df$Ystar +
                    x1y1_c1 * df$C1 + x1y1_c2 * df$C2 + x1y1_c3 * df$C3)

    denom <- (1 + p_x1y0 + p_x0y1 + p_x1y1)

    x0y0_pred <- 1 / denom
    x1y0_pred <- p_x1y0 / denom
    x0y1_pred <- p_x0y1 / denom
    x1y1_pred <- p_x1y1 / denom

    df_xy_pred <- data.frame(
      X0Y0 = x0y0_pred,
      X1Y0 = x1y0_pred,
      X0Y1 = x0y1_pred,
      X1Y1 = x1y1_pred
    )
    df_xy_pred4 <- bind_rows(df_xy_pred, df_xy_pred, df_xy_pred, df_xy_pred)

    combined <- bind_rows(df, df, df, df) %>%
      bind_cols(df_xy_pred4) %>%
      mutate(Xbar = rep(c(1, 0, 1, 0), each = n),
             Ybar = rep(c(1, 1, 0, 0), each = n),
             pXY = case_when(Xbar == 0 & Ybar == 0 ~ X0Y0,
                             Xbar == 1 & Ybar == 0 ~ X1Y0,
                             Xbar == 0 & Ybar == 1 ~ X0Y1,
                             Xbar == 1 & Ybar == 1 ~ X1Y1))

    suppressWarnings({
      final <- glm(
        Ybar ~ Xbar + C1 + C2 + C3,
        family = binomial(link = "logit"),
        weights = combined$pXY,
        data = combined
      )
    })

  } else if (len_c > 3) {

    stop("This function is currently not compatible with >3 confounders.")

  }

  est <- summary(final)$coef[2, 1]
  se <- summary(final)$coef[2, 2]
  alpha <- 1 - level

  estimate <- exp(est)
  ci <- c(exp(est + se * qnorm(alpha / 2)),
          exp(est + se * qnorm(1 - alpha / 2)))
  return(list(estimate = estimate, ci = ci))

}