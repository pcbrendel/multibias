#' Adust for uncontrolled confounding and outcome misclassification.
#'
#' \code{adjust_multinom_uc_omc} returns the exposure-outcome odds ratio and
#' confidence interval, adjusted for uncontrolled confounding and outcome
#' misclassificaiton.
#'
#' This function uses one bias model, a multinomial logistic regression model,
#' to predict the uncontrolled confounder (U) and outcome (Y). If separate bias
#' models for X and U are desired, use \code{adjust_uc_omc}.
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
#' @param u1y0_model_coefs The regression coefficients corresponding to the
#'  model: \ifelse{html}{\out{log(P(U=1,Y=0)/P(U=0,Y=0)) =
#'  &gamma;<sub>1,0</sub> + &gamma;<sub>1,1</sub>X +
#'  &gamma;<sub>1,2</sub>Y* + &gamma;<sub>1,2+j</sub>C<sub>j</sub>, }
#'  where U is the binary unmeasured confounder, Y is the binary true outcome,
#'  X is the binary exposure, Y* is the binary misclassified
#'  outcome, C represents the vector of binary measured confounders (if any),
#'  and j corresponds to the number of measured
#'  confounders.}{\eqn{log(P(U=1,Y=0)/P(U=0,Y=0)) =}}
#' @param u0y1_model_coefs The regression coefficients corresponding to the
#'  model: \ifelse{html}{\out{log(P(U=0,Y=1)/P(U=0,Y=0)) =
#'  &gamma;<sub>2,0</sub> + &gamma;<sub>2,1</sub>X +
#'  &gamma;<sub>2,2</sub>Y* + &gamma;<sub>2,2+j</sub>C<sub>j</sub>, }
#'  where U is the binary unmeasured confounder, Y is the binary true outcome,
#'  X is the binary exposure, Y* is the binary misclassified
#'  outcome, C represents the vector of binary measured confounders (if any),
#'  and j corresponds to the number of measured
#'  confounders.}{\eqn{log(P(U=0,Y=1)/P(U=0,Y=0)) =}}
#' @param u1y1_model_coefs The regression coefficients corresponding to the
#'  model: \ifelse{html}{\out{log(P(U=1,Y=1)/P(U=0,Y=0)) =
#'  &gamma;<sub>3,0</sub> + &gamma;<sub>3,1</sub>X +
#'  &gamma;<sub>3,2</sub>Y* + &gamma;<sub>3,2+j</sub>C<sub>j</sub>, }
#'  where U is the binary unmeasured confounder, Y is the binary true outcome,
#'  X is the binary exposure, Y* is the binary misclassified
#'  outcome, C represents the vector of binary measured confounders (if any),
#'  and j corresponds to the number of measured
#'  confounders.}{\eqn{log(P(U=1,Y=1)/P(U=0,Y=0)) =}}
#' @return A list where the first item is the odds ratio estimate of the
#'  effect of the exposure on the outcome and the second item is the
#'  confidence interval as the vector: (lower bound, upper bound).
#'
#' @examples
#' adjust_multinom_uc_omc(
#'   df_uc_omc,
#'   "X",
#'   "Ystar",
#'   "C1",
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
#'
#' @export

adjust_multinom_uc_omc <- function(
  data,
  exposure,
  outcome,
  confounders = NULL,
  u1y0_model_coefs,
  u0y1_model_coefs,
  u1y1_model_coefs,
  level = 0.95
) {

  n <- nrow(data)
  len_c <- length(confounders)
  len_u1y0_coefs <- length(u1y0_model_coefs)
  len_u0y1_coefs <- length(u0y1_model_coefs)
  len_u1y1_coefs <- length(u1y1_model_coefs)

  x <- data[, exposure]
  ystar <- data[, outcome]

  if (sum(x %in% c(0, 1)) != n) {
    stop("Exposure must be a binary integer.")
  }
  if (sum(ystar %in% c(0, 1)) != n) {
    stop("Outcome must be a binary integer.")
  }
  if (len_u1y0_coefs != 3 + len_c) {
    stop(
      paste0(
        "Incorrect length of U1Y0 model coefficients. ",
        "Length should equal 3 + number of confounders."
      )
    )
  }
  if (len_u0y1_coefs != 3 + len_c) {
    stop(
      paste0(
        "Incorrect length of U0Y1 model coefficients. ",
        "Length should equal 3 + number of confounders."
      )
    )
  }
  if (len_u1y1_coefs != 3 + len_c) {
    stop(
      paste0(
        "Incorrect length of U1Y1 model coefficients. ",
        "Length should equal 3 + number of confounders."
      )
    )
  }

  u1y0_0     <- u1y0_model_coefs[1]
  u1y0_x     <- u1y0_model_coefs[2]
  u1y0_ystar <- u1y0_model_coefs[3]

  u0y1_0     <- u0y1_model_coefs[1]
  u0y1_x     <- u0y1_model_coefs[2]
  u0y1_ystar <- u0y1_model_coefs[3]

  u1y1_0     <- u1y1_model_coefs[1]
  u1y1_x     <- u1y1_model_coefs[2]
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
      mutate(Ubar = rep(c(1, 0, 1, 0), each = n),
             Ybar = rep(c(1, 1, 0, 0), each = n),
             pUY = case_when(Ubar == 0 & Ybar == 0 ~ U0Y0,
                             Ubar == 1 & Ybar == 0 ~ U1Y0,
                             Ubar == 0 & Ybar == 1 ~ U0Y1,
                             Ubar == 1 & Ybar == 1 ~ U1Y1))
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

    p_u1y0 <- exp(u1y0_0 + u1y0_x * df$X + u1y0_ystar * df$Ystar +
                    u1y0_c1 * df$C1)
    p_u0y1 <- exp(u0y1_0 + u0y1_x * df$X + u0y1_ystar * df$Ystar +
                    u0y1_c1 * df$C1)
    p_u1y1 <- exp(u1y1_0 + u1y1_x * df$X + u1y1_ystar * df$Ystar +
                    u1y1_c1 * df$C1)

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
      mutate(Ubar = rep(c(1, 0, 1, 0), each = n),
             Ybar = rep(c(1, 1, 0, 0), each = n),
             pUY = case_when(Ubar == 0 & Ybar == 0 ~ U0Y0,
                             Ubar == 1 & Ybar == 0 ~ U1Y0,
                             Ubar == 0 & Ybar == 1 ~ U0Y1,
                             Ubar == 1 & Ybar == 1 ~ U1Y1))

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

    p_u1y0 <- exp(u1y0_0 + u1y0_x * df$X + u1y0_ystar * df$Ystar +
                    u1y0_c1 * df$C1 + u1y0_c2 * df$C2)
    p_u0y1 <- exp(u0y1_0 + u0y1_x * df$X + u0y1_ystar * df$Ystar +
                    u0y1_c1 * df$C1 + u0y1_c2 * df$C2)
    p_u1y1 <- exp(u1y1_0 + u1y1_x * df$X + u1y1_ystar * df$Ystar +
                    u1y1_c1 * df$C1 + u1y1_c2 * df$C2)

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
      mutate(Ubar = rep(c(1, 0, 1, 0), each = n),
             Ybar = rep(c(1, 1, 0, 0), each = n),
             pUY = case_when(Ubar == 0 & Ybar == 0 ~ U0Y0,
                             Ubar == 1 & Ybar == 0 ~ U1Y0,
                             Ubar == 0 & Ybar == 1 ~ U0Y1,
                             Ubar == 1 & Ybar == 1 ~ U1Y1))

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

    p_u1y0 <- exp(u1y0_0 + u1y0_x * df$X + u1y0_ystar * df$Ystar +
                    u1y0_c1 * df$C1 + u1y0_c2 * df$C2 + u1y0_c3 * df$C3)
    p_u0y1 <- exp(u0y1_0 + u0y1_x * df$X + u0y1_ystar * df$Ystar +
                    u0y1_c1 * df$C1 + u0y1_c2 * df$C2 + u0y1_c3 * df$C3)
    p_u1y1 <- exp(u1y1_0 + u1y1_x * df$X + u1y1_ystar * df$Ystar +
                    u1y1_c1 * df$C1 + u1y1_c2 * df$C2 + u1y1_c3 * df$C3)

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
      mutate(Ubar = rep(c(1, 0, 1, 0), each = n),
             Ybar = rep(c(1, 1, 0, 0), each = n),
             pUY = case_when(Ubar == 0 & Ybar == 0 ~ U0Y0,
                             Ubar == 1 & Ybar == 0 ~ U1Y0,
                             Ubar == 0 & Ybar == 1 ~ U0Y1,
                             Ubar == 1 & Ybar == 1 ~ U1Y1))

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

  est <- summary(final)$coef[2, 1]
  se <- summary(final)$coef[2, 2]
  alpha <- 1 - level

  estimate <- exp(est)
  ci <- c(exp(est + se * qnorm(alpha / 2)),
          exp(est + se * qnorm(1 - alpha / 2)))
  return(list(estimate = estimate, ci = ci))

}