#' Adust for uncontrolled confounding, outcome misclassification, and selection
#' bias.
#'
#' \code{adjust_multinom_uc_omc_sel} returns the exposure-outcome odds ratio and
#' confidence interval, adjusted for uncontrolled confounding, outcome
#' misclassificaiton, and selection bias.
#'
#' This function uses one bias model, a multinomial logistic regression model,
#' to predict the uncontrolled confounder (\emph{U}) and outcome (\emph{Y}).
#' If separate bias models for \emph{U} and \emph{Y} are desired, use
#' \code{adjust_uc_omc_sel}.
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
#' model:
#' \ifelse{html}{\out{log(P(U=1,Y=0)/P(U=0,Y=0)) = &gamma;<sub>1,0</sub> + &gamma;<sub>1,1</sub>X + &gamma;<sub>1,2</sub>Y* + &gamma;<sub>1,2+j</sub>C<sub>j</sub>, }}{\eqn{log(P(U=1,Y=0)/P(U=0,Y=0)) = \gamma_{1,0} + \gamma_{1,1} X + \gamma_{1,2} Y^* + \gamma_{1,2+j} C_j, }}
#' where \emph{U} is the binary unmeasured confounder,
#' \emph{Y} is the binary true outcome,
#' \emph{X} is the exposure, \emph{Y*} is the binary misclassified outcome,
#' \emph{C} represents the vector of measured confounders (if any), and
#' \emph{j} corresponds to the number of measured confounders.
#' @param u0y1_model_coefs The regression coefficients corresponding to the
#' model:
#' \ifelse{html}{\out{log(P(U=0,Y=1)/P(U=0,Y=0)) = &gamma;<sub>2,0</sub> + &gamma;<sub>2,1</sub>X + &gamma;<sub>2,2</sub>Y* + &gamma;<sub>2,2+j</sub>C<sub>j</sub>, }}{\eqn{log(P(U=0,Y=1)/P(U=0,Y=0)) = \gamma_{2,0} + \gamma_{2,1} X + \gamma_{2,2} Y^* + \gamma_{2,2+j} C_j, }}
#' where \emph{U} is the binary unmeasured confounder,
#' \emph{Y} is the binary true outcome,
#' \emph{X} is the exposure, \emph{Y*} is the binary misclassified outcome,
#' \emph{C} represents the vector of measured confounders (if any), and
#' \emph{j} corresponds to the number of measured confounders.
#' @param u1y1_model_coefs The regression coefficients corresponding to the
#' model:
#' \ifelse{html}{\out{log(P(U=1,Y=1)/P(U=0,Y=0)) = &gamma;<sub>3,0</sub> + &gamma;<sub>3,1</sub>X + &gamma;<sub>3,2</sub>Y* + &gamma;<sub>3,2+j</sub>C<sub>j</sub>, }}{\eqn{log(P(U=1,Y=1)/P(U=0,Y=0)) = \gamma_{3,0} + \gamma_{3,1} X + \gamma_{3,2} Y^* + \gamma_{3,2+j} C_j, }}
#' where \emph{U} is the binary unmeasured confounder,
#' \emph{Y} is the binary true outcome,
#' \emph{X} is the exposure, \emph{Y*} is the binary misclassified outcome,
#' \emph{C} represents the vector of measured confounders (if any), and
#' \emph{j} corresponds to the number of measured confounders.
#' @param s_model_coefs The regression coefficients corresponding to the model:
#' \ifelse{html}{\out{logit(P(S=1)) = &beta;<sub>0</sub> + &beta;<sub>1</sub>X + &beta;<sub>2</sub>Y* + &beta;<sub>2+j</sub>C<sub>j</sub>, }}{\eqn{logit(P(S=1)) = \beta_0 + \beta_1 X + \beta_2 Y^* + \beta_{2+j} C_j, }}
#' where \emph{S} represents binary selection,
#' \emph{X} is the exposure, \emph{Y*} is the binary misclassified outcome,
#' \emph{C} represents the vector of measured confounders (if any), and
#' \emph{j} corresponds to the number of measured confounders.
#' @return A list where the first item is the odds ratio estimate of the
#' effect of the exposure on the outcome and the second item is the
#' confidence interval as the vector: (lower bound, upper bound).
#'
#' @examples
#' adjust_multinom_uc_omc_sel(
#'   df_uc_omc_sel,
#'   exposure = "X",
#'   outcome = "Ystar",
#'   confounders = c("C1", "C2", "C3"),
#'   u1y0_model_coefs = c(-0.20, 0.62, 0.01, -0.08, 0.10, -0.15),
#'   u0y1_model_coefs = c(-3.28, 0.63, 1.65, 0.42, -0.85, 0.26),
#'   u1y1_model_coefs = c(-2.70, 1.22, 1.64, 0.32, -0.77, 0.09),
#'   s_model_coefs = c(0.00, 0.74, 0.19, 0.02, -0.06, 0.02)
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

adjust_multinom_uc_omc_sel <- function(
  data,
  exposure,
  outcome,
  confounders = NULL,
  u0y1_model_coefs,
  u1y0_model_coefs,
  u1y1_model_coefs,
  s_model_coefs,
  level = 0.95
) {

  n <- nrow(data)
  len_c <- length(confounders)
  len_u0y1_coefs <- length(u0y1_model_coefs)
  len_u1y0_coefs <- length(u1y0_model_coefs)
  len_u1y1_coefs <- length(u1y1_model_coefs)
  len_s_coefs <- length(s_model_coefs)

  x <- data[, exposure]
  ystar <- data[, outcome]

  if (sum(ystar %in% c(0, 1)) != n) {
    stop("Outcome must be a binary integer.")
  }
  if (len_u0y1_coefs != 3 + len_c) {
    stop(
      paste0(
        "Incorrect length of U0Y1 model coefficients. ",
        "Length should equal 3 + number of confounders."
      )
    )
  }
  if (len_u1y0_coefs != 3 + len_c) {
    stop(
      paste0(
        "Incorrect length of U1Y0 model coefficients. ",
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
  if (len_s_coefs != 3 + len_c) {
    stop(
      paste0(
        "Incorrect length of S model coefficients. ",
        "Length should equal 3 + number of confounders."
      )
    )
  }

  s1_0       <- s_model_coefs[1]
  s1_x       <- s_model_coefs[2]
  s1_ystar   <- s_model_coefs[3]

  u0y1_0     <- u0y1_model_coefs[1]
  u0y1_x     <- u0y1_model_coefs[2]
  u0y1_ystar <- u0y1_model_coefs[3]

  u1y0_0     <- u1y0_model_coefs[1]
  u1y0_x     <- u1y0_model_coefs[2]
  u1y0_ystar <- u1y0_model_coefs[3]

  u1y1_0     <- u1y1_model_coefs[1]
  u1y1_x     <- u1y1_model_coefs[2]
  u1y1_ystar <- u1y1_model_coefs[3]

  if (is.null(confounders)) {

    df <- data.frame(X = x, Ystar = ystar)

    p_u0y1 <- exp(u0y1_0 + u0y1_x * df$X + u0y1_ystar * df$Ystar)
    p_u1y0 <- exp(u1y0_0 + u1y0_x * df$X + u1y0_ystar * df$Ystar)
    p_u1y1 <- exp(u1y1_0 + u1y1_x * df$X + u1y1_ystar * df$Ystar)

    denom <- (1 + p_u0y1 + p_u1y0 + p_u1y1)

    u0y0_pred <- 1 / denom
    u0y1_pred <- p_u0y1 / denom
    u1y0_pred <- p_u1y0 / denom
    u1y1_pred <- p_u1y1 / denom

    df_uy_pred <- data.frame(
      U0Y0 = u0y0_pred,
      U0Y1 = u0y1_pred,
      U1Y0 = u1y0_pred,
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
                             Ubar == 1 & Ybar == 1 ~ U1Y1),
             pS = plogis(s1_0 + s1_x * .data$X + s1_ystar * .data$Ystar))

    suppressWarnings({
      final <- glm(
        Ybar ~ X + Ubar,
        family = binomial(link = "logit"),
        weights = (combined$pUY / combined$pS),
        data = combined
      )
    })

  } else if (len_c == 1) {

    c1 <- data[, confounders]
    df <- data.frame(X = x, Ystar = ystar, C1 = c1)

    s1_c1   <- s_model_coefs[4]
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
                             Ubar == 1 & Ybar == 1 ~ U1Y1),
             pS = plogis(s1_0 + s1_x * .data$X + s1_ystar * .data$Ystar +
                           s1_c1 * .data$C1))

    suppressWarnings({
      final <- glm(
        Ybar ~ X + C1 + Ubar,
        family = binomial(link = "logit"),
        weights = (combined$pUY / combined$pS),
        data = combined
      )
    })

  } else if (len_c == 2) {

    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]

    df <- data.frame(X = x, Ystar = ystar, C1 = c1, C2 = c2)

    s1_c1   <- s_model_coefs[4]
    s1_c2   <- s_model_coefs[5]

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
                             Ubar == 1 & Ybar == 1 ~ U1Y1),
             pS = plogis(s1_0 + s1_x * .data$X + s1_ystar * .data$Ystar +
                           s1_c1 * .data$C1 + s1_c2 * .data$C2))

    suppressWarnings({
      final <- glm(
        Ybar ~ X + C1 + C2 + Ubar,
        family = binomial(link = "logit"),
        weights = (combined$pUY / combined$pS),
        data = combined
      )
    })

  } else if (len_c == 3) {

    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]
    c3 <- data[, confounders[3]]

    df <- data.frame(X = x, Ystar = ystar, C1 = c1, C2 = c2, C3 = c3)

    s1_c1   <- s_model_coefs[4]
    s1_c2   <- s_model_coefs[5]
    s1_c3   <- s_model_coefs[6]

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
                             Ubar == 1 & Ybar == 1 ~ U1Y1),
             pS = plogis(s1_0 + s1_x * .data$X + s1_ystar * .data$Ystar +
                           s1_c1 * .data$C1 + s1_c2 * .data$C2 +
                           s1_c3 * .data$C3))

    suppressWarnings({
      final <- glm(
        Ybar ~ X + C1 + C2 + C3 + Ubar,
        family = binomial(link = "logit"),
        weights = (combined$pUY / combined$pS),
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