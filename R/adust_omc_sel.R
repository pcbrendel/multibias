#' Adust for outcome misclassification and selection bias.
#'
#' \code{adjust_omc_sel} returns the exposure-outcome odds ratio and confidence
#' interval, adjusted for outcome misclassification and selection bias.
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
#' @param y_model_coefs The regression coefficients corresponding to the model:
#'  \ifelse{html}{\out{logit(P(Y=1)) = &delta;<sub>0</sub> + &delta;<sub>1</sub>X + &delta;<sub>2</sub>Y* + &delta;<sub>2+j</sub>C<sub>j</sub>, }}{\eqn{logit(P(Y=1)) = \delta_0 + \delta_1 X + \delta_2 Y^* + \delta_{2+j} C_j, }}
#'  where Y represents the binary true outcome, X is the binary exposure,
#'  Y* is the binary misclassified outcome, C represents the vector of binary
#'  measured confounders (if any), and j corresponds to the number of measured
#'  confounders. The number of parameters is therefore 3 + j.
#' @param s_model_coefs The regression coefficients corresponding to the model:
#'  \ifelse{html}{\out{logit(P(S=1)) = &beta;<sub>0</sub> + &beta;<sub>1</sub>X + &beta;<sub>2</sub>Y* + &beta;<sub>2+j</sub>C<sub>j</sub>, }}{\eqn{logit(P(S=1)) = \beta_0 + \beta_1 X + \beta_2 Y^* + \beta_{2+j} C_j, }}
#'  where S represents binary selection,
#'  X is the binary exposure, Y* is the binary misclassified outcome,
#'  C represents the vector of binary measured confounders (if any), and j
#'  corresponds to the number of measured confounders. The number of
#'  parameters is therefore 3 + j.
#' @return A list where the first item is the odds ratio estimate of the
#'  effect of the exposure on the outcome and the second item is the
#'  confidence interval as the vector: (lower bound, upper bound).
#'
#' @examples
#' adjust_omc_sel(
#'   df_omc_sel,
#'   exposure = "X",
#'   outcome = "Ystar",
#'   confounders = "C1",
#'   y_model_coefs = c(-3.24, 0.58, 1.59, 0.45),
#'   s_model_coefs = c(0.03, 0.92, 0.12, 0.05)
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

adjust_omc_sel <- function(
  data,
  exposure,
  outcome,
  confounders = NULL,
  y_model_coefs,
  s_model_coefs,
  level = 0.95
) {

  n <- nrow(data)
  len_c <- length(confounders)
  len_y_coefs <- length(y_model_coefs)
  len_s_coefs <- length(s_model_coefs)

  x     <- data[, exposure]
  ystar <- data[, outcome]

  if (sum(x %in% c(0, 1)) != n) {
    stop("Exposure must be a binary integer.")
  }
  if (sum(ystar %in% c(0, 1)) != n) {
    stop("Outcome must be a binary integer.")
  }
  if (len_y_coefs != 3 + len_c) {
    stop(
      paste0(
        "Incorrect length of Y model coefficients. ",
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

  s1_0     <- s_model_coefs[1]
  s1_x     <- s_model_coefs[2]
  s1_ystar <- s_model_coefs[3]

  y1_0     <- y_model_coefs[1]
  y1_x     <- y_model_coefs[2]
  y1_ystar <- y_model_coefs[3]

  if (is.null(confounders)) {

    df <- data.frame(X = x, Ystar = ystar)

    y1_pred <- plogis(y1_0 + y1_x * x + y1_ystar * ystar)
    y1_pred <- rep(y1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(Ybar = rep(c(1, 0), each = n),
             pS = plogis(s1_0 + s1_x * .data$X + s1_ystar * .data$Ystar),
             pY = case_when(Ybar == 1 ~ y1_pred,
                            Ybar == 0 ~ 1 - y1_pred))

    suppressWarnings({
      final <- glm(
        Ybar ~ X,
        family = binomial(link = "logit"),
        weights = (combined$pY / combined$pS),
        data = combined
      )
    })

  } else if (len_c == 1) {

    c1 <- data[, confounders]
    df <- data.frame(X = x, Ystar = ystar, C1 = c1)
    y1_c1 <- y_model_coefs[4]
    s1_c1 <- s_model_coefs[4]

    y1_pred <- plogis(y1_0 + y1_x * x + y1_ystar * ystar + y1_c1 * c1)
    y1_pred <- rep(y1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(
        Ybar = rep(c(1, 0), each = n),
        pS = plogis(s1_0 + s1_x * .data$X + s1_ystar * .data$Ystar +
                      s1_c1 * .data$C1),
        pY = case_when(Ybar == 1 ~ y1_pred,
                       Ybar == 0 ~ 1 - y1_pred)
      )

    suppressWarnings({
      final <- glm(
        Ybar ~ X + C1,
        family = binomial(link = "logit"),
        weights = (combined$pY / combined$pS),
        data = combined
      )
    })

  } else if (len_c == 2) {

    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]

    df <- data.frame(X = x, Ystar = ystar, C1 = c1, C2 = c2)

    s1_c1 <- s_model_coefs[4]
    s1_c2 <- s_model_coefs[5]

    y1_c1 <- y_model_coefs[4]
    y1_c2 <- y_model_coefs[5]

    y1_pred <- plogis(y1_0 + y1_x * x +
                        y1_ystar * ystar + y1_c1 * c1 + y1_c2 * c2)
    y1_pred <- rep(y1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(
        Ybar = rep(c(1, 0), each = n),
        pS = plogis(s1_0 + s1_x * .data$X + s1_ystar * .data$Ystar +
                      s1_c1 * .data$C1 + s1_c2 * .data$C2),
        pY = case_when(Ybar == 1 ~ y1_pred,
                       Ybar == 0 ~ 1 - y1_pred)
      )

    suppressWarnings({
      final <- glm(
        Ybar ~ X + C1 + C2,
        family = binomial(link = "logit"),
        weights = (combined$pY / combined$pS),
        data = combined
      )
    })

  } else if (len_c == 3) {

    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]
    c3 <- data[, confounders[3]]

    df <- data.frame(X = x, Ystar = ystar, C1 = c1, C2 = c2, C3 = c3)

    s1_c1 <- s_model_coefs[4]
    s1_c2 <- s_model_coefs[5]
    s1_c3 <- s_model_coefs[6]

    y1_c1 <- y_model_coefs[4]
    y1_c2 <- y_model_coefs[5]
    y1_c3 <- y_model_coefs[6]

    y1_pred <- plogis(
      y1_0 + y1_x * x + y1_ystar * ystar + y1_c1 * c1 + y1_c2 * c2 + y1_c3 * c3
    )
    y1_pred <- rep(y1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(
        Ybar = rep(c(1, 0), each = n),
        pS = plogis(s1_0 + s1_x * .data$X + s1_ystar * .data$Ystar +
                      s1_c1 * .data$C1 + s1_c2 * .data$C2 + s1_c3 * .data$C3),
        pY = case_when(Ybar == 1 ~ y1_pred,
                       Ybar == 0 ~ 1 - y1_pred)
      )

    suppressWarnings({
      final <- glm(
        Ybar ~ X + C1 + C2 + C3,
        family = binomial(link = "logit"),
        weights = (combined$pY / combined$pS),
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