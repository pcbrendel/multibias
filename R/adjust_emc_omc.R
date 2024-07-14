#' Adust for exposure misclassification and outcome misclassification.
#'
#' \code{adjust_emc_omc} returns the exposure-outcome odds ratio and confidence
#' interval, adjusted for exposure misclassification and outcome
#' misclassification.
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
#' @inheritParams adjust_uc
#' @param x_model_coefs The regression coefficients corresponding to the model:
#'  \ifelse{html}{\out{logit(P(X=1)) = &delta;<sub>0</sub> + &delta;<sub>1</sub>X* + &delta;<sub>2</sub>Y* + &delta;<sub>2+j</sub>C<sub>j</sub>, }}{\eqn{logit(P(X=1)) = \delta_0 + \delta_1 X^* + \delta_2 Y^* + \delta{2+j} C_j, }}
#'  where \emph{X} represents the binary true exposure, \emph{X*} is the
#'  binary misclassified exposure, \emph{Y*} is the binary misclassified
#'  outcome, \emph{C} represents the vector of
#'  measured confounders (if any), and \emph{j} corresponds to the number
#'  of measured confounders. The number of parameters is therefore 3 + \emph{j}.
#' @param y_model_coefs The regression coefficients corresponding to the model:
#'  \ifelse{html}{\out{logit(P(Y=1)) = &beta;<sub>0</sub> + &beta;<sub>1</sub>X + &beta;<sub>2</sub>Y* + &beta;<sub>2+j</sub>C<sub>j</sub>, }}{\eqn{logit(P(Y=1)) = \beta_0 + \beta_1 X + \beta_2 Y^* + \beta{{2+j}} C_j, }}
#'  where \emph{Y} represents the binary true exposure,
#'  \emph{X} is the binary exposure, \emph{Y} is the binary
#'  misclassified outcome, \emph{C} represents the vector of measured
#'  confounders (if any), and \emph{j} corresponds to the number of measured
#'  confounders. The number of parameters is therefore 3 + \emph{j}.
#' @return A list where the first item is the odds ratio estimate of the
#'  effect of the exposure on the outcome and the second item is the
#'  confidence interval as the vector: (lower bound, upper bound).
#'
#' @examples
#' adjust_emc_omc(
#'   df_emc_omc,
#'   exposure = "Xstar",
#'   outcome = "Ystar",
#'   confounders = "C1",
#'   x_model_coefs = c(-2.15, 1.64, 0.35, 0.38),
#'   y_model_coefs = c(-3.10, 0.63, 1.60, 0.39)
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

adjust_emc_omc <- function(
  data,
  exposure,
  outcome,
  confounders = NULL,
  x_model_coefs,
  y_model_coefs,
  level = 0.95
) {

  n <- nrow(data)
  len_c <- length(confounders)
  len_x_coefs <- length(x_model_coefs)
  len_y_coefs <- length(y_model_coefs)

  xstar <- data[, exposure]
  ystar <- data[, outcome]

  if (sum(xstar %in% c(0, 1)) != n) {
    stop("Exposure must be a binary integer.")
  }
  if (sum(ystar %in% c(0, 1)) != n) {
    stop("Outcome must be a binary integer.")
  }
  if (len_x_coefs != 3 + len_c) {
    stop(
      paste0(
        "Incorrect length of X model coefficients. ",
        "Length should equal 3 + number of confounders."
      )
    )
  }
  if (len_y_coefs != 3 + len_c) {
    stop(
      paste0(
        "Incorrect length of Y model coefficients. ",
        "Length should equal 3 + number of confounders."
      )
    )
  }

  x1_0     <- x_model_coefs[1]
  x1_xstar <- x_model_coefs[2]
  x1_ystar <- x_model_coefs[3]

  y1_0     <- y_model_coefs[1]
  y1_x     <- y_model_coefs[2]
  y1_ystar <- y_model_coefs[3]

  if (is.null(confounders)) {

    df <- data.frame(Xstar = xstar, Ystar = ystar)
    df$Xpred <- rbinom(n, 1, plogis(x1_0 + x1_xstar * df$Xstar +
                                      x1_ystar * df$Ystar))
    df$Ypred <- rbinom(n, 1, plogis(y1_0 + y1_x * df$Xpred +
                                      y1_ystar * df$Ystar))

    final <- glm(
      Ypred ~ Xpred,
      family = binomial(link = "logit"),
      data = df
    )

  } else if (len_c == 1) {

    c1 <- data[, confounders]
    df <- data.frame(Xstar = xstar, Ystar = ystar, C1 = c1)

    x1_c1 <- x_model_coefs[4]
    y1_c1 <- y_model_coefs[4]

    df$Xpred <- rbinom(n, 1, plogis(x1_0 + x1_xstar * df$Xstar +
                                      x1_ystar * df$Ystar + x1_c1 * df$C1))
    df$Ypred <- rbinom(n, 1, plogis(y1_0 + y1_x * df$Xpred +
                                      y1_ystar * df$Ystar + y1_c1 * df$C1))

    final <- glm(
      Ypred ~ Xpred + C1,
      family = binomial(link = "logit"),
      data = df
    )

  } else if (len_c == 2) {

    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]

    df <- data.frame(Xstar = xstar, Ystar = ystar, C1 = c1, C2 = c2)

    x1_c1 <- x_model_coefs[4]
    x1_c2 <- x_model_coefs[5]

    y1_c1 <- y_model_coefs[4]
    y1_c2 <- y_model_coefs[5]

    df$Xpred <- rbinom(n, 1, plogis(x1_0 + x1_xstar * df$Xstar +
                                      x1_ystar * df$Ystar +
                                      x1_c1 * df$C1 + x1_c2 * df$C2))
    df$Ypred <- rbinom(n, 1, plogis(y1_0 + y1_x * df$Xpred +
                                      y1_ystar * df$Ystar +
                                      y1_c1 * df$C1 + y1_c2 * df$C2))

    final <- glm(
      Ypred ~ Xpred + C1 + C2,
      family = binomial(link = "logit"),
      data = df
    )

  } else if (len_c == 3) {

    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]
    c3 <- data[, confounders[3]]

    df <- data.frame(Xstar = xstar, Ystar = ystar, C1 = c1, C2 = c2, C3 = c3)

    x1_c1 <- x_model_coefs[4]
    x1_c2 <- x_model_coefs[5]
    x1_c3 <- x_model_coefs[6]

    y1_c1 <- y_model_coefs[4]
    y1_c2 <- y_model_coefs[5]
    y1_c3 <- y_model_coefs[6]

    df$Xpred <- rbinom(n, 1, plogis(x1_0 + x1_xstar * df$Xstar +
                                      x1_ystar * df$Ystar + x1_c1 * df$C1 +
                                      x1_c2 * df$C2 + x1_c3 * df$C3))
    df$Ypred <- rbinom(n, 1, plogis(y1_0 + y1_x * df$Xpred +
                                      y1_ystar * df$Ystar +
                                      y1_c1 * df$C1 + y1_c2 * df$C2))

    final <- glm(
      Ypred ~ Xpred + C1 + C2 + C3,
      family = binomial(link = "logit"),
      data = df
    )

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