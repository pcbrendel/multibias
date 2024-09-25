#' Adust for exposure misclassification.
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `adjust_emc()` was renamed to `adjust_em()`
#' @keywords internal
#'
#' @export
adjust_emc <- function(
    data_observed,
    x_model_coefs,
    level = 0.95) {
  lifecycle::deprecate_warn("1.5.3", "adjust_emc()", "adjust_em()")
  adjust_em(data_observed, x_model_coefs, level)
}


#' Adust for exposure misclassification.
#'
#' `adjust_em` returns the exposure-outcome odds ratio and confidence
#' interval, adjusted for exposure misclassificaiton.
#'
#' Values for the regression coefficients can be applied as
#' fixed values or as single draws from a probability
#' distribution (ex: `rnorm(1, mean = 2, sd = 1)`). The latter has
#' the advantage of allowing the researcher to capture the uncertainty
#' in the bias parameter estimates. To incorporate this uncertainty in the
#' estimate and confidence interval, this function should be run in loop across
#' bootstrap samples of the dataframe for analysis. The estimate and
#' confidence interval would then be obtained from the median and quantiles
#' of the distribution of odds ratio estimates.
#'
#' @param data_observed Object of class `data_observed` corresponding to the
#' data to perform bias analysis on.
#' @param x_model_coefs The regression coefficients corresponding to the model:
#' \ifelse{html}{\out{logit(P(X=1)) = &delta;<sub>0</sub> + &delta;<sub>1</sub>X* + &delta;<sub>2</sub>Y + &delta;<sub>2+j</sub>C<sub>j</sub>, }}{\eqn{logit(P(X=1)) = \delta_0 + \delta_1 X^* + \delta_2 Y + \delta_{2+j} C_j, }}
#' where *X* represents the binary true exposure, *X** is the binary
#' misclassified exposure, *Y* is the outcome, *C* represents
#' the vector of measured confounders (if any),
#' and *j* corresponds to the number of measured confounders. The number
#' of parameters is therefore 3 + *j*.
#' @param level Value from 0-1 representing the full range of the confidence
#'  interval. Default is 0.95.
#'
#' @return A list where the first item is the odds ratio estimate of the
#' effect of the exposure on the outcome and the second item is the
#' confidence interval as the vector: (lower bound, upper bound).
#'
#' @examples
#' df <- data_observed(
#'   data = df_em,
#'   exposure = "Xstar",
#'   outcome = "Y",
#'   confounders = "C1"
#' )
#' adjust_em(
#'   df,
#'   x_model_coefs = c(-2.10, 1.62, 0.63, 0.35)
#' )
#'
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom stats binomial
#' @importFrom stats glm
#' @importFrom stats lm
#' @importFrom stats qnorm
#' @importFrom stats rbinom
#' @importFrom stats plogis
#' @importFrom rlang .data
#'
#' @export

adjust_em <- function(
    data_observed,
    x_model_coefs,
    level = 0.95) {
  data <- data_observed$data
  n <- nrow(data)
  confounders <- data_observed$confounders
  len_c <- length(confounders)
  len_x_coefs <- length(x_model_coefs)

  xstar <- data[, data_observed$exposure]
  y <- data[, data_observed$outcome]

  force_binary(xstar, "Exposure must be a binary integer.")
  force_len(
    len_x_coefs,
    3 + len_c,
    paste0(
      "Incorrect length of X model coefficients. ",
      "Length should equal 3 + number of confounders."
    )
  )

  if (all(y %in% 0:1)) {
    y_binary <- TRUE
  } else {
    y_binary <- FALSE
  }

  x1_0 <- x_model_coefs[1]
  x1_xstar <- x_model_coefs[2]
  x1_y <- x_model_coefs[3]

  if (is.null(confounders)) {
    df <- data.frame(Xstar = xstar, Y = y)
    df$Xpred <- rbinom(n, 1, plogis(x1_0 + x1_xstar * df$Xstar + x1_y * df$Y))

    if (y_binary) {
      final <- glm(
        Y ~ Xpred,
        family = binomial(link = "logit"),
        data = df
      )
    } else {
      final <- lm(
        Y ~ Xpred,
        data = df
      )
    }
  } else if (len_c == 1) {
    c1 <- data[, confounders]
    df <- data.frame(Xstar = xstar, Y = y, C1 = c1)

    x1_c1 <- x_model_coefs[4]

    df$Xpred <- rbinom(
      n, 1, plogis(
        x1_0 + x1_xstar * df$Xstar +
          x1_y * df$Y + x1_c1 * df$C1
      )
    )

    if (y_binary) {
      final <- glm(
        Y ~ Xpred + C1,
        family = binomial(link = "logit"),
        data = df
      )
    } else {
      final <- lm(
        Y ~ Xpred + C1,
        data = df
      )
    }
  } else if (len_c == 2) {
    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]

    df <- data.frame(Xstar = xstar, Y = y, C1 = c1, C2 = c2)

    x1_c1 <- x_model_coefs[4]
    x1_c2 <- x_model_coefs[5]

    df$Xpred <- rbinom(
      n, 1, plogis(
        x1_0 + x1_xstar * df$Xstar + x1_y * df$Y +
          x1_c1 * df$C1 + x1_c2 * df$C2
      )
    )

    if (y_binary) {
      final <- glm(
        Y ~ Xpred + C1 + C2,
        family = binomial(link = "logit"),
        data = df
      )
    } else {
      final <- lm(
        Y ~ Xpred + C1 + C2,
        data = df
      )
    }
  } else if (len_c == 3) {
    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]
    c3 <- data[, confounders[3]]

    df <- data.frame(Xstar = xstar, Y = y, C1 = c1, C2 = c2, C3 = c3)

    x1_c1 <- x_model_coefs[4]
    x1_c2 <- x_model_coefs[5]
    x1_c3 <- x_model_coefs[6]

    df$Xpred <- rbinom(
      n, 1,
      plogis(
        x1_0 + x1_xstar * df$Xstar + x1_y * df$Y +
          x1_c1 * df$C1 + x1_c2 * df$C2 + x1_c3 * df$C3
      )
    )

    if (y_binary) {
      final <- glm(
        Y ~ Xpred + C1 + C2 + C3,
        family = binomial(link = "logit"),
        data = df
      )
    } else {
      final <- lm(
        Y ~ Xpred + C1 + C2 + C3,
        data = df
      )
    }
  } else if (len_c > 3) {
    stop("This function is currently not compatible with >3 confounders.")
  }

  est <- summary(final)$coef[2, 1]
  se <- summary(final)$coef[2, 2]
  alpha <- 1 - level

  if (y_binary) {
    estimate <- exp(est)
    ci <- c(
      exp(est + se * qnorm(alpha / 2)),
      exp(est + se * qnorm(1 - alpha / 2))
    )
  } else {
    estimate <- est
    ci <- c(
      est + se * qnorm(alpha / 2),
      est + se * qnorm(1 - alpha / 2)
    )
  }

  return(list(estimate = estimate, ci = ci))
}
