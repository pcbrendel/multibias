#' Adust for uncontrolled confounding.
#'
#' \code{adjust_uc} returns the exposure-outcome odds ratio and confidence
#' interval, adjusted for uncontrolled confounding from a binary confounder.
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
#' @inheritParams adjust_em_sel
#' @param u_model_coefs The regression coefficients corresponding to the model:
#' \ifelse{html}{\out{logit(P(U=1)) = &alpha;<sub>0</sub> + &alpha;<sub>1</sub>X + &alpha;<sub>2</sub>Y + &alpha;<sub>2+j</sub>C<sub>j</sub>, }}{\eqn{logit(P(U=1)) = \alpha_0 + \alpha_1 X + \alpha_2 Y + \alpha_{2+j} C_j, }}
#' where \emph{U} is the binary unmeasured confounder, \emph{X} is the
#' exposure, \emph{Y} is the outcome, \emph{C} represents the vector of
#' measured confounders (if any),
#' and \emph{j} corresponds to the number of measured confounders.
#' The number of parameters therefore equals 3 + \emph{j}.
#' @return A list where the first item is the odds ratio estimate of the
#' effect of the exposure on the outcome and the second item is the
#' confidence interval as the vector: (lower bound, upper bound).
#'
#' @examples
#' adjust_uc(
#'   evans,
#'   exposure = "SMK",
#'   outcome = "CHD",
#'   confounders = "HPT",
#'   u_model_coefs = c(qlogis(0.25), log(0.5), log(2.5), log(2)),
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

adjust_uc <- function(
    data,
    exposure,
    outcome,
    confounders = NULL,
    u_model_coefs,
    level = 0.95) {
  n <- nrow(data)
  len_c <- length(confounders)
  len_u_coefs <- length(u_model_coefs)

  x <- data[, exposure]
  y <- data[, outcome]

  if (sum(y %in% c(0, 1)) == n) {
    y_binary <- TRUE
  } else {
    y_binary <- FALSE
  }

  if (len_u_coefs != 3 + len_c) {
    stop(
      paste0(
        "Incorrect length of U model coefficients. ",
        "Length should equal 3 + number of confounders."
      )
    )
  }

  u1_0 <- u_model_coefs[1]
  u1_x <- u_model_coefs[2]
  u1_y <- u_model_coefs[3]

  if (is.null(confounders)) {
    df <- data.frame(X = x, Y = y)
    df$Upred <- rbinom(n, 1, plogis(u1_0 + u1_x * df$X + u1_y * df$Y))

    if (y_binary) {
      final <- glm(
        Y ~ X + Upred,
        family = binomial(link = "logit"),
        data = df
      )
    } else {
      final <- lm(
        Y ~ X + Upred,
        data = df
      )
    }
  } else if (len_c == 1) {
    c1 <- data[, confounders]
    df <- data.frame(X = x, Y = y, C1 = c1)

    u1_c1 <- u_model_coefs[4]

    df$Upred <- rbinom(
      n, 1, plogis(u1_0 + u1_x * df$X + u1_y * df$Y + u1_c1 * df$C1)
    )

    if (y_binary) {
      final <- glm(
        Y ~ X + C1 + Upred,
        family = binomial(link = "logit"),
        data = df
      )
    } else {
      final <- lm(
        Y ~ X + C1 + Upred,
        data = df
      )
    }
  } else if (len_c == 2) {
    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]

    df <- data.frame(X = x, Y = y, C1 = c1, C2 = c2)

    u1_c1 <- u_model_coefs[4]
    u1_c2 <- u_model_coefs[5]

    df$Upred <- rbinom(
      n, 1, plogis(
        u1_0 + u1_x * df$X + u1_y * df$Y + u1_c1 * df$C1 + u1_c2 * df$C2
      )
    )

    if (y_binary) {
      final <- glm(
        Y ~ X + C1 + C2 + Upred,
        family = binomial(link = "logit"),
        data = df
      )
    } else {
      final <- lm(
        Y ~ X + C1 + C2 + Upred,
        data = df
      )
    }
  } else if (len_c == 3) {
    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]
    c3 <- data[, confounders[3]]

    df <- data.frame(X = x, Y = y, C1 = c1, C2 = c2, C3 = c3)

    u1_c1 <- u_model_coefs[4]
    u1_c2 <- u_model_coefs[5]
    u1_c3 <- u_model_coefs[6]

    df$Upred <- rbinom(
      n, 1,
      plogis(
        u1_0 + u1_x * df$X + u1_y * df$Y +
          u1_c1 * df$C1 + u1_c2 * df$C2 + u1_c3 * df$C3
      )
    )

    if (y_binary) {
      final <- glm(
        Y ~ X + C1 + C2 + C3 + Upred,
        family = binomial(link = "logit"),
        data = df
      )
    } else {
      final <- lm(
        Y ~ X + C1 + C2 + C3 + Upred,
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
