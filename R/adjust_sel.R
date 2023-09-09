#' Adust for selection bias.
#'
#' \code{adjust_sel} returns the exposure-outcome odds ratio and confidence
#' interval, adjusted for selection bias.
#'
#' Details
#'
#' @inheritParams adjust_emc_sel
#' @param s_model_coefs The regression coefficients corresponding to the model:
#'  \ifelse{html}{\out{logit(P(S=1)) = &beta;<sub>0</sub> +
#'  &beta;<sub>1</sub>X + &beta;<sub>2</sub>Y, } where S represents
#'  (binary) selection, X is the (binary) exposure, Y is the (binary)
#'  outcome. The number of parameters is therefore 3.}{\eqn{logit(P(S=1)) =}}
#'
#' @examples
#' adjust_sel(
#'  evans,
#'  exposure = "SMK",
#'  outcome = "CHD",
#'  confounders = "AGE",
#'  s_model_coefs = c(0.05, 0.20, 0.20)
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

adjust_sel <- function(
  data,
  exposure,
  outcome,
  confounders = NULL,
  s_model_coefs,
  level = 0.95
) {

  n <- nrow(data)
  len_c <- length(confounders)
  len_s_coefs <- length(s_model_coefs)

  x <- data[, exposure]
  y <- data[, outcome]

  if (sum(x %in% c(0, 1)) != n) {
    stop("Exposure must be a binary integer.")
  }
  if (sum(y %in% c(0, 1)) != n) {
    stop("Outcome must be a binary integer.")
  }
  if (len_s_coefs != 3) {
    stop(
      paste0(
        "Incorrect length of S model coefficients. ",
        "Length should equal 3."
      )
    )
  }

  s1_0 <- s_model_coefs[1]
  s1_x <- s_model_coefs[2]
  s1_y <- s_model_coefs[3]

  if (is.null(confounders)) {

    df <- data.frame(X = x, Y = y)
    df$pS <- plogis(s1_0 + s1_x * df$X + s1_y * df$Y)

    final <- glm(
      Y ~ X,
      family = binomial(link = "logit"),
      weights = df$pS,
      data = df
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
    df$pS <- plogis(s1_0 + s1_x * df$X + s1_y * df$Y)

    final <- glm(
      Y ~ X + C1,
      family = binomial(link = "logit"),
      weights = df$pS,
      data = df
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
    df$pS <- plogis(s1_0 + s1_x * df$X + s1_y * df$Y)

    final <- glm(
      Y ~ X + C1 + C2,
      family = binomial(link = "logit"),
      weights = df$pS,
      data = df
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
    df$pS <- plogis(s1_0 + s1_x * df$X + s1_y * df$Y)

    final <- glm(
      Y ~ X + C1 + C2 + C3,
      family = binomial(link = "logit"),
      weights = df$pS,
      data = df
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