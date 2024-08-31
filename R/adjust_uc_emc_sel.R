#' Adust for uncontrolled confounding, exposure misclassification, and selection
#' bias.
#'
#' \code{adjust_uc_emc_sel} returns the exposure-outcome odds ratio and
#' confidence interval, adjusted for uncontrolled confounding, exposure
#' misclassificaiton, and selection bias.
#'
#' This function uses two separate logistic regression models to predict the
#' uncontrolled confounder (\emph{U}) and exposure (\emph{X}).
#' If a single bias model for jointly modeling \emph{X} and \emph{U} is desired
#' use \code{adjust_multinom_uc_emc_sel}.
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
#' @param u_model_coefs The regression coefficients corresponding to the model:
#' \ifelse{html}{\out{logit(P(U=1)) = &alpha;<sub>0</sub> + &alpha;<sub>1</sub>X + &alpha;<sub>2</sub>Y, }}{\eqn{logit(P(U=1)) = \alpha_0 + \alpha_1 X + \alpha_2 Y, }}
#' where \emph{U} is the binary unmeasured confounder, \emph{X} is the
#' binary true exposure, and \emph{Y} is the outcome.
#' The number of parameters therefore equals 3.
#' @param x_model_coefs The regression coefficients corresponding to the model:
#' \ifelse{html}{\out{logit(P(X=1)) = &delta;<sub>0</sub> + &delta;<sub>1</sub>X* + &delta;<sub>2</sub>Y + &delta;<sub>2+j</sub>C<sub>j</sub>, }}{\eqn{logit(P(X=1)) = \delta_0 + \delta_1 X^* + \delta_2 Y + \delta_{2+j} C_j, }}
#' where \emph{X} represents binary true exposure, \emph{X*} is the
#' binary misclassified exposure, \emph{Y} is the outcome, \emph{C}
#' represents the vector of measured confounders (if any), and
#' \emph{j} corresponds to the number of measured confounders.
#' The number of parameters therefore equals 3 + \emph{j}.
#' @param s_model_coefs The regression coefficients corresponding to the model:
#' \ifelse{html}{\out{logit(P(S=1)) = &beta;<sub>0</sub> + &beta;<sub>1</sub>X* + &beta;<sub>2</sub>Y + &beta;<sub>2+j</sub>C<sub>2+j</sub>, }}{\eqn{logit(P(S=1)) = \beta_0 + \beta_1 X^* + \beta_2 Y + \beta_{2+j} C_j, }}
#' where \emph{S} represents binary selection, \emph{X*} is the
#' binary misclassified exposure, \emph{Y} is the outcome,
#' \emph{C} represents the vector of measured confounders (if any),
#' and \emph{j} corresponds to the number of measured confounders.
#' The number of parameters therefore equals 3 + \emph{j}.
#' @return A list where the first item is the odds ratio estimate of the
#' effect of the exposure on the outcome and the second item is the
#' confidence interval as the vector: (lower bound, upper bound).
#'
#' @examples
#' adjust_uc_emc_sel(
#'   df_uc_emc_sel,
#'   exposure = "Xstar",
#'   outcome = "Y",
#'   confounders = c("C1", "C2", "C3"),
#'   u_model_coefs = c(-0.32, 0.59, 0.69),
#'   x_model_coefs = c(-2.44, 1.62, 0.72, 0.32, -0.15, 0.85),
#'   s_model_coefs = c(0.00, 0.26, 0.78, 0.03, -0.02, 0.10)
#' )
#'
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom stats binomial
#' @importFrom stats glm
#' @importFrom stats qnorm
#' @importFrom stats rbinom
#' @importFrom stats plogis
#' @importFrom rlang .data
#'
#' @export

adjust_uc_emc_sel <- function(
  data,
  exposure,
  outcome,
  confounders = NULL,
  u_model_coefs,
  x_model_coefs,
  s_model_coefs,
  level = 0.95
) {

  n <- nrow(data)
  len_c <- length(confounders)

  len_u_coefs <- length(u_model_coefs)
  len_x_coefs <- length(x_model_coefs)
  len_s_coefs <- length(s_model_coefs)

  xstar <- data[, exposure]
  y <- data[, outcome]

  if (sum(xstar %in% c(0, 1)) != n) {
    stop("Exposure must be a binary integer.")
  }
  if (sum(y %in% c(0, 1)) == n) {
    y_binary <- TRUE
  } else {
    y_binary <- FALSE
  }

  if (len_u_coefs != 3) {
    stop("Incorrect length of U model coefficients. Length should equal 3.")
  }
  if (len_x_coefs != 3 + len_c) {
    stop(
      paste0(
        "Incorrect length of X model coefficients. ",
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

  u1_0     <- u_model_coefs[1]
  u1_x     <- u_model_coefs[2]
  u1_y     <- u_model_coefs[3]

  x1_0     <- x_model_coefs[1]
  x1_xstar <- x_model_coefs[2]
  x1_y     <- x_model_coefs[3]

  s1_0     <- s_model_coefs[1]
  s1_xstar <- s_model_coefs[2]
  s1_y     <- s_model_coefs[3]

  if (is.null(confounders)) {

    df <- data.frame(Xstar = xstar, Y = y)

    df <- df %>%
      mutate(
        Xpred = rbinom(
          n, 1,
          plogis(x1_0 + x1_xstar * .data$Xstar + x1_y * .data$Y)
        ),
        Upred = rbinom(
          n, 1,
          plogis(u1_0 + u1_x * .data$Xpred + u1_y * .data$Y)
        ),
        pS = plogis(s1_0 + s1_xstar * .data$Xstar + s1_y * .data$Y)
      )

    if (y_binary) {
      suppressWarnings({
        final <- glm(
          Y ~ Xpred + Upred,
          family = binomial(link = "logit"),
          weights = (1 / df$pS),
          data = df
        )
      })
    } else {
      final <- lm(
        Y ~ Xpred + Upred,
        weights = (1 / df$pS),
        data = df
      )
    }

  } else if (len_c == 1) {

    c1 <- data[, confounders]
    df <- data.frame(Xstar = xstar, Y = y, C1 = c1)

    x1_c1 <- x_model_coefs[4]
    s1_c1 <- s_model_coefs[4]

    df <- df %>%
      mutate(
        Xpred = rbinom(n, 1, plogis(x1_0 + x1_xstar * .data$Xstar +
                                      x1_y * .data$Y + x1_c1 * .data$C1)),
        Upred = rbinom(n, 1, plogis(u1_0 + u1_x * .data$Xpred +
                                      u1_y * .data$Y)),
        pS = plogis(s1_0 + s1_xstar * .data$Xstar + s1_y * .data$Y +
                      s1_c1 * .data$C1)
      )

    if (y_binary) {
      suppressWarnings({
        final <- glm(
          Y ~ Xpred + C1 + Upred,
          family = binomial(link = "logit"),
          weights = (1 / df$pS),
          data = df
        )
      })
    } else {
      final <- lm(
        Y ~ Xpred + C1 + Upred,
        weights = (1 / df$pS),
        data = df
      )
    }

  } else if (len_c == 2) {

    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]

    df <- data.frame(Xstar = xstar, Y = y, C1 = c1, C2 = c2)

    x1_c1 <- x_model_coefs[4]
    x1_c2 <- x_model_coefs[5]

    s1_c1 <- s_model_coefs[4]
    s1_c2 <- s_model_coefs[5]

    df <- df %>%
      mutate(
        Xpred = rbinom(
          n, 1,
          plogis(
            x1_0 + x1_xstar * .data$Xstar + x1_y * .data$Y +
              x1_c1 * .data$C1 + x1_c2 * .data$C2
          )
        ),
        Upred = rbinom(n, 1, plogis(u1_0 + u1_x * .data$Xpred +
                                      u1_y * .data$Y)),
        pS = plogis(s1_0 + s1_xstar * .data$Xstar + s1_y * .data$Y +
                      s1_c1 * .data$C1 + s1_c2 * .data$C2)
      )

    if (y_binary) {
      suppressWarnings({
        final <- glm(
          Y ~ Xpred + C1 + C2 + Upred,
          family = binomial(link = "logit"),
          weights = (1 / df$pS),
          data = df
        )
      })
    } else {
      final <- lm(
        Y ~ Xpred + C1 + C2 + Upred,
        weights = (1 / df$pS),
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

    s1_c1 <- s_model_coefs[4]
    s1_c2 <- s_model_coefs[5]
    s1_c3 <- s_model_coefs[6]

    df <- df %>%
      mutate(
        Xpred = rbinom(
          n, 1,
          plogis(x1_0 + x1_xstar * .data$Xstar + x1_y * .data$Y +
                   x1_c1 * .data$C1 + x1_c2 * .data$C2 + x1_c3 * .data$C3)
        ),
        Upred = rbinom(n, 1, plogis(u1_0 + u1_x * .data$Xpred +
                                      u1_y * .data$Y)),
        pS = plogis(s1_0 + s1_xstar * .data$Xstar + s1_y * .data$Y +
                      s1_c1 * .data$C1 + s1_c2 * .data$C2 + s1_c3 * .data$C3)
      )

    if (y_binary) {
      suppressWarnings({
        final <- glm(
          Y ~ Xpred + C1 + C2 + C3 + Upred,
          family = binomial(link = "logit"),
          weights = (1 / df$pS),
          data = df
        )
      })
    } else {
      final <- lm(
        Y ~ Xpred + C1 + C2 + C3 + Upred,
        weights = (1 / df$pS),
        data = df
      )
    }

  }

  if (len_c > 3) {

    stop("This function is currently not compatible with >3 confounders.")

  }

  est <- summary(final)$coef[2, 1]
  se <- summary(final)$coef[2, 2]
  alpha <- 1 - level

  if (y_binary) {
    estimate <- exp(est)
    ci <- c(exp(est + se * qnorm(alpha / 2)),
            exp(est + se * qnorm(1 - alpha / 2)))
  } else {
    estimate <- est
    ci <- c(est + se * qnorm(alpha / 2),
            est + se * qnorm(1 - alpha / 2))
  }

  return(list(estimate = estimate, ci = ci))

}