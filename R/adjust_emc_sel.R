#' Adust for exposure misclassification and selection bias.
#'
#' \code{adjust_emc_sel} returns the exposure-outcome odds ratio and confidence
#' interval, adjusted for exposure misclassification and selection bias.
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
#' @param data Dataframe for analysis.
#' @param exposure String name of the exposure variable.
#' @param outcome String name of the outcome variable.
#' @param confounders String name(s) of the confounder(s).
#'  A maximum of three confounders is allowed.
#' @param x_model_coefs The regression coefficients corresponding to the model:
#'  \ifelse{html}{\out{logit(P(X=1)) = &delta;<sub>0</sub> + &delta;<sub>1</sub>X* + &delta;<sub>2</sub>Y + &delta;<sub>2+j</sub>C<sub>j</sub>, }}{\eqn{logit(P(X=1)) = \delta_0 + \delta_1 X^* + \delta_2 Y + \delta{2+j} C_j, }}
#'  where \emph{X} represents the binary true exposure, \emph{X*} is the
#'  binary misclassified exposure, \emph{Y} is the outcome,
#'  \emph{C} represents the vector of
#'  measured confounders (if any), and \emph{j} corresponds to the number of
#'  measured confounders. The number of parameters is therefore 3 + \emph{j}.
#' @param s_model_coefs The regression coefficients corresponding to the model:
#'  \ifelse{html}{\out{logit(P(S=1)) = &beta;<sub>0</sub> + &beta;<sub>1</sub>X* + &beta;<sub>2</sub>Y + &beta;<sub>2+j</sub>C<sub>j</sub>, }}{\eqn{logit(P(S=1)) = \beta_0 + \beta_1 X^* + \beta_2 Y + \beta{{2+j}} C_j, }}
#'  where \emph{S} represents binary selection, \emph{X*} is the
#'  binary misclassified exposure,
#'  \emph{Y} is the outcome, \emph{C} represents the vector of
#'  measured confounders (if any), and \emph{j} corresponds to the number of
#'  measured confounders. The number of parameters is therefore 3 + \emph{j}.
#' @param level Value from 0-1 representing the full range of the confidence
#'  interval. Default is 0.95.
#' @return A list where the first item is the odds ratio estimate of the
#'  effect of the exposure on the outcome and the second item is the
#'  confidence interval as the vector: (lower bound, upper bound).
#'
#' @examples
#' adjust_emc_sel(
#'   df_emc_sel,
#'   exposure = "Xstar",
#'   outcome = "Y",
#'   confounders = "C1",
#'   x_model_coefs = c(-2.78, 1.62, 0.58, 0.34),
#'   s_model_coefs = c(0.04, 0.18, 0.92, 0.05)
#' )
#'
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom stats binomial
#' @importFrom stats glm
#' @importFrom stats lm
#' @importFrom stats qnorm
#' @importFrom stats plogis
#' @importFrom rlang .data
#'
#' @export

adjust_emc_sel <- function(
  data,
  exposure,
  outcome,
  confounders = NULL,
  x_model_coefs,
  s_model_coefs,
  level = 0.95
) {

  n <- nrow(data)
  len_c <- length(confounders)
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

  s1_0     <- s_model_coefs[1]
  s1_xstar <- s_model_coefs[2]
  s1_y     <- s_model_coefs[3]

  x1_0     <- x_model_coefs[1]
  x1_xstar <- x_model_coefs[2]
  x1_y     <- x_model_coefs[3]

  if (is.null(confounders)) {

    df <- data.frame(Xstar = xstar, Y = y)

    x1_pred <- plogis(x1_0 + x1_xstar * xstar + x1_y * y)
    x1_pred <- rep(x1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(Xbar = rep(c(1, 0), each = n),
             pS = plogis(s1_0 + s1_xstar * .data$Xstar + s1_y * .data$Y),
             pX = case_when(Xbar == 1 ~ x1_pred,
                            Xbar == 0 ~ 1 - x1_pred))

    if (y_binary) {
      suppressWarnings({
        final <- glm(
          Y ~ Xbar,
          family = binomial(link = "logit"),
          weights = (combined$pX / combined$pS),
          data = combined
        )
      })
    } else {
      suppressWarnings({
        final <- lm(
          Y ~ Xbar,
          weights = (combined$pX / combined$pS),
          data = combined
        )
      })
    }

  } else if (len_c == 1) {

    c1 <- data[, confounders]
    df <- data.frame(Xstar = xstar, Y = y, C1 = c1)
    x1_c1 <- x_model_coefs[4]
    s1_c1 <- s_model_coefs[4]

    x1_pred <- plogis(x1_0 + x1_xstar * xstar + x1_y * y + x1_c1 * c1)
    x1_pred <- rep(x1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(
        Xbar = rep(c(1, 0), each = n),
        pS = plogis(s1_0 + s1_xstar * .data$Xstar + s1_y * .data$Y +
                      s1_c1 * .data$C1),
        pX = case_when(Xbar == 1 ~ x1_pred,
                       Xbar == 0 ~ 1 - x1_pred)
      )

    if (y_binary) {
      suppressWarnings({
        final <- glm(
          Y ~ Xbar + C1,
          family = binomial(link = "logit"),
          weights = (combined$pX / combined$pS),
          data = combined
        )
      })
    } else {
      suppressWarnings({
        final <- lm(
          Y ~ Xbar + C1,
          weights = (combined$pX / combined$pS),
          data = combined
        )
      })
    }

  } else if (len_c == 2) {

    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]

    df <- data.frame(Xstar = xstar, Y = y, C1 = c1, C2 = c2)

    s1_c1 <- s_model_coefs[4]
    s1_c2 <- s_model_coefs[5]

    x1_c1 <- x_model_coefs[4]
    x1_c2 <- x_model_coefs[5]

    x1_pred <- plogis(x1_0 + x1_xstar * xstar +
                        x1_y * y + x1_c1 * c1 + x1_c2 * c2)
    x1_pred <- rep(x1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(
        Xbar = rep(c(1, 0), each = n),
        pS = plogis(s1_0 + s1_xstar * .data$Xstar + s1_y * .data$Y +
                      s1_c1 * .data$C1 + s1_c2 * .data$C2),
        pX = case_when(Xbar == 1 ~ x1_pred,
                       Xbar == 0 ~ 1 - x1_pred)
      )

    if (y_binary) {
      suppressWarnings({
        final <- glm(
          Y ~ Xbar + C1 + C2,
          family = binomial(link = "logit"),
          weights = (combined$pX / combined$pS),
          data = combined
        )
      })
    } else {
      suppressWarnings({
        final <- lm(
          Y ~ Xbar + C1 + C2,
          weights = (combined$pX / combined$pS),
          data = combined
        )
      })
    }

  } else if (len_c == 3) {

    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]
    c3 <- data[, confounders[3]]

    df <- data.frame(Xstar = xstar, Y = y, C1 = c1, C2 = c2, C3 = c3)

    s1_c1 <- s_model_coefs[4]
    s1_c2 <- s_model_coefs[5]
    s1_c3 <- s_model_coefs[6]

    x1_c1 <- x_model_coefs[4]
    x1_c2 <- x_model_coefs[5]
    x1_c3 <- x_model_coefs[6]

    x1_pred <- plogis(
      x1_0 + x1_xstar * xstar + x1_y * y + x1_c1 * c1 + x1_c2 * c2 + x1_c3 * c3
    )
    x1_pred <- rep(x1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(
        Xbar = rep(c(1, 0), each = n),
        pS = plogis(s1_0 + s1_xstar * .data$Xstar + s1_y * .data$Y +
                      s1_c1 * .data$C1 + s1_c2 * .data$C2 + s1_c3 * .data$C3),
        pX = case_when(Xbar == 1 ~ x1_pred,
                       Xbar == 0 ~ 1 - x1_pred)
      )

    if (y_binary) {
      suppressWarnings({
        final <- glm(
          Y ~ Xbar + C1 + C2 + C3,
          family = binomial(link = "logit"),
          weights = (combined$pX / combined$pS),
          data = combined
        )
      })
    } else {
      suppressWarnings({
        final <- lm(
          Y ~ Xbar + C1 + C2 + C3,
          weights = (combined$pX / combined$pS),
          data = combined
        )
      })
    }

  } else if (len_c > 3) {

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