adjust_om_val <- function(
    data_observed,
    data_validation) {
  if (!all(data_observed$confounders %in% data_validation$confounders)) {
    stop("All confounders in observed data must be present in validation data.")
  }

  if (is.null(data_validation$misclassified_outcome)) {
    stop(
      paste0(
        "This function is adjusting for a misclassified outcome.",
        "\n",
        "Validation data must have a true and misclassified outcome specified."
      )
    )
  }

  n <- nrow(data_observed$data)

  df <- data.frame(
    X = data_observed$data[, data_observed$exposure],
    Ystar = data_observed$data[, data_observed$outcome]
  )
  df <- bind_cols(df, data_observed$data[, data_observed$confounders])

  df_val <- data.frame(
    X = data_validation$data[, data_validation$true_exposure],
    Y = data_validation$data[, data_validation$true_outcome],
    Ystar = data_validation$data[, data_validation$misclassified_outcome]
  )
  df_val <- bind_cols(
    df_val,
    data_validation$data[, data_validation$confounders]
  )

  if (all(df$X %in% 0:1)) {
    if (!all(df_val$X %in% 0:1)) {
      stop("Exposures from both datasets must both be binary or both be continuous.")
    }
  }
  if (all(df$Y %in% 0:1)) {
    if (!all(df_val$Y %in% 0:1) || !all(df_val$Ystar %in% 0:1)) {
      stop("Outcomes from both datasets must both be binary or both be continuous.")
    }
  }
  if (!all(df$Ystar %in% 0:1)) {
    stop("Misclassified outcome in observed data must be a binary integer.")
  }
  if (!all(df_val$Ystar %in% 0:1)) {
    stop("Misclassified outcome in validation data must be a binary integer.")
  }
  if (!all(df_val$Y %in% 0:1)) {
    stop("True outcome in validation data must be a binary integer.")
  }

  y_mod <- glm(Y ~ X + Ystar + .,
    family = binomial(link = "logit"),
    data = df_val
  )

  y_mod_coefs <- coef(y_mod)
  y_pred <- y_mod_coefs[1]

  for (i in 2:length(y_mod_coefs)) {
    var_name <- names(y_mod_coefs)[i]
    y_pred <- y_pred + df[[var_name]] * y_mod_coefs[i]
  }

  df$Ypred <- rbinom(n, 1, plogis(y_pred))

  final <- glm(
    Ypred ~ X + . - Ystar,
    family = binomial(link = "logit"),
    data = df
  )

  return(final)
}


adjust_om_coef <- function(
    data_observed,
    y_model_coefs) {
  data <- data_observed$data
  n <- nrow(data)
  confounders <- data_observed$confounders
  len_c <- length(confounders)
  len_y_coefs <- length(y_model_coefs)

  x <- data[, data_observed$exposure]
  ystar <- data[, data_observed$outcome]

  force_binary(ystar, "Outcome must be a binary integer.")
  force_len(
    len_y_coefs,
    3 + len_c,
    paste0(
      "Incorrect length of Y model coefficients. ",
      "Length should equal 3 + number of confounders."
    )
  )

  y1_0 <- y_model_coefs[1]
  y1_x <- y_model_coefs[2]
  y1_ystar <- y_model_coefs[3]

  if (is.null(confounders)) {
    df <- data.frame(X = x, Ystar = ystar)
    df$Ypred <- rbinom(n, 1, plogis(y1_0 + y1_x * df$X + y1_ystar * df$Ystar))

    final <- glm(
      Ypred ~ X,
      family = binomial(link = "logit"),
      data = df
    )
  } else if (len_c == 1) {
    c1 <- data[, confounders]
    df <- data.frame(X = x, Ystar = ystar, C1 = c1)

    y1_c1 <- y_model_coefs[4]

    df$Ypred <- rbinom(
      n, 1, plogis(
        y1_0 + y1_x * df$X + y1_ystar * df$Ystar +
          y1_c1 * df$C1
      )
    )

    final <- glm(
      Ypred ~ X + C1,
      family = binomial(link = "logit"),
      data = df
    )
  } else if (len_c == 2) {
    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]

    df <- data.frame(X = x, Ystar = ystar, C1 = c1, C2 = c2)

    y1_c1 <- y_model_coefs[4]
    y1_c2 <- y_model_coefs[5]

    df$Ypred <- rbinom(
      n, 1, plogis(
        y1_0 + y1_x * df$X + y1_ystar * df$Ystar +
          y1_c1 * df$C1 + y1_c2 * df$C2
      )
    )

    final <- glm(
      Ypred ~ X + C1 + C2,
      family = binomial(link = "logit"),
      data = df
    )
  } else if (len_c == 3) {
    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]
    c3 <- data[, confounders[3]]

    df <- data.frame(X = x, Ystar = ystar, C1 = c1, C2 = c2, C3 = c3)

    y1_c1 <- y_model_coefs[4]
    y1_c2 <- y_model_coefs[5]
    y1_c3 <- y_model_coefs[6]

    df$Ypred <- rbinom(
      n, 1, plogis(
        y1_0 + y1_x * df$X + y1_ystar * df$Ystar +
          y1_c1 * df$C1 + y1_c2 * df$C2 + y1_c3 * df$C3
      )
    )

    final <- glm(
      Ypred ~ X + C1 + C2 + C3,
      family = binomial(link = "logit"),
      data = df
    )
  } else if (len_c > 3) {
    stop("This function is currently not compatible with >3 confounders.")
  }

  return(final)
}


#' Adust for outcome misclassification.
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `adjust_omc()` was renamed to `adjust_om()`
#' @keywords internal
#'
#' @export
adjust_omc <- function(
    data_observed,
    y_model_coefs,
    level = 0.95) {
  lifecycle::deprecate_warn("1.5.3", "adjust_omc()", "adjust_om()")
  adjust_om(data_observed, y_model_coefs, level)
}


#' Adust for outcome misclassification.
#'
#' `adjust_om` returns the exposure-outcome odds ratio and confidence
#' interval, adjusted for outcome misclassificaiton.
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
#' @param data_validation Object of class `data_validation` corresponding to
#' the validation data used to adjust for bias in the observed data. Here, the
#' validation data should have data for the same variables as in the observed
#' data, plus data for the true and misclassified outcome corresponding to the
#' observed outcome in `data_observed`.
#' @param y_model_coefs The regression coefficients corresponding to the model:
#' \ifelse{html}{\out{logit(P(Y=1)) = &delta;<sub>0</sub> + &delta;<sub>1</sub>X + &delta;<sub>2</sub>Y* + &delta;<sub>2+j</sub>C<sub>j</sub>, }}{\eqn{logit(P(Y=1)) = \_delta_0 + \_delta_1 X + \_delta_2 Y^* + \_delta_{2+j} C_j, }}
#' where *Y* represents the binary true outcome, *X* is the exposure,
#' *Y** is the binary misclassified outcome,
#' *C* represents the vector of measured confounders (if any),
#' and *j* corresponds to the number of measured confounders. The number
#' of parameters is therefore 3 + *j*.
#' @param level Value from 0-1 representing the full range of the confidence
#' interval. Default is 0.95.
#'
#' @return A list where the first item is the odds ratio estimate of the
#' effect of the exposure on the outcome and the second item is the
#' confidence interval as the vector: (lower bound, upper bound).
#'
#' @examples
#' df_observed <- data_observed(
#'   data = df_om,
#'   exposure = "X",
#'   outcome = "Ystar",
#'   confounders = "C1"
#' )
#' # Using validation data -----------------------------------------------------
#' df_validation <- data_validation(
#'   data = df_om_source,
#'   true_exposure = "X",
#'   true_outcome = "Y",
#'   confounders = "C1",
#'   misclassified_outcome = "Ystar"
#' )
#'
#' adjust_om(
#'   data_observed = df_observed,
#'   data_validation = df_validation
#' )
#'
#' # Using y_model_coefs -------------------------------------------------------
#' adjust_om(
#'   data_observed = df_observed,
#'   y_model_coefs = c(-3.1, 0.6, 1.6, 0.4)
#' )
#'
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom stats binomial
#' @importFrom stats glm
#' @importFrom stats qnorm
#' @importFrom stats rbinom
#' @importFrom stats plogis
#' @importFrom stats coef
#' @importFrom rlang .data
#'
#' @export

adjust_om <- function(
    data_observed,
    data_validation = NULL,
    y_model_coefs = NULL,
    level = 0.95) {
  if (
    (!is.null(data_validation) && !is.null(y_model_coefs)) ||
      (is.null(data_validation) && is.null(y_model_coefs))
  ) {
    stop("One of data_validation or y_model_coefs must be non-null.")
  }
  data <- data_observed$data
  n <- nrow(data)
  confounders <- data_observed$confounders
  len_c <- length(confounders)

  if (!is.null(data_validation)) {
    final <- adjust_om_val(
      data_observed,
      data_validation
    )
  } else if (!is.null(y_model_coefs)) {
    final <- adjust_om_coef(
      data_observed,
      y_model_coefs
    )
  }

  est <- summary(final)$coef[2, 1]
  se <- summary(final)$coef[2, 2]
  alpha <- 1 - level

  estimate <- exp(est)
  ci <- c(
    exp(est + se * qnorm(alpha / 2)),
    exp(est + se * qnorm(1 - alpha / 2))
  )

  return(list(estimate = estimate, ci = ci))
}
