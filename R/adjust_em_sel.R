adjust_em_sel_val <- function(
    data_observed,
    data_validation) {
  if (!all(data_observed$confounders %in% data_validation$confounders)) {
    stop("All confounders in observed data must be present in validation data.", call. = FALSE)
  }

  if (is.null(data_validation$misclassified_exposure)) {
    stop(
      paste0(
        "This function is adjusting for a misclassified exposure.",
        "\n",
        "Validation data must have a true and misclassified exposure specified."
      ),
      call. = FALSE
    )
  }
  if (is.null(data_validation$selection)) {
    stop(
      paste0(
        "This function is adjusting for selection bias.",
        "\n",
        "Validation data must have a selection indicator column specified."
      ),
      call. = FALSE
    )
  }

  n <- nrow(data_observed$data)

  df <- data.frame(
    Xstar = data_observed$data[, data_observed$exposure],
    Y = data_observed$data[, data_observed$outcome]
  )
  df <- bind_cols(
    df,
    data_observed$data %>%
      select(all_of(data_observed$confounders))
  )

  if (all(df$Y %in% 0:1)) {
    y_binary <- TRUE
  } else {
    y_binary <- FALSE
  }

  df_val <- data.frame(
    X = data_validation$data[, data_validation$true_exposure],
    Y = data_validation$data[, data_validation$true_outcome],
    Xstar = data_validation$data[, data_validation$misclassified_exposure],
    S = data_validation$data[, data_validation$selection]
  )
  df_val <- bind_cols(
    df_val,
    data_validation$data %>%
      select(all_of(data_validation$confounders))
  )

  force_match(
    df$Y,
    df_val$Y,
    "Outcomes from both datasets must both be binary or both be continuous."
  )
  force_binary(
    df$Xstar,
    "Exposure in observed data must be a binary integer."
  )
  force_binary(
    df_val$Xstar,
    "Misclassified exposure in validation data must be a binary integer."
  )
  force_binary(
    df_val$X,
    "True exposure in validation data must be a binary integer."
  )
  force_binary(
    df_val$S,
    "Selection indicator in validation data must be a binary integer."
  )

  x_mod <- glm(X ~ Xstar + Y + . - S,
    family = binomial(link = "logit"),
    data = df_val
  )

  x_mod_coefs <- coef(x_mod)
  x_pred <- x_mod_coefs[1]

  for (i in 2:length(x_mod_coefs)) {
    var_name <- names(x_mod_coefs)[i]
    x_pred <- x_pred + df[[var_name]] * x_mod_coefs[i]
  }

  df$Xpred <- rbinom(n, 1, plogis(x_pred))

  s_mod <- glm(S ~ Xstar + Y + . - X,
    family = binomial(link = "logit"),
    data = df_val
  )

  s_mod_coefs <- coef(s_mod)
  s_pred <- s_mod_coefs[1]

  for (i in 2:length(s_mod_coefs)) {
    var_name <- names(s_mod_coefs)[i]
    s_pred <- s_pred + df[[var_name]] * s_mod_coefs[i]
  }

  df$Spred <- plogis(s_pred)

  if (y_binary) {
    suppressWarnings({
      final <- glm(
        Y ~ Xpred + . - Xstar - Spred,
        family = binomial(link = "logit"),
        weights = (1 / df$Spred),
        data = df
      )
    })
  } else {
    suppressWarnings({
      final <- lm(
        Y ~ Xpred + . - Xstar - Spred,
        weights = (1 / df$Spred),
        data = df
      )
    })
  }

  return(final)
}


adjust_em_sel_coef <- function(
    data_observed,
    x_model_coefs,
    s_model_coefs,
    level = 0.95) {
  data <- data_observed$data
  n <- nrow(data)
  confounders <- data_observed$confounders
  len_c <- length(confounders)
  len_x_coefs <- length(x_model_coefs)
  len_s_coefs <- length(s_model_coefs)

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
  force_len(
    len_s_coefs,
    3 + len_c,
    paste0(
      "Incorrect length of S model coefficients. ",
      "Length should equal 3 + number of confounders."
    )
  )

  if (all(y %in% 0:1)) {
    y_binary <- TRUE
  } else {
    y_binary <- FALSE
  }

  s1_0 <- s_model_coefs[1]
  s1_xstar <- s_model_coefs[2]
  s1_y <- s_model_coefs[3]

  x1_0 <- x_model_coefs[1]
  x1_xstar <- x_model_coefs[2]
  x1_y <- x_model_coefs[3]

  if (is.null(confounders)) {
    df <- data.frame(Xstar = xstar, Y = y)

    x1_pred <- plogis(x1_0 + x1_xstar * xstar + x1_y * y)
    x1_pred <- rep(x1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(
        Xbar = rep(c(1, 0), each = n),
        pS = plogis(s1_0 + s1_xstar * .data$Xstar + s1_y * .data$Y),
        pX = case_when(
          Xbar == 1 ~ x1_pred,
          Xbar == 0 ~ 1 - x1_pred
        )
      )

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
        pS = plogis(
          s1_0 + s1_xstar * .data$Xstar + s1_y * .data$Y +
            s1_c1 * .data$C1
        ),
        pX = case_when(
          Xbar == 1 ~ x1_pred,
          Xbar == 0 ~ 1 - x1_pred
        )
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

    x1_pred <- plogis(
      x1_0 + x1_xstar * xstar +
        x1_y * y + x1_c1 * c1 + x1_c2 * c2
    )
    x1_pred <- rep(x1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(
        Xbar = rep(c(1, 0), each = n),
        pS = plogis(
          s1_0 + s1_xstar * .data$Xstar + s1_y * .data$Y +
            s1_c1 * .data$C1 + s1_c2 * .data$C2
        ),
        pX = case_when(
          Xbar == 1 ~ x1_pred,
          Xbar == 0 ~ 1 - x1_pred
        )
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
        pS = plogis(
          s1_0 + s1_xstar * .data$Xstar + s1_y * .data$Y +
            s1_c1 * .data$C1 + s1_c2 * .data$C2 + s1_c3 * .data$C3
        ),
        pX = case_when(
          Xbar == 1 ~ x1_pred,
          Xbar == 0 ~ 1 - x1_pred
        )
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
    stop(
      "This function is currently not compatible with >3 confounders.",
      call. = FALSE
    )
  }

  return(final)
}


#' Adust for exposure misclassification and selection bias.
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `adjust_emc_sel()` was renamed to `adjust_em_sel()`
#' @keywords internal
#'
#' @export
adjust_emc_sel <- function(
    data_observed,
    x_model_coefs,
    s_model_coefs,
    level = 0.95) {
  lifecycle::deprecate_warn("1.5.3", "adjust_emc_sel()", "adjust_em_sel()")
  adjust_em_sel(
    data_observed,
    x_model_coefs,
    s_model_coefs,
    level
  )
}


#' Adust for exposure misclassification and selection bias.
#'
#' `adjust_em_sel` returns the exposure-outcome odds ratio and confidence
#' interval, adjusted for exposure misclassification and selection bias.
#'
#' Bias adjustment can be performed by inputting either a validation dataset or
#' the necessary bias parameters. Values for the bias parameters
#' can be applied as fixed values or as single draws from a probability
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
#' data, plus data for the true and misclassified exposure,
#' corresponding to the observed exposure in `data_observed`. There should also
#' be a selection indicator representing whether the observation in
#' `data_validation` was selected in `data_observed`.
#' @param x_model_coefs The regression coefficients corresponding to the model:
#' \ifelse{html}{\out{logit(P(X=1)) = &delta;<sub>0</sub> + &delta;<sub>1</sub>X* + &delta;<sub>2</sub>Y + &delta;<sub>2+j</sub>C<sub>j</sub>, }}{\eqn{logit(P(X=1)) = \delta_0 + \delta_1 X^* + \delta_2 Y + \delta{2+j} C_j, }}
#' where *X* represents the binary true exposure, *X** is the
#' binary misclassified exposure, *Y* is the outcome,
#' *C* represents the vector of
#' measured confounders (if any), and *j* corresponds to the number of
#' measured confounders. The number of parameters is therefore 3 + *j*.
#' @param s_model_coefs The regression coefficients corresponding to the model:
#' \ifelse{html}{\out{logit(P(S=1)) = &beta;<sub>0</sub> + &beta;<sub>1</sub>X* + &beta;<sub>2</sub>Y + &beta;<sub>2+j</sub>C<sub>j</sub>, }}{\eqn{logit(P(S=1)) = \beta_0 + \beta_1 X^* + \beta_2 Y + \beta{{2+j}} C_j, }}
#' where *S* represents binary selection, *X** is the
#' binary misclassified exposure,
#' *Y* is the outcome, *C* represents the vector of
#' measured confounders (if any), and *j* corresponds to the number of
#' measured confounders. The number of parameters is therefore 3 + *j*.
#' @param level Value from 0-1 representing the full range of the confidence
#' interval. Default is 0.95.
#'
#' @return A list where the first item is the odds ratio estimate of the
#' effect of the exposure on the outcome and the second item is the
#' confidence interval as the vector: (lower bound, upper bound).
#'
#' @examples
#' df_observed <- data_observed(
#'   data = df_em_sel,
#'   exposure = "Xstar",
#'   outcome = "Y",
#'   confounders = "C1"
#' )
#'
#' # Using validation data -----------------------------------------------------
#' df_validation <- data_validation(
#'   data = df_em_sel_source,
#'   true_exposure = "X",
#'   true_outcome = "Y",
#'   confounders = "C1",
#'   misclassified_exposure = "Xstar",
#'   selection = "S"
#' )
#'
#' adjust_em_sel(
#'   data_observed = df_observed,
#'   data_validation = df_validation
#' )
#'
#' # Using x_model_coefs and s_model_coefs -------------------------------------
#' adjust_em_sel(
#'   data_observed = df_observed,
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

adjust_em_sel <- function(
    data_observed,
    data_validation = NULL,
    x_model_coefs = NULL,
    s_model_coefs = NULL,
    level = 0.95) {
  check_inputs2(
    data_validation,
    list(x_model_coefs, s_model_coefs)
  )

  data <- data_observed$data
  xstar <- data[, data_observed$exposure]
  y <- data[, data_observed$outcome]

  if (all(y %in% 0:1)) {
    y_binary <- TRUE
  } else {
    y_binary <- FALSE
  }

  if (!is.null(data_validation)) {
    final <- adjust_em_sel_val(
      data_observed,
      data_validation
    )
  } else if (!is.null(x_model_coefs)) {
    final <- adjust_em_sel_coef(
      data_observed,
      x_model_coefs,
      s_model_coefs
    )
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
