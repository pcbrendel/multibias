adjust_sel_val <- function(
    data_observed,
    data_validation) {
  if (!all(data_observed$confounders %in% data_validation$confounders)) {
    stop(
      "All confounders in observed data must be present in validation data.",
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
    X = data_observed$data[, data_observed$exposure],
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
    S = data_validation$data[, data_validation$selection]
  )

  force_match(
    df$X,
    df_val$X,
    "Exposures from both datasets must both be binary or both be continuous."
  )
  force_match(
    df$Y,
    df_val$Y,
    "Outcomes from both datasets must both be binary or both be continuous."
  )
  force_binary(
    df_val$S,
    "Selection indicator in validation data must be a binary integer."
  )

  s_mod <- glm(S ~ X + Y,
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
        Y ~ X + . - Spred,
        family = binomial(link = "logit"),
        weights = (1 / df$Spred),
        data = df
      )
    })
  } else {
    suppressWarnings({
      final <- lm(
        Y ~ X + . - Spred,
        weights = (1 / df$Spred),
        data = df
      )
    })
  }

  return(final)
}


adjust_sel_coef <- function(
    data_observed,
    s_model_coefs) {
  data <- data_observed$data
  n <- nrow(data)
  confounders <- data_observed$confounders
  len_c <- length(confounders)
  len_s_coefs <- length(s_model_coefs)

  x <- data[, data_observed$exposure]
  y <- data[, data_observed$outcome]

  force_len(
    len_s_coefs,
    3,
    paste0(
      "Incorrect length of S model coefficients. ",
      "Length should equal 3."
    )
  )

  if (all(y %in% 0:1)) {
    y_binary <- TRUE
  } else {
    y_binary <- FALSE
  }

  s1_0 <- s_model_coefs[1]
  s1_x <- s_model_coefs[2]
  s1_y <- s_model_coefs[3]

  if (is.null(confounders)) {
    df <- data.frame(X = x, Y = y)
    df$pS <- plogis(s1_0 + s1_x * df$X + s1_y * df$Y)

    if (y_binary) {
      suppressWarnings({
        final <- glm(
          Y ~ X,
          family = binomial(link = "logit"),
          weights = (1 / df$pS),
          data = df
        )
      })
    } else {
      suppressWarnings({
        final <- lm(
          Y ~ X,
          weights = (1 / df$pS),
          data = df
        )
      })
    }
  } else if (len_c == 1) {
    c1 <- data[, confounders]
    df <- data.frame(X = x, Y = y, C1 = c1)
    df$pS <- plogis(s1_0 + s1_x * df$X + s1_y * df$Y)

    if (y_binary) {
      suppressWarnings({
        final <- glm(
          Y ~ X + C1,
          family = binomial(link = "logit"),
          weights = (1 / df$pS),
          data = df
        )
      })
    } else {
      suppressWarnings({
        final <- lm(
          Y ~ X + C1,
          weights = (1 / df$pS),
          data = df
        )
      })
    }
  } else if (len_c == 2) {
    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]

    df <- data.frame(X = x, Y = y, C1 = c1, C2 = c2)
    df$pS <- plogis(s1_0 + s1_x * df$X + s1_y * df$Y)

    if (y_binary) {
      suppressWarnings({
        final <- glm(
          Y ~ X + C1 + C2,
          family = binomial(link = "logit"),
          weights = (1 / df$pS),
          data = df
        )
      })
    } else {
      suppressWarnings({
        final <- lm(
          Y ~ X + C1 + C2,
          weights = (1 / df$pS),
          data = df
        )
      })
    }
  } else if (len_c == 3) {
    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]
    c3 <- data[, confounders[3]]

    df <- data.frame(X = x, Y = y, C1 = c1, C2 = c2, C3 = c3)
    df$pS <- plogis(s1_0 + s1_x * df$X + s1_y * df$Y)

    if (y_binary) {
      suppressWarnings({
        final <- glm(
          Y ~ X + C1 + C2 + C3,
          family = binomial(link = "logit"),
          weights = (1 / df$pS),
          data = df
        )
      })
    } else {
      suppressWarnings({
        final <- lm(
          Y ~ X + C1 + C2 + C3,
          weights = (1 / df$pS),
          data = df
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


#' Adust for selection bias.
#'
#' `adjust_sel` returns the exposure-outcome odds ratio and confidence
#' interval, adjusted for selection bias.
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
#' data, plus data for the selection indicator representing whether the
#' observation was selected in `data_observed`.
#' @param bias_params Object of class 'bias_params'
#' @param level Value from 0-1 representing the full range of the confidence
#' interval. Default is 0.95.
#'
#' @return A list where the first item is the odds ratio estimate of the
#' effect of the exposure on the outcome and the second item is the
#' confidence interval as the vector: (lower bound, upper bound).
#'
#' @examples
#' df_observed <- data_observed(
#'   data = df_sel,
#'   bias = "sel",
#'   exposure = "X",
#'   outcome = "Y",
#'   confounders = "C1"
#' )
#'
#' # Using validation data -----------------------------------------------------
#' df_validation <- data_validation(
#'   data = df_sel_source,
#'   true_exposure = "X",
#'   true_outcome = "Y",
#'   confounders = "C1",
#'   selection = "S"
#' )
#'
#' adjust_sel(
#'   data_observed = df_observed,
#'   data_validation = df_validation
#' )
#'
#' # Using bias_params ---------------------------------------------------------
#' bp <- bias_params(coef_list = list(s = c(0, 0.9, 0.9)))
#'
#' adjust_sel(
#'   data_observed = df_observed,
#'   s_model_coefs = c(0, 0.9, 0.9)
#' )
#'
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom stats binomial
#' @importFrom stats glm
#' @importFrom stats lm
#' @importFrom stats qnorm
#' @importFrom stats plogis
#' @importFrom stats coef
#' @importFrom rlang .data

adjust_sel <- function(
    data_observed,
    data_validation = NULL,
    bias_params = NULL,
    level = 0.95) {
  if (
    (!is.null(data_validation) && !is.null(bias_params)) ||
      (is.null(data_validation) && is.null(bias_params))
  ) {
    stop(
      "One of data_validation or bias_params must be non-null.",
      call. = FALSE
    )
  }
  data <- data_observed$data

  x <- data[, data_observed$exposure]
  y <- data[, data_observed$outcome]

  if (all(y %in% 0:1)) {
    y_binary <- TRUE
  } else {
    y_binary <- FALSE
  }

  if (!is.null(data_validation)) {
    final <- adjust_sel_val(
      data_observed,
      data_validation
    )
  } else if (!is.null(bias_params)) {
    if (is.null(bias_params$coef_list$s)) {
      stop(
        "bias_params must specify parameters for s to adjust for selection bias",
        call. = FALSE
      )
    }
    final <- adjust_sel_coef(
      data_observed,
      s_model_coefs = bias_params$coef_list$s
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
