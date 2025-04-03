adjust_uc_val <- function(
    data_observed,
    data_validation) {
  if (!all(data_observed$confounders %in% data_validation$confounders)) {
    stop(
      "All confounders in observed data must be present in validation data.",
      call. = FALSE
    )
  }

  if (
    length(data_validation$confounders) - length(data_observed$confounders) != 1
  ) {
    stop(
      paste0(
        "This function adjusts for unobserved confounding from one confounder.",
        "\n",
        "Validation data must have one more confounder than the observed data."
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
    Y = data_validation$data[, data_validation$true_outcome]
  )

  uc <- setdiff(data_validation$confounders, data_observed$confounders)
  df_val$U <- data_validation$data[, uc]
  df_val <- bind_cols(
    df_val,
    data_validation$data %>%
      select(all_of(data_observed$confounders))
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
    df_val$U,
    "Uncontrolled confounder in validation data must be a binary integer."
  )

  u_mod <- glm(U ~ X + Y + .,
    family = binomial(link = "logit"),
    data = df_val
  )

  u_mod_coefs <- coef(u_mod)
  u_pred <- u_mod_coefs[1]

  for (i in 2:length(u_mod_coefs)) {
    var_name <- names(u_mod_coefs)[i]
    u_pred <- u_pred + df[[var_name]] * u_mod_coefs[i]
  }

  df$Upred <- rbinom(n, 1, plogis(u_pred))

  if (y_binary) {
    final <- glm(
      Y ~ X + Upred + .,
      family = binomial(link = "logit"),
      data = df
    )
  } else {
    final <- lm(
      Y ~ X + Upred + .,
      data = df
    )
  }

  return(final)
}


adjust_uc_coef <- function(
    data_observed,
    u_model_coefs) {
  data <- data_observed$data
  n <- nrow(data)
  confounders <- data_observed$confounders
  len_c <- length(confounders)
  len_u_coefs <- length(u_model_coefs)

  x <- data[, data_observed$exposure]
  y <- data[, data_observed$outcome]

  force_len(
    len_u_coefs,
    3 + len_c,
    paste0(
      "Incorrect length of U model coefficients. ",
      "Length should equal 3 + number of confounders."
    )
  )

  if (all(y %in% 0:1)) {
    y_binary <- TRUE
  } else {
    y_binary <- FALSE
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
    stop(
      "This function is currently not compatible with >3 confounders.",
      call. = FALSE
    )
  }

  return(final)
}


#' Adust for uncontrolled confounding.
#'
#' `adjust_uc` returns the exposure-outcome odds ratio and confidence
#' interval, adjusted for uncontrolled confounding from a binary confounder.
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
#' data, plus data for the confounder missing in `data_observed`.
#' @param bias_params Object of class 'bias_params'
#' @param level Value from 0-1 representing the full range of the confidence
#' interval. Default is 0.95.
#'
#' @return A list where the first item is the odds ratio estimate of the
#' effect of the exposure on the outcome and the second item is the
#' confidence interval as the vector: (lower bound, upper bound).
#'
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom stats binomial
#' @importFrom stats glm
#' @importFrom stats lm
#' @importFrom stats qnorm
#' @importFrom stats rbinom
#' @importFrom stats plogis
#' @importFrom stats coef
#' @importFrom rlang .data

adjust_uc <- function(
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
    final <- adjust_uc_val(
      data_observed,
      data_validation
    )
  } else if (!is.null(bias_params)) {
    if (is.null(bias_params$coef_list$u)) {
      stop(
        "bias_params must specify parameters for u to adjust for uncontrolled confounding",
        call. = FALSE
      )
    }
    final <- adjust_uc_coef(
      data_observed,
      u_model_coefs = bias_params$coef_list$u
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
