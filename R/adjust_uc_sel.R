adjust_uc_sel_val <- function(
    data_observed,
    data_validation) {
  if (!all(data_observed$confounders %in% data_validation$confounders)) {
    stop("All confounders in observed data must be present in validation data.")
  }

  if (
    length(data_validation$confounders) - length(data_observed$confounders) != 1
  ) {
    stop(
      paste0(
        "This function adjusts for unobserved confounding from one confounder.",
        "\n",
        "Validation data must have one more confounder than the observed data."
      )
    )
  }
  if (is.null(data_validation$selection)) {
    stop(
      paste0(
        "This function is adjusting for selection bias.",
        "\n",
        "Validation data must have a selection indicator column specified."
      )
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
  force_binary(
    df_val$S,
    "Selection indicator in validation data must be a binary integer."
  )

  u_mod <- glm(U ~ X + Y + . - S,
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
        Y ~ X + Upred + . - Spred,
        family = binomial(link = "logit"),
        weights = (1 / df$Spred),
        data = df
      )
    })
  } else {
    suppressWarnings({
      final <- lm(
        Y ~ X + Upred + . - Spred,
        weights = (1 / df$Spred),
        data = df
      )
    })
  }

  return(final)
}


adjust_uc_sel_coef <- function(
    data_observed,
    u_model_coefs,
    s_model_coefs) {
  data <- data_observed$data
  n <- nrow(data)
  confounders <- data_observed$confounders
  len_c <- length(confounders)
  len_u_coefs <- length(u_model_coefs)
  len_s_coefs <- length(s_model_coefs)

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

  u1_0 <- u_model_coefs[1]
  u1_x <- u_model_coefs[2]
  u1_y <- u_model_coefs[3]

  if (is.null(confounders)) {
    df <- data.frame(X = x, Y = y)

    u1_pred <- plogis(u1_0 + u1_x * x + u1_y * y)
    u1_pred <- rep(u1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(
        Ubar = rep(c(1, 0), each = n),
        pS = plogis(s1_0 + s1_x * .data$X + s1_y * .data$Y),
        pU = case_when(
          Ubar == 1 ~ u1_pred,
          Ubar == 0 ~ 1 - u1_pred
        )
      )

    if (y_binary) {
      suppressWarnings({
        final <- glm(
          Y ~ X + Ubar,
          family = binomial(link = "logit"),
          weights = (combined$pU / combined$pS),
          data = combined
        )
      })
    } else {
      suppressWarnings({
        final <- lm(
          Y ~ X + Ubar,
          weights = (combined$pU / combined$pS),
          data = combined
        )
      })
    }
  } else if (len_c == 1) {
    c1 <- data[, confounders]
    df <- data.frame(X = x, Y = y, C1 = c1)
    u1_c1 <- u_model_coefs[4]

    u1_pred <- plogis(u1_0 + u1_x * x + u1_y * y + u1_c1 * c1)
    u1_pred <- rep(u1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(
        Ubar = rep(c(1, 0), each = n),
        pS = plogis(s1_0 + s1_x * .data$X + s1_y * .data$Y),
        pU = case_when(
          Ubar == 1 ~ u1_pred,
          Ubar == 0 ~ 1 - u1_pred
        )
      )

    if (y_binary) {
      suppressWarnings({
        final <- glm(
          Y ~ X + C1 + Ubar,
          family = binomial(link = "logit"),
          weights = (combined$pU / combined$pS),
          data = combined
        )
      })
    } else {
      suppressWarnings({
        final <- lm(
          Y ~ X + C1 + Ubar,
          weights = (combined$pU / combined$pS),
          data = combined
        )
      })
    }
  } else if (len_c == 2) {
    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]

    df <- data.frame(X = x, Y = y, C1 = c1, C2 = c2)

    u1_c1 <- u_model_coefs[4]
    u1_c2 <- u_model_coefs[5]

    u1_pred <- plogis(u1_0 + u1_x * x + u1_y * y + u1_c1 * c1 + u1_c2 * c2)
    u1_pred <- rep(u1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(
        Ubar = rep(c(1, 0), each = n),
        pS = plogis(s1_0 + s1_x * .data$X + s1_y * .data$Y),
        pU = case_when(
          Ubar == 1 ~ u1_pred,
          Ubar == 0 ~ 1 - u1_pred
        )
      )

    if (y_binary) {
      suppressWarnings({
        final <- glm(
          Y ~ X + C1 + C2 + Ubar,
          family = binomial(link = "logit"),
          weights = (combined$pU / combined$pS),
          data = combined
        )
      })
    } else {
      suppressWarnings({
        final <- lm(
          Y ~ X + C1 + C2 + Ubar,
          weights = (combined$pU / combined$pS),
          data = combined
        )
      })
    }
  } else if (len_c == 3) {
    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]
    c3 <- data[, confounders[3]]

    df <- data.frame(X = x, Y = y, C1 = c1, C2 = c2, C3 = c3)

    u1_c1 <- u_model_coefs[4]
    u1_c2 <- u_model_coefs[5]
    u1_c3 <- u_model_coefs[6]

    u1_pred <- plogis(
      u1_0 + u1_x * x + u1_y * y +
        u1_c1 * c1 + u1_c2 * c2 + u1_c3 * c3
    )
    u1_pred <- rep(u1_pred, times = 2)

    combined <- bind_rows(df, df) %>%
      mutate(
        Ubar = rep(c(1, 0), each = n),
        pS = plogis(s1_0 + s1_x * .data$X + s1_y * .data$Y),
        pU = case_when(
          Ubar == 1 ~ u1_pred,
          Ubar == 0 ~ 1 - u1_pred
        )
      )

    if (y_binary) {
      suppressWarnings({
        final <- glm(
          Y ~ X + C1 + C2 + C3 + Ubar,
          family = binomial(link = "logit"),
          weights = (combined$pU / combined$pS),
          data = combined
        )
      })
    } else {
      suppressWarnings({
        final <- lm(
          Y ~ X + C1 + C2 + C3 + Ubar,
          weights = (combined$pU / combined$pS),
          data = combined
        )
      })
    }
  } else if (len_c > 3) {
    stop("This function is currently not compatible with >3 confounders.")
  }

  return(final)
}


#' Adust for uncontrolled confounding and selection bias.
#'
#' `adjust_uc_sel` returns the exposure-outcome odds ratio and confidence
#' interval, adjusted for uncontrolled confounding and exposure
#' misclassificaiton.
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
#' data, plus data for the confounder missing in `data_observed`. There
#' should also be a selection indicator representing whether the observation in
#' `data_validation` was selected in `data_observed`.
#' @param u_model_coefs The regression coefficients corresponding to the model:
#' \ifelse{html}{\out{logit(P(U=1)) = &alpha;<sub>0</sub> + &alpha;<sub>1</sub>X + &alpha;<sub>2</sub>Y + &alpha;<sub>2+j</sub>C<sub>j</sub>, }}{\eqn{logit(P(U=1)) = \alpha_0 + \alpha_1 X + \alpha_2 Y + \alpha_{2+j} C_j, }}
#' where *U* is the binary unmeasured
#' confounder, *X* is the exposure, *Y* is the outcome, *C*
#' represents the vector of measured confounders (if any), and *j*
#' corresponds to the number of measured confounders. The number of parameters
#' therefore equals 3 + *j*.
#' @param s_model_coefs The regression coefficients corresponding to the model:
#' \ifelse{html}{\out{logit(P(S=1)) = &beta;<sub>0</sub> + &beta;<sub>1</sub>X + &beta;<sub>2</sub>Y, }}{\eqn{logit(P(S=1)) = \beta_0 + \beta_1 X + \beta_2 Y, }}
#' where *S* represents binary selection, *X* is the exposure,
#' and *Y* is the outcome. The number of parameters therefore equals 3.
#' @param level Value from 0-1 representing the full range of the confidence
#' interval. Default is 0.95.
#
#' @return A list where the first item is the odds ratio estimate of the
#' effect of the exposure on the outcome and the second item is the
#' confidence interval as the vector: (lower bound, upper bound).
#'
#' @examples
#' df_observed <- data_observed(
#'   data = df_uc_sel,
#'   exposure = "X",
#'   outcome = "Y",
#'   confounders = c("C1", "C2", "C3")
#' )
#' # Using validation data -----------------------------------------------------
#' df_validation <- data_validation(
#'   data = df_uc_sel_source,
#'   true_exposure = "X",
#'   true_outcome = "Y",
#'   confounders = c("C1", "C2", "C3", "U"),
#'   selection = "S"
#' )
#'
#' adjust_uc_sel(
#'   data_observed = df_observed,
#'   data_validation = df_validation
#' )
#'
#' # Using u_model_coefs and s_model_coefs -------------------------------------
#' adjust_uc_sel(
#'   data_observed = df_observed,
#'   u_model_coefs = c(-0.19, 0.61, 0.72, -0.09, 0.10, -0.15),
#'   s_model_coefs = c(-0.01, 0.92, 0.94)
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

adjust_uc_sel <- function(
    data_observed,
    data_validation = NULL,
    u_model_coefs = NULL,
    s_model_coefs = NULL,
    level = 0.95) {
  if (!is.null(data_validation)) {
    if (!all(is.null(u_model_coefs), is.null(s_model_coefs))) {
      stop("No bias parameters should be specified when 'data_validation' is used.")
    }
  } else if (!is.null(u_model_coefs) && !is.null(s_model_coefs)) {
    if (!is.null(data_validation)) {
      stop("No other bias-adjusting inputs should be specified when 'u_model_coefs' and 's_model_coefs' are used.")
    }
  } else {
    stop(
      paste(
        "One of:",
        "1. data_validation",
        "2. (u_model_coefs & s_model_coefs)",
        "must be non-null.",
        sep = "\n"
      )
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
    final <- adjust_uc_sel_val(
      data_observed,
      data_validation
    )
  } else if (!is.null(u_model_coefs)) {
    final <- adjust_uc_sel_coef(
      data_observed,
      u_model_coefs,
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
