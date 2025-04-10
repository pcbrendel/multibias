# Adjust for uncontrolled confounding and selection bias

# the following functions feed into adjust_uc_sel():
# adjust_uc_sel_val() (data_validation input),
# adjust_uc_sel_coef() (bias_params input)

adjust_uc_sel_val <- function(
    data_observed,
    data_validation) {
  if (!all(data_observed$confounders %in% data_validation$confounders)) {
    stop(
      "All confounders in observed data must be present in validation data.",
      call. = FALSE
    )
  }

  if (
    (length(data_validation$confounders) - length(data_observed$confounders) != 1) ||
      (is.null(data_validation$selection))
  ) {
    stop(
      paste0(
        "Attempting to adjust for unobserved confounding from one confounder and selection bias.",
        "\n",
        "Validation data must have: 1) one more confounder than the observed data, 2) a selection indicator column specified."
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

  # Extract model coefficients
  s1_0 <- s_model_coefs[1]
  s1_x <- s_model_coefs[2]
  s1_y <- s_model_coefs[3]

  u1_0 <- u_model_coefs[1]
  u1_x <- u_model_coefs[2]
  u1_y <- u_model_coefs[3]
  u_coefs_c <- u_model_coefs[4:len_u_coefs]

  # Create base dataframe
  df <- data.frame(X = x, Y = y)

  # Add confounders if they exist
  if (!is.null(confounders)) {
    for (i in seq_along(confounders)) {
      df[[paste0("C", i)]] <- data[, confounders[i]]
    }
  }

  # Construct U prediction formula dynamically
  u_formula <- "u1_0 + u1_x * x + u1_y * y"
  if (!is.null(confounders)) {
    for (i in seq_along(confounders)) {
      u_formula <- paste0(u_formula, " + u_coefs_c[", i, "] * df$C", i)
    }
  }

  # Calculate U predictions
  u1_pred <- plogis(eval(parse(text = u_formula)))
  u1_pred <- rep(u1_pred, times = 2)

  # Create combined dataframe with both U=0 and U=1 scenarios
  combined <- bind_rows(df, df) %>%
    mutate(
      Ubar = rep(c(1, 0), each = n),
      pU = case_when(
        Ubar == 1 ~ u1_pred,
        Ubar == 0 ~ 1 - u1_pred
      ),
      pS = plogis(s1_0 + s1_x * X + s1_y * Y)
    )

  # Construct final model formula
  model_terms <- c("X", "Ubar")
  if (!is.null(confounders)) {
    model_terms <- c(model_terms, paste0("C", seq_along(confounders)))
  }
  model_formula <- as.formula(paste("Y ~", paste(model_terms, collapse = " + ")))

  # Fit final model with weights
  if (y_binary) {
    suppressWarnings({
      final <- glm(
        model_formula,
        family = binomial(link = "logit"),
        weights = (combined$pU / combined$pS),
        data = combined
      )
    })
  } else {
    suppressWarnings({
      final <- lm(
        model_formula,
        weights = (combined$pU / combined$pS),
        data = combined
      )
    })
  }

  return(final)
}


adjust_uc_sel <- function(
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
    final <- adjust_uc_sel_val(
      data_observed,
      data_validation
    )
  } else if (!is.null(bias_params)) {
    if (is.null(bias_params$coef_list$u) && is.null(bias_params$coef_list$s)) {
      stop(
        paste0(
          "bias_params must specify parameters for uncontrolled ",
          "confounding and selection bias"
        ),
        call. = FALSE
      )
    }
    final <- adjust_uc_sel_coef(
      data_observed,
      u_model_coefs = bias_params$coef_list$u,
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
