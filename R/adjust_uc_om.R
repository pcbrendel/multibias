# Adjust for uncontrolled confounding and outcome misclassification

# the following functions feed into adjust_uc_om():
# adjust_uc_om_val() (data_validation input),
# adjust_uc_om_coef_single() (bias_params input),
# adjust_uc_om_coef_multinom() (bias_params input)

adjust_uc_om_val <- function(
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
      (is.null(data_validation$misclassified_outcome))
  ) {
    stop(
      paste0(
        "Attempting to adjust for unobserved confounding from one confounder and outcome misclassification.",
        "\n",
        "Validation data must have: 1) one more confounder than the observed data, 2) a true and misclassified outcome specified."
      ),
      call. = FALSE
    )
  }

  n <- nrow(data_observed$data)

  df <- data.frame(
    X = data_observed$data[, data_observed$exposure],
    Ystar = data_observed$data[, data_observed$outcome]
  )
  df <- bind_cols(
    df,
    data_observed$data %>%
      select(all_of(data_observed$confounders))
  )

  df_val <- data.frame(
    X = data_validation$data[, data_validation$true_exposure],
    Y = data_validation$data[, data_validation$true_outcome],
    Ystar = data_validation$data[, data_validation$misclassified_outcome]
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
    "Outcomes from both datasets must both be binary or both be continuous."
  )
  force_binary(
    df_val$U,
    "Uncontrolled confounder in validation data must be a binary integer."
  )
  force_binary(
    df$Ystar,
    "Outcome in observed data must be a binary integer."
  )
  force_binary(
    df_val$Ystar,
    "Misclassified outcome in validation data must be a binary integer."
  )
  force_binary(
    df_val$Y,
    "True outcome in validation data must be a binary integer."
  )

  y_mod <- glm(Y ~ X + Ystar + . - U,
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

  u_mod <- glm(U ~ X + Y,
    family = binomial(link = "logit"),
    data = df_val
  )

  u_mod_coefs <- coef(u_mod)
  u_pred <- u_mod_coefs[1]

  for (i in 2:length(u_mod_coefs)) {
    var_name <- names(u_mod_coefs)[i]
    var_name <- gsub("Y", "Ypred", var_name) # col Y is not in df
    u_pred <- u_pred + df[[var_name]] * u_mod_coefs[i]
  }

  df$Upred <- rbinom(n, 1, plogis(u_pred))

  final <- glm(
    Ypred ~ X + Upred + . - Ystar,
    family = binomial(link = "logit"),
    data = df
  )

  return(final)
}


adjust_uc_om_coef_single <- function(
    data_observed,
    u_model_coefs,
    y_model_coefs) {
  data <- data_observed$data
  n <- nrow(data)
  confounders <- data_observed$confounders
  len_c <- length(confounders)
  len_u_coefs <- length(u_model_coefs)
  len_y_coefs <- length(y_model_coefs)

  x <- data[, data_observed$exposure]
  ystar <- data[, data_observed$outcome]

  force_binary(ystar, "Outcome must be a binary integer.")
  force_len(
    len_u_coefs,
    3,
    paste0(
      "Incorrect length of U model coefficients. ",
      "Length should equal 3."
    )
  )
  force_len(
    len_y_coefs,
    3 + len_c,
    paste0(
      "Incorrect length of Y model coefficients. ",
      "Length should equal 3 + number of confounders."
    )
  )

  # Extract model coefficients
  u1_0 <- u_model_coefs[1]
  u1_x <- u_model_coefs[2]
  u1_y <- u_model_coefs[3]

  y1_0 <- y_model_coefs[1]
  y1_x <- y_model_coefs[2]
  y1_ystar <- y_model_coefs[3]
  y_coefs_c <- y_model_coefs[4:len_y_coefs]

  # Create base dataframe
  df <- data.frame(X = x, Ystar = ystar)

  # Add confounders if they exist
  if (!is.null(confounders)) {
    for (i in seq_along(confounders)) {
      df[[paste0("C", i)]] <- data[, confounders[i]]
    }
  }

  # Construct Y prediction formula dynamically
  y_formula <- "y1_0 + y1_x * df$X + y1_ystar * df$Ystar"
  if (!is.null(confounders)) {
    for (i in seq_along(confounders)) {
      y_formula <- paste0(y_formula, " + y_coefs_c[", i, "] * df$C", i)
    }
  }

  # Calculate Y predictions
  df$Ypred <- rbinom(n, 1, plogis(eval(parse(text = y_formula))))

  # Construct U prediction formula
  u_formula <- "u1_0 + u1_x * df$X + u1_y * df$Ypred"

  # Calculate U predictions
  df$Upred <- rbinom(n, 1, plogis(eval(parse(text = u_formula))))

  # Construct final model formula
  model_terms <- c("X", "Upred")
  if (!is.null(confounders)) {
    model_terms <- c(model_terms, paste0("C", seq_along(confounders)))
  }
  model_formula <- as.formula(paste("Ypred ~", paste(model_terms, collapse = " + ")))

  # Fit final model
  final <- glm(
    model_formula,
    family = binomial(link = "logit"),
    data = df
  )

  return(final)
}


adjust_uc_om_coef_multinom <- function(
    data_observed,
    u1y0_model_coefs,
    u0y1_model_coefs,
    u1y1_model_coefs) {
  data <- data_observed$data
  n <- nrow(data)
  confounders <- data_observed$confounders
  len_c <- length(confounders)
  len_u1y0_coefs <- length(u1y0_model_coefs)
  len_u0y1_coefs <- length(u0y1_model_coefs)
  len_u1y1_coefs <- length(u1y1_model_coefs)

  x <- data[, data_observed$exposure]
  ystar <- data[, data_observed$outcome]

  force_binary(ystar, "Outcome must be a binary integer.")
  force_len(
    len_u1y0_coefs,
    3 + len_c,
    paste0(
      "Incorrect length of U1Y0 model coefficients. ",
      "Length should equal 3 + number of confounders."
    )
  )
  force_len(
    len_u0y1_coefs,
    3 + len_c,
    paste0(
      "Incorrect length of U0Y1 model coefficients. ",
      "Length should equal 3 + number of confounders."
    )
  )
  force_len(
    len_u1y1_coefs,
    3 + len_c,
    paste0(
      "Incorrect length of U1Y1 model coefficients. ",
      "Length should equal 3 + number of confounders."
    )
  )

  # Extract model coefficients
  u1y0_0 <- u1y0_model_coefs[1]
  u1y0_x <- u1y0_model_coefs[2]
  u1y0_ystar <- u1y0_model_coefs[3]
  u1y0_coefs_c <- u1y0_model_coefs[4:len_u1y0_coefs]

  u0y1_0 <- u0y1_model_coefs[1]
  u0y1_x <- u0y1_model_coefs[2]
  u0y1_ystar <- u0y1_model_coefs[3]
  u0y1_coefs_c <- u0y1_model_coefs[4:len_u0y1_coefs]

  u1y1_0 <- u1y1_model_coefs[1]
  u1y1_x <- u1y1_model_coefs[2]
  u1y1_ystar <- u1y1_model_coefs[3]
  u1y1_coefs_c <- u1y1_model_coefs[4:len_u1y1_coefs]

  # Create base dataframe
  df <- data.frame(X = x, Ystar = ystar)

  # Add confounders if they exist
  if (!is.null(confounders)) {
    for (i in seq_along(confounders)) {
      df[[paste0("C", i)]] <- data[, confounders[i]]
    }
  }

  # Construct prediction formulas dynamically
  u1y0_formula <- "u1y0_0 + u1y0_x * df$X + u1y0_ystar * df$Ystar"
  u0y1_formula <- "u0y1_0 + u0y1_x * df$X + u0y1_ystar * df$Ystar"
  u1y1_formula <- "u1y1_0 + u1y1_x * df$X + u1y1_ystar * df$Ystar"

  if (!is.null(confounders)) {
    for (i in seq_along(confounders)) {
      u1y0_formula <- paste0(u1y0_formula, " + u1y0_coefs_c[", i, "] * df$C", i)
      u0y1_formula <- paste0(u0y1_formula, " + u0y1_coefs_c[", i, "] * df$C", i)
      u1y1_formula <- paste0(u1y1_formula, " + u1y1_coefs_c[", i, "] * df$C", i)
    }
  }

  # Calculate predictions
  p_u1y0 <- exp(eval(parse(text = u1y0_formula)))
  p_u0y1 <- exp(eval(parse(text = u0y1_formula)))
  p_u1y1 <- exp(eval(parse(text = u1y1_formula)))

  denom <- (1 + p_u1y0 + p_u0y1 + p_u1y1)

  u0y0_pred <- 1 / denom
  u1y0_pred <- p_u1y0 / denom
  u0y1_pred <- p_u0y1 / denom
  u1y1_pred <- p_u1y1 / denom

  # Create prediction dataframe
  df_uy_pred <- data.frame(
    U0Y0 = u0y0_pred,
    U1Y0 = u1y0_pred,
    U0Y1 = u0y1_pred,
    U1Y1 = u1y1_pred
  )
  df_uy_pred4 <- bind_rows(df_uy_pred, df_uy_pred, df_uy_pred, df_uy_pred)

  # Create combined dataframe with all scenarios
  combined <- bind_rows(df, df, df, df) %>%
    bind_cols(df_uy_pred4) %>%
    mutate(
      Ubar = rep(c(1, 0, 1, 0), each = n),
      Ybar = rep(c(1, 1, 0, 0), each = n),
      pUY = case_when(
        Ubar == 0 & Ybar == 0 ~ U0Y0,
        Ubar == 1 & Ybar == 0 ~ U1Y0,
        Ubar == 0 & Ybar == 1 ~ U0Y1,
        Ubar == 1 & Ybar == 1 ~ U1Y1
      )
    )

  # Construct final model formula
  model_terms <- c("X", "Ubar")
  if (!is.null(confounders)) {
    model_terms <- c(model_terms, paste0("C", seq_along(confounders)))
  }
  model_formula <- as.formula(paste("Ybar ~", paste(model_terms, collapse = " + ")))

  # Fit final model with weights
  suppressWarnings({
    final <- glm(
      model_formula,
      family = binomial(link = "logit"),
      weights = combined$pUY,
      data = combined
    )
  })

  return(final)
}


adjust_uc_om <- function(
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

  if (!is.null(data_validation)) {
    final <- adjust_uc_om_val(
      data_observed,
      data_validation
    )
  } else if (!is.null(bias_params)) {
    if (all(c("u", "y") %in% names(bias_params$coef_list))) {
      final <- adjust_uc_om_coef_single(
        data_observed,
        u_model_coefs = bias_params$coef_list$u,
        y_model_coefs = bias_params$coef_list$y
      )
    } else if (
      all(
        c("u1y0", "u0y1", "u1y1") %in%
          names(bias_params$coef_list)
      )
    ) {
      final <- adjust_uc_om_coef_multinom(
        data_observed,
        u1y0_model_coefs = bias_params$coef_list$u1y0,
        u0y1_model_coefs = bias_params$coef_list$u0y1,
        u1y1_model_coefs = bias_params$coef_list$u1y1
      )
    } else {
      (
        stop(
          paste0(
            "bias_params must specify parameters for ",
            "uncontrolled confounding and outcome misclassification"
          ),
          call. = FALSE
        )
      )
    }
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
