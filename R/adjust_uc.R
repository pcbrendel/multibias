# Adjust for uncontrolled confounding

# the following functions feed into adjust_uc():
# adjust_uc_val() (data_validation input, method: imputation),
# adjust_uc_coef() (bias_params input, method: imputation)

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
        "Attempting to adjust for unobserved confounding from one confounder.",
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
  u_mod_se <- sqrt(diag(vcov(u_mod)))

  # Sample coefficients independently using their standard errors
  u_mod_coefs_sampled <- rnorm(length(u_mod_coefs), u_mod_coefs, u_mod_se)
  u_pred <- u_mod_coefs_sampled[1]

  for (i in 2:length(u_mod_coefs_sampled)) {
    var_name <- names(u_mod_coefs)[i]
    u_pred <- u_pred + df[[var_name]] * u_mod_coefs_sampled[i]
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

  # Extract U model coefficients
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
  u_formula <- "u1_0 + u1_x * df$X + u1_y * df$Y"
  if (!is.null(confounders)) {
    for (i in seq_along(confounders)) {
      u_formula <- paste0(u_formula, " + u_coefs_c[", i, "] * df$C", i)
    }
  }

  # Calculate U predictions
  df$Upred <- rbinom(n, 1, plogis(eval(parse(text = u_formula))))

  # Construct final model formula
  model_terms <- c("X", "Upred")
  if (!is.null(confounders)) {
    model_terms <- c(model_terms, paste0("C", seq_along(confounders)))
  }
  model_formula <- as.formula(
    paste("Y ~", paste(model_terms, collapse = " + "))
  )

  # Fit final model
  if (y_binary) {
    final <- glm(
      model_formula,
      family = binomial(link = "logit"),
      data = df
    )
  } else {
    final <- lm(
      model_formula,
      data = df
    )
  }

  return(final)
}


adjust_uc <- function(
    data_observed,
    data_validation = NULL,
    bias_params = NULL,
    level = 0.95) {
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

  calculate_results(final, level, y_binary)
}
