# Adjust for exposure misclassification

# the following functions feed into adjust_em():
# adjust_em_val() (data_validation input, method: imputation),
# adjust_em_coef() (bias_params input, method: imputation)

adjust_em_val <- function(
    data_observed,
    data_validation) {
  if (!all(data_observed$confounders %in% data_validation$confounders)) {
    stop(
      "All confounders in observed data must be present in validation data.",
      call. = FALSE
    )
  }

  if (is.null(data_validation$misclassified_exposure)) {
    stop(
      paste0(
        "Attempting to adjust for a misclassified exposure.",
        "\n",
        "Validation data must have a true and misclassified exposure specified."
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
    Xstar = data_validation$data[, data_validation$misclassified_exposure]
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

  x_mod <- glm(X ~ Xstar + Y + .,
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

  if (y_binary) {
    final <- glm(
      Y ~ Xpred + . - Xstar,
      family = binomial(link = "logit"),
      data = df
    )
  } else {
    final <- lm(
      Y ~ Xpred + . - Xstar,
      data = df
    )
  }

  return(final)
}


adjust_em_coef <- function(
    data_observed,
    x_model_coefs) {
  data <- data_observed$data
  n <- nrow(data)
  confounders <- data_observed$confounders
  len_c <- length(confounders)
  len_x_coefs <- length(x_model_coefs)

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

  if (all(y %in% 0:1)) {
    y_binary <- TRUE
  } else {
    y_binary <- FALSE
  }

  # Extract X model coefficients
  x1_0 <- x_model_coefs[1]
  x1_xstar <- x_model_coefs[2]
  x1_y <- x_model_coefs[3]
  x_coefs_c <- x_model_coefs[4:len_x_coefs]

  # Create base dataframe
  df <- data.frame(Xstar = xstar, Y = y)

  # Add confounders if they exist
  if (!is.null(confounders)) {
    for (i in seq_along(confounders)) {
      df[[paste0("C", i)]] <- data[, confounders[i]]
    }
  }

  # Construct X prediction formula dynamically
  x_formula <- "x1_0 + x1_xstar * df$Xstar + x1_y * df$Y"
  if (!is.null(confounders)) {
    for (i in seq_along(confounders)) {
      x_formula <- paste0(x_formula, " + x_coefs_c[", i, "] * df$C", i)
    }
  }

  # Calculate X predictions
  df$Xpred <- rbinom(n, 1, plogis(eval(parse(text = x_formula))))

  # Construct final model formula
  model_terms <- c("Xpred")
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


adjust_em <- function(
    data_observed,
    data_validation = NULL,
    bias_params = NULL,
    level = 0.95) {
  data <- data_observed$data
  xstar <- data[, data_observed$exposure]
  y <- data[, data_observed$outcome]

  if (all(y %in% 0:1)) {
    y_binary <- TRUE
  } else {
    y_binary <- FALSE
  }

  if (!is.null(data_validation)) {
    final <- adjust_em_val(
      data_observed,
      data_validation
    )
  } else if (!is.null(bias_params)) {
    if (is.null(bias_params$coef_list$x)) {
      stop(
        "bias_params must specify parameters for x to adjust for exposure misclassification",
        call. = FALSE
      )
    }
    final <- adjust_em_coef(
      data_observed,
      x_model_coefs = bias_params$coef_list$x
    )
  }

  calculate_results(final, level, y_binary)
}
