# Adjust for outcome misclassification

# the following functions feed into adjust_om():
# adjust_om_val() (data_validation input),
# adjust_om_coef() (bias_params input)

adjust_om_val <- function(
    data_observed,
    data_validation) {
  if (!all(data_observed$confounders %in% data_validation$confounders)) {
    stop(
      "All confounders in observed data must be present in validation data.",
      call. = FALSE
    )
  }

  if (is.null(data_validation$misclassified_outcome)) {
    stop(
      paste0(
        "Attempting to adjust for a misclassified outcome.",
        "\n",
        "Validation data must have a true and misclassified outcome specified."
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
  df_val <- bind_cols(
    df_val,
    data_validation$data %>%
      select(all_of(data_validation$confounders))
  )

  force_match(
    df$X,
    df_val$X,
    "Exposures from both datasets must both be binary or both be continuous."
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

  # Extract Y model coefficients
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

  # Construct final model formula
  model_terms <- c("X")
  if (!is.null(confounders)) {
    model_terms <- c(model_terms, paste0("C", seq_along(confounders)))
  }
  model_formula <- as.formula(
    paste("Ypred ~", paste(model_terms, collapse = " + "))
  )

  # Fit final model
  final <- glm(
    model_formula,
    family = binomial(link = "logit"),
    data = df
  )

  return(final)
}


adjust_om <- function(
    data_observed,
    data_validation = NULL,
    bias_params = NULL,
    level = 0.95) {
  if (!is.null(data_validation)) {
    final <- adjust_om_val(
      data_observed,
      data_validation
    )
  } else if (!is.null(bias_params)) {
    if (is.null(bias_params$coef_list$y)) {
      stop(
        "bias_params must specify parameters for y to adjust for outcome misclassification",
        call. = FALSE
      )
    }
    final <- adjust_om_coef(
      data_observed,
      y_model_coefs = bias_params$coef_list$y
    )
  }

  calculate_results(final, level, y_binary = TRUE)
}
