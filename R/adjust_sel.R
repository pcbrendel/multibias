# Adjust for selection bias

# the following functions feed into adjust_sel():
# adjust_sel_val() (data_validation input, method: weighting),
# adjust_sel_coef() (bias_params input, method: weighting)

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
        "Attempting to adjust for selection bias.",
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

  # Extract S model coefficients
  s1_0 <- s_model_coefs[1]
  s1_x <- s_model_coefs[2]
  s1_y <- s_model_coefs[3]

  # Create base dataframe
  df <- data.frame(X = x, Y = y)

  # Add confounders if they exist
  if (!is.null(confounders)) {
    for (i in seq_along(confounders)) {
      df[[paste0("C", i)]] <- data[, confounders[i]]
    }
  }

  # Calculate selection probabilities
  df$pS <- plogis(s1_0 + s1_x * df$X + s1_y * df$Y)

  # Construct final model formula
  model_terms <- c("X")
  if (!is.null(confounders)) {
    model_terms <- c(model_terms, paste0("C", seq_along(confounders)))
  }
  model_formula <- as.formula(
    paste("Y ~", paste(model_terms, collapse = " + "))
  )

  # Fit final model with weights
  suppressWarnings({
    if (y_binary) {
      final <- glm(
        model_formula,
        family = binomial(link = "logit"),
        weights = (1 / df$pS),
        data = df
      )
    } else {
      final <- lm(
        model_formula,
        weights = (1 / df$pS),
        data = df
      )
    }
  })

  return(final)
}


adjust_sel <- function(
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

  calculate_results(final, level, y_binary)
}
