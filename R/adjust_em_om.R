# Adjust for exposure misclassification and outcome misclassification

# the following functions feed into adjust_em_om():
# adjust_em_om_val() (data_validation input, method: imputation),
# adjust_em_om_coef_single() (bias_params input, method: imputation),
# adjust_em_om_coef_multinom() (bias_params input, method: weighting)


adjust_em_om_val <- function(
    data_observed,
    data_validation) {
  if (!all(data_observed$confounders %in% data_validation$confounders)) {
    stop("All confounders in observed data must be present in validation data.", call. = FALSE)
  }

  if (
    is.null(data_validation$misclassified_exposure) ||
      is.null(data_validation$misclassified_outcome)
  ) {
    stop(
      paste0(
        "Attempting to adjust for a misclassified exposure and misclassified outcome.",
        "\n",
        "Validation data must include a true exposure, misclassified exposure, true outcome, and misclassified outcome."
      ),
      call. = FALSE
    )
  }

  n <- nrow(data_observed$data)

  df <- data.frame(
    Xstar = data_observed$data[, data_observed$exposure],
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
    Xstar = data_validation$data[, data_validation$misclassified_exposure],
    Ystar = data_validation$data[, data_validation$misclassified_outcome]
  )
  df_val <- bind_cols(
    df_val,
    data_validation$data %>%
      select(all_of(data_validation$confounders))
  )

  force_binary(
    df$Xstar,
    "Exposure in observed data must be a binary integer."
  )
  force_binary(
    df$Ystar,
    "Outcome in observed data must be a binary integer."
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
    df_val$Ystar,
    "Misclassified outcome in validation data must be a binary integer."
  )
  force_binary(
    df_val$Y,
    "True outcome in validation data must be a binary integer."
  )

  x_mod <- glm(X ~ Xstar + Ystar + . - Y,
    family = binomial(link = "logit"),
    data = df_val
  )

  x_mod_coefs <- coef(x_mod)
  x_mod_se <- sqrt(diag(vcov(x_mod)))

  # Sample coefficients independently using their standard errors
  x_mod_coefs_sampled <- rnorm(length(x_mod_coefs), x_mod_coefs, x_mod_se)
  x_pred <- x_mod_coefs_sampled[1]

  for (i in 2:length(x_mod_coefs_sampled)) {
    var_name <- names(x_mod_coefs)[i]
    x_pred <- x_pred + df[[var_name]] * x_mod_coefs_sampled[i]
  }

  df$Xpred <- rbinom(n, 1, plogis(x_pred))

  y_mod <- glm(Y ~ X + Ystar + . - Xstar,
    family = binomial(link = "logit"),
    data = df_val
  )

  y_mod_coefs <- coef(y_mod)
  y_mod_se <- sqrt(diag(vcov(y_mod)))

  # Sample coefficients independently using their standard errors
  y_mod_coefs_sampled <- rnorm(length(y_mod_coefs), y_mod_coefs, y_mod_se)
  y_pred <- y_mod_coefs_sampled[1]

  for (i in 2:length(y_mod_coefs_sampled)) {
    var_name <- names(y_mod_coefs)[i]
    var_name <- gsub("X", "Xpred", var_name) # col X is not in df
    y_pred <- y_pred + df[[var_name]] * y_mod_coefs_sampled[i]
  }

  df$Ypred <- rbinom(n, 1, plogis(y_pred))

  final <- glm(
    Ypred ~ Xpred + . - Xstar - Ystar,
    family = binomial(link = "logit"),
    data = df
  )

  return(final)
}


adjust_em_om_coef_single <- function(
    data_observed,
    x_model_coefs,
    y_model_coefs) {
  data <- data_observed$data
  n <- nrow(data)
  confounders <- data_observed$confounders
  len_c <- length(confounders)
  len_x_coefs <- length(x_model_coefs)
  len_y_coefs <- length(y_model_coefs)

  xstar <- data[, data_observed$exposure]
  ystar <- data[, data_observed$outcome]

  force_binary(xstar, "Exposure must be a binary integer.")
  force_binary(ystar, "Outcome must be a binary integer.")
  force_len(
    len_x_coefs,
    3 + len_c,
    paste0(
      "Incorrect length of X model coefficients. ",
      "Length should equal 3 + number of confounders."
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
  x1_0 <- x_model_coefs[1]
  x1_xstar <- x_model_coefs[2]
  x1_ystar <- x_model_coefs[3]
  x_coefs_c <- x_model_coefs[4:len_x_coefs]

  y1_0 <- y_model_coefs[1]
  y1_x <- y_model_coefs[2]
  y1_ystar <- y_model_coefs[3]
  y_coefs_c <- y_model_coefs[4:len_y_coefs]

  # Create base dataframe
  df <- data.frame(Xstar = xstar, Ystar = ystar)

  # Add confounders if they exist
  if (!is.null(confounders)) {
    for (i in seq_along(confounders)) {
      df[[paste0("C", i)]] <- data[, confounders[i]]
    }
  }

  # Construct X prediction formula dynamically
  x_formula <- "x1_0 + x1_xstar * df$Xstar + x1_ystar * df$Ystar"
  if (!is.null(confounders)) {
    for (i in seq_along(confounders)) {
      x_formula <- paste0(x_formula, " + x_coefs_c[", i, "] * df$C", i)
    }
  }

  # Calculate X predictions
  df$Xpred <- rbinom(n, 1, plogis(eval(parse(text = x_formula))))

  # Construct Y prediction formula dynamically
  y_formula <- "y1_0 + y1_x * df$Xpred + y1_ystar * df$Ystar"
  if (!is.null(confounders)) {
    for (i in seq_along(confounders)) {
      y_formula <- paste0(y_formula, " + y_coefs_c[", i, "] * df$C", i)
    }
  }

  # Calculate Y predictions
  df$Ypred <- rbinom(n, 1, plogis(eval(parse(text = y_formula))))

  # Construct final model formula
  model_terms <- c("Xpred")
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


adjust_em_om_coef_multinom <- function(
    data_observed,
    x1y0_model_coefs,
    x0y1_model_coefs,
    x1y1_model_coefs) {
  data <- data_observed$data
  n <- nrow(data)
  confounders <- data_observed$confounders
  len_c <- length(confounders)
  len_x1y0_coefs <- length(x1y0_model_coefs)
  len_x0y1_coefs <- length(x0y1_model_coefs)
  len_x1y1_coefs <- length(x1y1_model_coefs)

  xstar <- data[, data_observed$exposure]
  ystar <- data[, data_observed$outcome]

  force_binary(xstar, "Exposure must be a binary integer.")
  force_binary(ystar, "Outcome must be a binary integer.")
  force_len(
    len_x1y0_coefs,
    3 + len_c,
    paste0(
      "Incorrect length of X1Y0 model coefficients. ",
      "Length should equal 3 + number of confounders."
    )
  )
  force_len(
    len_x0y1_coefs,
    3 + len_c,
    paste0(
      "Incorrect length of X0Y1 model coefficients. ",
      "Length should equal 3 + number of confounders."
    )
  )
  force_len(
    len_x1y1_coefs,
    3 + len_c,
    paste0(
      "Incorrect length of X1Y1 model coefficients. ",
      "Length should equal 3 + number of confounders."
    )
  )

  # Extract model coefficients
  x1y0_0 <- x1y0_model_coefs[1]
  x1y0_xstar <- x1y0_model_coefs[2]
  x1y0_ystar <- x1y0_model_coefs[3]
  x1y0_coefs_c <- x1y0_model_coefs[4:len_x1y0_coefs]

  x0y1_0 <- x0y1_model_coefs[1]
  x0y1_xstar <- x0y1_model_coefs[2]
  x0y1_ystar <- x0y1_model_coefs[3]
  x0y1_coefs_c <- x0y1_model_coefs[4:len_x0y1_coefs]

  x1y1_0 <- x1y1_model_coefs[1]
  x1y1_xstar <- x1y1_model_coefs[2]
  x1y1_ystar <- x1y1_model_coefs[3]
  x1y1_coefs_c <- x1y1_model_coefs[4:len_x1y1_coefs]

  # Create base dataframe
  df <- data.frame(Xstar = xstar, Ystar = ystar)

  # Add confounders if they exist
  if (!is.null(confounders)) {
    for (i in seq_along(confounders)) {
      df[[paste0("C", i)]] <- data[, confounders[i]]
    }
  }

  # Construct prediction formulas dynamically
  x1y0_formula <- "x1y0_0 + x1y0_xstar * df$Xstar + x1y0_ystar * df$Ystar"
  x0y1_formula <- "x0y1_0 + x0y1_xstar * df$Xstar + x0y1_ystar * df$Ystar"
  x1y1_formula <- "x1y1_0 + x1y1_xstar * df$Xstar + x1y1_ystar * df$Ystar"

  if (!is.null(confounders)) {
    for (i in seq_along(confounders)) {
      x1y0_formula <- paste0(x1y0_formula, " + x1y0_coefs_c[", i, "] * df$C", i)
      x0y1_formula <- paste0(x0y1_formula, " + x0y1_coefs_c[", i, "] * df$C", i)
      x1y1_formula <- paste0(x1y1_formula, " + x1y1_coefs_c[", i, "] * df$C", i)
    }
  }

  # Calculate predictions
  p_x1y0 <- exp(eval(parse(text = x1y0_formula)))
  p_x0y1 <- exp(eval(parse(text = x0y1_formula)))
  p_x1y1 <- exp(eval(parse(text = x1y1_formula)))

  denom <- (1 + p_x1y0 + p_x0y1 + p_x1y1)

  x0y0_pred <- 1 / denom
  x1y0_pred <- p_x1y0 / denom
  x0y1_pred <- p_x0y1 / denom
  x1y1_pred <- p_x1y1 / denom

  # Create prediction dataframe
  df_xy_pred <- data.frame(
    X0Y0 = x0y0_pred,
    X1Y0 = x1y0_pred,
    X0Y1 = x0y1_pred,
    X1Y1 = x1y1_pred
  )
  df_xy_pred4 <- bind_rows(df_xy_pred, df_xy_pred, df_xy_pred, df_xy_pred)

  # Create combined dataframe with all scenarios
  combined <- bind_rows(df, df, df, df) %>%
    bind_cols(df_xy_pred4) %>%
    mutate(
      Xbar = rep(c(1, 0, 1, 0), each = n),
      Ybar = rep(c(1, 1, 0, 0), each = n),
      pXY = case_when(
        Xbar == 0 & Ybar == 0 ~ X0Y0,
        Xbar == 1 & Ybar == 0 ~ X1Y0,
        Xbar == 0 & Ybar == 1 ~ X0Y1,
        Xbar == 1 & Ybar == 1 ~ X1Y1
      )
    )

  # Construct final model formula
  model_terms <- c("Xbar")
  if (!is.null(confounders)) {
    model_terms <- c(model_terms, paste0("C", seq_along(confounders)))
  }
  model_formula <- as.formula(
    paste("Ybar ~", paste(model_terms, collapse = " + "))
  )

  # Fit final model with weights
  suppressWarnings({
    final <- glm(
      model_formula,
      family = binomial(link = "logit"),
      weights = combined$pXY,
      data = combined
    )
  })

  return(final)
}


adjust_em_om <- function(
    data_observed,
    data_validation = NULL,
    bias_params = NULL,
    level = 0.95) {
  if (!is.null(data_validation)) {
    final <- adjust_em_om_val(
      data_observed,
      data_validation
    )
  } else if (!is.null(bias_params)) {
    if (all(c("x", "y") %in% names(bias_params$coef_list))) {
      final <- adjust_em_om_coef_single(
        data_observed,
        x_model_coefs = bias_params$coef_list$x,
        y_model_coefs = bias_params$coef_list$y
      )
    } else if (
      all(
        c("x1y0", "x0y1", "x1y1") %in%
          names(bias_params$coef_list)
      )
    ) {
      final <- adjust_em_om_coef_multinom(
        data_observed,
        x1y0_model_coefs = bias_params$coef_list$x1y0,
        x0y1_model_coefs = bias_params$coef_list$x0y1,
        x1y1_model_coefs = bias_params$coef_list$x1y1
      )
    } else {
      (
        stop(
          paste0(
            "bias_params must specify parameters for ",
            "exposure misclassification and outcome misclassification"
          ),
          call. = FALSE
        )
      )
    }
  }

  calculate_results(final, level, y_binary = TRUE)
}
