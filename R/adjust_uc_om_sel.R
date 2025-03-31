adjust_uc_om_sel_val <- function(
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
  if (is.null(data_validation$misclassified_outcome)) {
    stop(
      paste0(
        "This function is adjusting for a misclassified outcome.",
        "\n",
        "Validation data must have a true and misclassified outcome specified."
      ),
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
    Ystar = data_validation$data[, data_validation$misclassified_outcome],
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
  force_binary(
    df_val$S,
    "Selection indicator in validation data must be a binary integer."
  )

  y_mod <- glm(Y ~ X + Ystar + . - U - S,
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

  s_mod <- glm(S ~ X + Ystar + . - Y - U,
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

  suppressWarnings({
    final <- glm(
      Ypred ~ X + Upred + . - Ystar - Spred,
      family = binomial(link = "logit"),
      weights = (1 / df$Spred),
      data = df
    )
  })

  return(final)
}


adjust_uc_om_sel_coef_single <- function(
    data_observed,
    u_model_coefs,
    y_model_coefs,
    s_model_coefs) {
  data <- data_observed$data
  n <- nrow(data)
  confounders <- data_observed$confounders
  len_c <- length(confounders)
  len_u_coefs <- length(u_model_coefs)
  len_y_coefs <- length(y_model_coefs)
  len_s_coefs <- length(s_model_coefs)

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
  force_len(
    len_s_coefs,
    3 + len_c,
    paste0(
      "Incorrect length of S model coefficients. ",
      "Length should equal 3 + number of confounders."
    )
  )

  u1_0 <- u_model_coefs[1]
  u1_x <- u_model_coefs[2]
  u1_y <- u_model_coefs[3]

  y1_0 <- y_model_coefs[1]
  y1_x <- y_model_coefs[2]
  y1_ystar <- y_model_coefs[3]

  s1_0 <- s_model_coefs[1]
  s1_x <- s_model_coefs[2]
  s1_ystar <- s_model_coefs[3]

  if (is.null(confounders)) {
    df <- data.frame(X = x, Ystar = ystar)

    df <- df %>%
      mutate(
        Ypred = rbinom(
          n, 1,
          plogis(y1_0 + y1_x * .data$X + y1_ystar * .data$Ystar)
        ),
        Upred = rbinom(
          n, 1,
          plogis(u1_0 + u1_x * .data$X + u1_y * .data$Ypred)
        ),
        pS = plogis(s1_0 + s1_x * .data$X + s1_ystar * .data$Ystar)
      )

    suppressWarnings({
      final <- glm(
        Ypred ~ X + Upred,
        family = binomial(link = "logit"),
        weights = (1 / df$pS),
        data = df
      )
    })
  } else if (len_c == 1) {
    c1 <- data[, confounders]
    df <- data.frame(X = x, Ystar = ystar, C1 = c1)

    y1_c1 <- y_model_coefs[4]
    s1_c1 <- s_model_coefs[4]

    df <- df %>%
      mutate(
        Ypred = rbinom(
          n, 1, plogis(
            y1_0 + y1_x * .data$X + y1_ystar * .data$Ystar +
              y1_c1 * .data$C1
          )
        ),
        Upred = rbinom(
          n, 1, plogis(u1_0 + u1_x * .data$X + u1_y * .data$Ypred)
        ),
        pS = plogis(
          s1_0 + s1_x * .data$X + s1_ystar * .data$Ystar +
            s1_c1 * .data$C1
        )
      )

    suppressWarnings({
      final <- glm(
        Ypred ~ X + C1 + Upred,
        family = binomial(link = "logit"),
        weights = (1 / df$pS),
        data = df
      )
    })
  } else if (len_c == 2) {
    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]

    df <- data.frame(X = x, Ystar = ystar, C1 = c1, C2 = c2)

    y1_c1 <- y_model_coefs[4]
    y1_c2 <- y_model_coefs[5]

    s1_c1 <- s_model_coefs[4]
    s1_c2 <- s_model_coefs[5]

    df <- df %>%
      mutate(
        Ypred = rbinom(
          n, 1,
          plogis(
            y1_0 + y1_x * .data$X + y1_ystar * .data$Ystar +
              y1_c1 * .data$C1 + y1_c2 * .data$C2
          )
        ),
        Upred = rbinom(
          n, 1, plogis(
            u1_0 + u1_x * .data$X +
              u1_y * .data$Ypred
          )
        ),
        pS = plogis(
          s1_0 + s1_x * .data$X + s1_ystar * .data$Ystar +
            s1_c1 * .data$C1 + s1_c2 * .data$C2
        )
      )

    suppressWarnings({
      final <- glm(
        Ypred ~ X + C1 + C2 + Upred,
        family = binomial(link = "logit"),
        weights = (1 / df$pS),
        data = df
      )
    })
  } else if (len_c == 3) {
    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]
    c3 <- data[, confounders[3]]

    df <- data.frame(X = x, Ystar = ystar, C1 = c1, C2 = c2, C3 = c3)

    y1_c1 <- y_model_coefs[4]
    y1_c2 <- y_model_coefs[5]
    y1_c3 <- y_model_coefs[6]

    s1_c1 <- s_model_coefs[4]
    s1_c2 <- s_model_coefs[5]
    s1_c3 <- s_model_coefs[6]

    df <- df %>%
      mutate(
        Ypred = rbinom(
          n, 1, plogis(
            y1_0 + y1_x * .data$X + y1_ystar * .data$Ystar +
              y1_c1 * .data$C1 + y1_c2 * .data$C2 + y1_c3 * .data$C3
          )
        ),
        Upred = rbinom(
          n, 1, plogis(u1_0 + u1_x * .data$X + u1_y * .data$Ypred)
        ),
        pS = plogis(
          s1_0 + s1_x * .data$X + s1_ystar * .data$Ystar +
            s1_c1 * .data$C1 + s1_c2 * .data$C2 + s1_c3 * .data$C3
        )
      )

    suppressWarnings({
      final <- glm(
        Ypred ~ X + C1 + C2 + C3 + Upred,
        family = binomial(link = "logit"),
        weights = (1 / df$pS),
        data = df
      )
    })
  } else if (len_c > 3) {
    stop(
      "This function is currently not compatible with >3 confounders.",
      call. = FALSE
    )
  }

  return(final)
}


adjust_uc_om_sel_coef_multinom <- function(
    data_observed,
    u1y0_model_coefs,
    u0y1_model_coefs,
    u1y1_model_coefs,
    s_model_coefs) {
  data <- data_observed$data
  n <- nrow(data)
  confounders <- data_observed$confounders
  len_c <- length(confounders)
  len_u0y1_coefs <- length(u0y1_model_coefs)
  len_u1y0_coefs <- length(u1y0_model_coefs)
  len_u1y1_coefs <- length(u1y1_model_coefs)
  len_s_coefs <- length(s_model_coefs)

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

  s1_0 <- s_model_coefs[1]
  s1_x <- s_model_coefs[2]
  s1_ystar <- s_model_coefs[3]

  u0y1_0 <- u0y1_model_coefs[1]
  u0y1_x <- u0y1_model_coefs[2]
  u0y1_ystar <- u0y1_model_coefs[3]

  u1y0_0 <- u1y0_model_coefs[1]
  u1y0_x <- u1y0_model_coefs[2]
  u1y0_ystar <- u1y0_model_coefs[3]

  u1y1_0 <- u1y1_model_coefs[1]
  u1y1_x <- u1y1_model_coefs[2]
  u1y1_ystar <- u1y1_model_coefs[3]

  if (is.null(confounders)) {
    df <- data.frame(X = x, Ystar = ystar)

    p_u0y1 <- exp(u0y1_0 + u0y1_x * df$X + u0y1_ystar * df$Ystar)
    p_u1y0 <- exp(u1y0_0 + u1y0_x * df$X + u1y0_ystar * df$Ystar)
    p_u1y1 <- exp(u1y1_0 + u1y1_x * df$X + u1y1_ystar * df$Ystar)

    denom <- (1 + p_u0y1 + p_u1y0 + p_u1y1)

    u0y0_pred <- 1 / denom
    u0y1_pred <- p_u0y1 / denom
    u1y0_pred <- p_u1y0 / denom
    u1y1_pred <- p_u1y1 / denom

    df_uy_pred <- data.frame(
      U0Y0 = u0y0_pred,
      U0Y1 = u0y1_pred,
      U1Y0 = u1y0_pred,
      U1Y1 = u1y1_pred
    )
    df_uy_pred4 <- bind_rows(df_uy_pred, df_uy_pred, df_uy_pred, df_uy_pred)

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
        ),
        pS = plogis(s1_0 + s1_x * .data$X + s1_ystar * .data$Ystar)
      )

    suppressWarnings({
      final <- glm(
        Ybar ~ X + Ubar,
        family = binomial(link = "logit"),
        weights = (combined$pUY / combined$pS),
        data = combined
      )
    })
  } else if (len_c == 1) {
    c1 <- data[, confounders]
    df <- data.frame(X = x, Ystar = ystar, C1 = c1)

    s1_c1 <- s_model_coefs[4]
    u1y0_c1 <- u1y0_model_coefs[4]
    u0y1_c1 <- u0y1_model_coefs[4]
    u1y1_c1 <- u1y1_model_coefs[4]

    p_u1y0 <- exp(
      u1y0_0 + u1y0_x * df$X + u1y0_ystar * df$Ystar +
        u1y0_c1 * df$C1
    )
    p_u0y1 <- exp(
      u0y1_0 + u0y1_x * df$X + u0y1_ystar * df$Ystar +
        u0y1_c1 * df$C1
    )
    p_u1y1 <- exp(
      u1y1_0 + u1y1_x * df$X + u1y1_ystar * df$Ystar +
        u1y1_c1 * df$C1
    )

    denom <- (1 + p_u1y0 + p_u0y1 + p_u1y1)

    u0y0_pred <- 1 / denom
    u1y0_pred <- p_u1y0 / denom
    u0y1_pred <- p_u0y1 / denom
    u1y1_pred <- p_u1y1 / denom

    df_uy_pred <- data.frame(
      U0Y0 = u0y0_pred,
      U1Y0 = u1y0_pred,
      U0Y1 = u0y1_pred,
      U1Y1 = u1y1_pred
    )
    df_uy_pred4 <- bind_rows(df_uy_pred, df_uy_pred, df_uy_pred, df_uy_pred)

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
        ),
        pS = plogis(
          s1_0 + s1_x * .data$X + s1_ystar * .data$Ystar +
            s1_c1 * .data$C1
        )
      )

    suppressWarnings({
      final <- glm(
        Ybar ~ X + C1 + Ubar,
        family = binomial(link = "logit"),
        weights = (combined$pUY / combined$pS),
        data = combined
      )
    })
  } else if (len_c == 2) {
    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]

    df <- data.frame(X = x, Ystar = ystar, C1 = c1, C2 = c2)

    s1_c1 <- s_model_coefs[4]
    s1_c2 <- s_model_coefs[5]

    u1y0_c1 <- u1y0_model_coefs[4]
    u1y0_c2 <- u1y0_model_coefs[5]

    u0y1_c1 <- u0y1_model_coefs[4]
    u0y1_c2 <- u0y1_model_coefs[5]

    u1y1_c1 <- u1y1_model_coefs[4]
    u1y1_c2 <- u1y1_model_coefs[5]

    p_u1y0 <- exp(
      u1y0_0 + u1y0_x * df$X + u1y0_ystar * df$Ystar +
        u1y0_c1 * df$C1 + u1y0_c2 * df$C2
    )
    p_u0y1 <- exp(
      u0y1_0 + u0y1_x * df$X + u0y1_ystar * df$Ystar +
        u0y1_c1 * df$C1 + u0y1_c2 * df$C2
    )
    p_u1y1 <- exp(
      u1y1_0 + u1y1_x * df$X + u1y1_ystar * df$Ystar +
        u1y1_c1 * df$C1 + u1y1_c2 * df$C2
    )

    denom <- (1 + p_u1y0 + p_u0y1 + p_u1y1)

    u0y0_pred <- 1 / denom
    u1y0_pred <- p_u1y0 / denom
    u0y1_pred <- p_u0y1 / denom
    u1y1_pred <- p_u1y1 / denom

    df_uy_pred <- data.frame(
      U0Y0 = u0y0_pred,
      U1Y0 = u1y0_pred,
      U0Y1 = u0y1_pred,
      U1Y1 = u1y1_pred
    )
    df_uy_pred4 <- bind_rows(df_uy_pred, df_uy_pred, df_uy_pred, df_uy_pred)

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
        ),
        pS = plogis(
          s1_0 + s1_x * .data$X + s1_ystar * .data$Ystar +
            s1_c1 * .data$C1 + s1_c2 * .data$C2
        )
      )

    suppressWarnings({
      final <- glm(
        Ybar ~ X + C1 + C2 + Ubar,
        family = binomial(link = "logit"),
        weights = (combined$pUY / combined$pS),
        data = combined
      )
    })
  } else if (len_c == 3) {
    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]
    c3 <- data[, confounders[3]]

    df <- data.frame(X = x, Ystar = ystar, C1 = c1, C2 = c2, C3 = c3)

    s1_c1 <- s_model_coefs[4]
    s1_c2 <- s_model_coefs[5]
    s1_c3 <- s_model_coefs[6]

    u1y0_c1 <- u1y0_model_coefs[4]
    u1y0_c2 <- u1y0_model_coefs[5]
    u1y0_c3 <- u1y0_model_coefs[6]

    u0y1_c1 <- u0y1_model_coefs[4]
    u0y1_c2 <- u0y1_model_coefs[5]
    u0y1_c3 <- u0y1_model_coefs[6]

    u1y1_c1 <- u1y1_model_coefs[4]
    u1y1_c2 <- u1y1_model_coefs[5]
    u1y1_c3 <- u1y1_model_coefs[6]

    p_u1y0 <- exp(
      u1y0_0 + u1y0_x * df$X + u1y0_ystar * df$Ystar +
        u1y0_c1 * df$C1 + u1y0_c2 * df$C2 + u1y0_c3 * df$C3
    )
    p_u0y1 <- exp(
      u0y1_0 + u0y1_x * df$X + u0y1_ystar * df$Ystar +
        u0y1_c1 * df$C1 + u0y1_c2 * df$C2 + u0y1_c3 * df$C3
    )
    p_u1y1 <- exp(
      u1y1_0 + u1y1_x * df$X + u1y1_ystar * df$Ystar +
        u1y1_c1 * df$C1 + u1y1_c2 * df$C2 + u1y1_c3 * df$C3
    )

    denom <- (1 + p_u1y0 + p_u0y1 + p_u1y1)

    u0y0_pred <- 1 / denom
    u1y0_pred <- p_u1y0 / denom
    u0y1_pred <- p_u0y1 / denom
    u1y1_pred <- p_u1y1 / denom

    df_uy_pred <- data.frame(
      U0Y0 = u0y0_pred,
      U1Y0 = u1y0_pred,
      U0Y1 = u0y1_pred,
      U1Y1 = u1y1_pred
    )
    df_uy_pred4 <- bind_rows(df_uy_pred, df_uy_pred, df_uy_pred, df_uy_pred)

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
        ),
        pS = plogis(
          s1_0 + s1_x * .data$X + s1_ystar * .data$Ystar +
            s1_c1 * .data$C1 + s1_c2 * .data$C2 + s1_c3 * .data$C3
        )
      )

    suppressWarnings({
      final <- glm(
        Ybar ~ X + C1 + C2 + C3 + Ubar,
        family = binomial(link = "logit"),
        weights = (combined$pUY / combined$pS),
        data = combined
      )
    })
  } else if (len_c > 3) {
    stop(
      "This function is currently not compatible with >3 confounders.",
      call. = FALSE
    )
  }

  return(final)
}


#' Adust for uncontrolled confounding, outcome misclassification, and selection
#' bias.
#'
#' `adjust_uc_om_sel` returns the exposure-outcome odds ratio and
#' confidence interval, adjusted for uncontrolled confounding, outcome
#' misclassificaiton, and selection bias.
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
#' data, plus data for: 1) the true and misclassified outcome corresponding
#' to the observed outcome in `data_observed`, 2) the confounder missing in
#' `data_observed`, 3) a selection indicator representing whether the
#' observation in `data_validation` was selected in `data_observed`.
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
#'   data = df_uc_om_sel,
#'   exposure = "X",
#'   outcome = "Ystar",
#'   confounders = c("C1", "C2", "C3")
#' )
#'
#' # Using validation data -----------------------------------------------------
#' df_validation <- data_validation(
#'   data = df_uc_om_sel_source,
#'   true_exposure = "X",
#'   true_outcome = "Y",
#'   confounders = c("C1", "C2", "C3", "U"),
#'   misclassified_outcome = "Ystar",
#'   selection = "S"
#' )
#'
#' adjust_uc_om_sel(
#'   data_observed = df_observed,
#'   data_validation = df_validation
#' )
#'
#' # Using bias_params ---------------------------------------------------------
#' bp1 <- bias_params(
#'   coef_list = list(
#'     u = c(-0.32, 0.59, 0.69),
#'     y = c(-2.85, 0.71, 1.63, 0.40, -0.85, 0.22),
#'     s = c(0.00, 0.74, 0.19, 0.02, -0.06, 0.02)
#'   )
#' )
#'
#' adjust_uc_om_sel(
#'   data_observed = df_observed,
#'   bias_params = bp1
#' )
#'
#' bp2 <- bias_params(
#'   coef_list = list(
#'     u1y0 = c(-0.20, 0.62, 0.01, -0.08, 0.10, -0.15),
#'     u0y1 = c(-3.28, 0.63, 1.65, 0.42, -0.85, 0.26),
#'     u1y1 = c(-2.70, 1.22, 1.64, 0.32, -0.77, 0.09),
#'     s = c(0.00, 0.74, 0.19, 0.02, -0.06, 0.02)
#'   )
#' )
#'
#' adjust_uc_om_sel(
#'   data_observed = df_observed,
#'   bias_params = bp2
#' )
#'
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom stats binomial
#' @importFrom stats glm
#' @importFrom stats qnorm
#' @importFrom stats rbinom
#' @importFrom stats plogis
#' @importFrom rlang .data
#'
#' @export

adjust_uc_om_sel <- function(
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
    final <- adjust_uc_om_sel_val(
      data_observed,
      data_validation
    )
  } else if (!is.null(bias_params)) {
    if (all(c("u", "y", "s") %in% names(bias_params$coef_list))) {
      final <- adjust_uc_om_sel_coef_single(
        data_observed,
        u_model_coefs = bias_params$coef_list$u,
        y_model_coefs = bias_params$coef_list$y,
        s_model_coefs = bias_params$coef_list$s
      )
    } else if (
      all(
        c("u1y0", "u0y1", "u1y1", "s") %in%
          names(bias_params$coef_list)
      )
    ) {
      final <- adjust_uc_om_sel_coef_multinom(
        data_observed,
        u1y0_model_coefs = bias_params$coef_list$u1y0,
        u0y1_model_coefs = bias_params$coef_list$u0y1,
        u1y1_model_coefs = bias_params$coef_list$u1y1,
        s_model_coefs = bias_params$coef_list$s
      )
    } else {
      (
        stop(
          paste0(
            "bias_params must specify parameters for uncontrolled ",
            "confounding, outcome misclassification, and selection bias"
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
