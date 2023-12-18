adjust_multinom_emc_omc <- function(
  data,
  exposure,
  outcome,
  confounders = NULL,
  x1y0_model_coefs,
  x0y1_model_coefs,
  x1y1_model_coefs,
  level = 0.95
) {

  n <- nrow(data)
  len_c <- length(confounders)
  len_x1y0_coefs <- length(x1y0_model_coefs)
  len_x0y1_coefs <- length(x0y1_model_coefs)
  len_x1y1_coefs <- length(x1y1_model_coefs)

  xstar <- data[, exposure]
  ystar <- data[, outcome]

  if (sum(xstar %in% c(0, 1)) != n) {
    stop("Exposure must be a binary integer.")
  }
  if (sum(ystar %in% c(0, 1)) != n) {
    stop("Outcome must be a binary integer.")
  }
  if (len_x1y0_coefs != 3 + len_c) {
    stop(
      paste0(
        "Incorrect length of X1Y0 model coefficients. ",
        "Length should equal 3 + number of confounders."
      )
    )
  }
  if (len_x0y1_coefs != 3 + len_c) {
    stop(
      paste0(
        "Incorrect length of X0Y1 model coefficients. ",
        "Length should equal 3 + number of confounders."
      )
    )
  }
  if (len_x1y1_coefs != 3 + len_c) {
    stop(
      paste0(
        "Incorrect length of X1Y1 model coefficients. ",
        "Length should equal 3 + number of confounders."
      )
    )
  }

  x1y0_0     <- x1y0_model_coefs[1]
  x1y0_xstar <- x1y0_model_coefs[2]
  x1y0_ystar <- x1y0_model_coefs[3]

  x0y1_0     <- x0y1_model_coefs[1]
  x0y1_xstar <- x0y1_model_coefs[2]
  x0y1_ystar <- x0y1_model_coefs[3]

  x1y1_0     <- x1y1_model_coefs[1]
  x1y1_xstar <- x1y1_model_coefs[2]
  x1y1_ystar <- x1y1_model_coefs[3]

  if (is.null(confounders)) {

    df <- data.frame(Xstar = xstar, Ystar = ystar)

    p_x1y0 <- exp(x1y0_0 + x1y0_xstar * df$Xstar + x1y0_ystar * df$Ystar)
    p_x0y1 <- exp(x0y1_0 + x0y1_xstar * df$Xstar + x0y1_ystar * df$Ystar)
    p_x1y1 <- exp(x1y1_0 + x1y1_xstar * df$Xstar + x1y1_ystar * df$Ystar)

    denom <- (1 + p_x1y0 + p_x0y1 + p_x1y1)

    x0y0_pred <- 1 / denom
    x1y0_pred <- p_x1y0 / denom
    x0y1_pred <- p_x0y1 / denom
    x1y1_pred <- p_x1y1 / denom

    df_xy_pred <- data.frame(
      X0Y0 = x0y0_pred,
      X1Y0 = x1y0_pred,
      X0Y1 = x0y1_pred,
      X1Y1 = x1y1_pred
    )
    df_xy_pred4 <- bind_rows(df_xy_pred, df_xy_pred, df_xy_pred, df_xy_pred)

    combined <- bind_rows(df, df, df, df) %>%
      bind_cols(df_xy_pred4) %>%
      mutate(Xbar = rep(c(1, 0, 1, 0), each = n),
             Ybar = rep(c(1, 1, 0, 0), each = n),
             pXY = case_when(Xbar == 0 & Ybar == 0 ~ X0Y0,
                             Xbar == 1 & Ybar == 0 ~ X1Y0,
                             Xbar == 0 & Ybar == 1 ~ X0Y1,
                             Xbar == 1 & Ybar == 1 ~ X1Y1))
    suppressWarnings({
      final <- glm(
        Ybar ~ Xbar,
        family = binomial(link = "logit"),
        weights = combined$pXY,
        data = combined
      )
    })

  } else if (len_c == 1) {

    c1 <- data[, confounders]

    df <- data.frame(Xstar = xstar, Ystar = ystar, C1 = c1)

    x1y0_c1 <- x1y0_model_coefs[4]
    x0y1_c1 <- x0y1_model_coefs[4]
    x1y1_c1 <- x1y1_model_coefs[4]

    p_x1y0 <- exp(x1y0_0 + x1y0_xstar * df$Xstar + x1y0_ystar * df$Ystar +
                    x1y0_c1 * df$C1)
    p_x0y1 <- exp(x0y1_0 + x0y1_xstar * df$Xstar + x0y1_ystar * df$Ystar +
                    x0y1_c1 * df$C1)
    p_x1y1 <- exp(x1y1_0 + x1y1_xstar * df$Xstar + x1y1_ystar * df$Ystar +
                    x1y1_c1 * df$C1)

    denom <- (1 + p_x1y0 + p_x0y1 + p_x1y1)

    x0y0_pred <- 1 / denom
    x1y0_pred <- p_x1y0 / denom
    x0y1_pred <- p_x0y1 / denom
    x1y1_pred <- p_x1y1 / denom

    df_xy_pred <- data.frame(
      X0Y0 = x0y0_pred,
      X1Y0 = x1y0_pred,
      X0Y1 = x0y1_pred,
      X1Y1 = x1y1_pred
    )
    df_xy_pred4 <- bind_rows(df_xy_pred, df_xy_pred, df_xy_pred, df_xy_pred)

    combined <- bind_rows(df, df, df, df) %>%
      bind_cols(df_xy_pred4) %>%
      mutate(Xbar = rep(c(1, 0, 1, 0), each = n),
             Ybar = rep(c(1, 1, 0, 0), each = n),
             pXY = case_when(Xbar == 0 & Ybar == 0 ~ X0Y0,
                             Xbar == 1 & Ybar == 0 ~ X1Y0,
                             Xbar == 0 & Ybar == 1 ~ X0Y1,
                             Xbar == 1 & Ybar == 1 ~ X1Y1))

    suppressWarnings({
      final <- glm(
        Ybar ~ Xbar + C1,
        family = binomial(link = "logit"),
        weights = combined$pXY,
        data = combined
      )
    })

  } else if (len_c == 2) {

    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]

    df <- data.frame(Xstar = xstar, Ystar = ystar, C1 = c1, C2 = c2)

    x1y0_c1 <- x1y0_model_coefs[4]
    x1y0_c2 <- x1y0_model_coefs[5]

    x0y1_c1 <- x0y1_model_coefs[4]
    x0y1_c2 <- x0y1_model_coefs[5]

    x1y1_c1 <- x1y1_model_coefs[4]
    x1y1_c2 <- x1y1_model_coefs[5]

    p_x1y0 <- exp(x1y0_0 + x1y0_xstar * df$Xstar + x1y0_ystar * df$Ystar +
                    x1y0_c1 * df$C1 + x1y0_c2 * df$C2)
    p_x0y1 <- exp(x0y1_0 + x0y1_xstar * df$Xstar + x0y1_ystar * df$Ystar +
                    x0y1_c1 * df$C1 + x0y1_c2 * df$C2)
    p_x1y1 <- exp(x1y1_0 + x1y1_xstar * df$Xstar + x1y1_ystar * df$Ystar +
                    x1y1_c1 * df$C1 + x1y1_c2 * df$C2)

    denom <- (1 + p_x1y0 + p_x0y1 + p_x1y1)

    x0y0_pred <- 1 / denom
    x1y0_pred <- p_x1y0 / denom
    x0y1_pred <- p_x0y1 / denom
    x1y1_pred <- p_x1y1 / denom

    df_xy_pred <- data.frame(
      X0Y0 = x0y0_pred,
      X1Y0 = x1y0_pred,
      X0Y1 = x0y1_pred,
      X1Y1 = x1y1_pred
    )
    df_xy_pred4 <- bind_rows(df_xy_pred, df_xy_pred, df_xy_pred, df_xy_pred)

    combined <- bind_rows(df, df, df, df) %>%
      bind_cols(df_xy_pred4) %>%
      mutate(Xbar = rep(c(1, 0, 1, 0), each = n),
             Ybar = rep(c(1, 1, 0, 0), each = n),
             pXY = case_when(Xbar == 0 & Ybar == 0 ~ X0Y0,
                             Xbar == 1 & Ybar == 0 ~ X1Y0,
                             Xbar == 0 & Ybar == 1 ~ X0Y1,
                             Xbar == 1 & Ybar == 1 ~ X1Y1))

    suppressWarnings({
      final <- glm(
        Ybar ~ Xbar + C1 + C2,
        family = binomial(link = "logit"),
        weights = combined$pXY,
        data = combined
      )
    })

  } else if (len_c == 3) {

    c1 <- data[, confounders[1]]
    c2 <- data[, confounders[2]]
    c3 <- data[, confounders[3]]

    df <- data.frame(Xstar = xstar, Ystar = ystar, C1 = c1, C2 = c2, C3 = c3)

    x1y0_c1 <- x1y0_model_coefs[4]
    x1y0_c2 <- x1y0_model_coefs[5]
    x1y0_c3 <- x1y0_model_coefs[6]

    x0y1_c1 <- x0y1_model_coefs[4]
    x0y1_c2 <- x0y1_model_coefs[5]
    x0y1_c3 <- x0y1_model_coefs[6]

    x1y1_c1 <- x1y1_model_coefs[4]
    x1y1_c2 <- x1y1_model_coefs[5]
    x1y1_c3 <- x1y1_model_coefs[6]

    p_x1y0 <- exp(x1y0_0 + x1y0_xstar * df$Xstar + x1y0_ystar * df$Ystar +
                    x1y0_c1 * df$C1 + x1y0_c2 * df$C2 + x1y0_c3 * df$C3)
    p_x0y1 <- exp(x0y1_0 + x0y1_xstar * df$Xstar + x0y1_ystar * df$Ystar +
                    x0y1_c1 * df$C1 + x0y1_c2 * df$C2 + x0y1_c3 * df$C3)
    p_x1y1 <- exp(x1y1_0 + x1y1_xstar * df$Xstar + x1y1_ystar * df$Ystar +
                    x1y1_c1 * df$C1 + x1y1_c2 * df$C2 + x1y1_c3 * df$C3)

    denom <- (1 + p_x1y0 + p_x0y1 + p_x1y1)

    x0y0_pred <- 1 / denom
    x1y0_pred <- p_x1y0 / denom
    x0y1_pred <- p_x0y1 / denom
    x1y1_pred <- p_x1y1 / denom

    df_xy_pred <- data.frame(
      X0Y0 = x0y0_pred,
      X1Y0 = x1y0_pred,
      X0Y1 = x0y1_pred,
      X1Y1 = x1y1_pred
    )
    df_xy_pred4 <- bind_rows(df_xy_pred, df_xy_pred, df_xy_pred, df_xy_pred)

    combined <- bind_rows(df, df, df, df) %>%
      bind_cols(df_xy_pred4) %>%
      mutate(Xbar = rep(c(1, 0, 1, 0), each = n),
             Ybar = rep(c(1, 1, 0, 0), each = n),
             pXY = case_when(Xbar == 0 & Ybar == 0 ~ X0Y0,
                             Xbar == 1 & Ybar == 0 ~ X1Y0,
                             Xbar == 0 & Ybar == 1 ~ X0Y1,
                             Xbar == 1 & Ybar == 1 ~ X1Y1))

    suppressWarnings({
      final <- glm(
        Ybar ~ Xbar + C1 + C2 + C3,
        family = binomial(link = "logit"),
        weights = combined$pXY,
        data = combined
      )
    })

  } else if (len_c > 3) {

    stop("This function is currently not compatible with >3 confounders.")

  }

  est <- summary(final)$coef[2, 1]
  se <- summary(final)$coef[2, 2]
  alpha <- 1 - level

  estimate <- exp(est)
  ci <- c(exp(est + se * qnorm(alpha / 2)),
          exp(est + se * qnorm(1 - alpha / 2)))
  return(list(estimate = estimate, ci = ci))

}