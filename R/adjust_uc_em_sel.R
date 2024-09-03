#' Adust for uncontrolled confounding, exposure misclassification, and selection
#' bias.
#'
#' \code{adjust_uc_em_sel} returns the exposure-outcome odds ratio and
#' confidence interval, adjusted for uncontrolled confounding, exposure
#' misclassificaiton, and selection bias. Two different options for the bias
#' parameters are availale here: 1) parameters from separate models
#' of \emph{U} and \emph{X} (\code{u_model_coefs} and \code{x_model_coefs})
#' or 2) parameters from a joint model of \emph{U} and \emph{X}
#' (\code{x1u0_model_coefs}, \code{x0u1_model_coefs}, and
#' \code{x1u1_model_coefs}). Both approaches require \code{s_model_coefs}.
#'
#' Values for the regression coefficients can be applied as
#' fixed values or as single draws from a probability
#' distribution (ex: \code{rnorm(1, mean = 2, sd = 1)}). The latter has
#' the advantage of allowing the researcher to capture the uncertainty
#' in the bias parameter estimates. To incorporate this uncertainty in the
#' estimate and confidence interval, this function should be run in loop across
#' bootstrap samples of the dataframe for analysis. The estimate and
#' confidence interval would then be obtained from the median and quantiles
#' of the distribution of odds ratio estimates.
#'
#' @inheritParams adjust_em_sel
#' @param u_model_coefs The regression coefficients corresponding to the model:
#' \ifelse{html}{\out{logit(P(U=1)) = &alpha;<sub>0</sub> + &alpha;<sub>1</sub>X + &alpha;<sub>2</sub>Y, }}{\eqn{logit(P(U=1)) = \alpha_0 + \alpha_1 X + \alpha_2 Y, }}
#' where \emph{U} is the binary unmeasured confounder, \emph{X} is the
#' binary true exposure, and \emph{Y} is the outcome.
#' The number of parameters therefore equals 3.
#' @param x_model_coefs The regression coefficients corresponding to the model:
#' \ifelse{html}{\out{logit(P(X=1)) = &delta;<sub>0</sub> + &delta;<sub>1</sub>X* + &delta;<sub>2</sub>Y + &delta;<sub>2+j</sub>C<sub>j</sub>, }}{\eqn{logit(P(X=1)) = \delta_0 + \delta_1 X^* + \delta_2 Y + \delta_{2+j} C_j, }}
#' where \emph{X} represents binary true exposure, \emph{X*} is the
#' binary misclassified exposure, \emph{Y} is the outcome, \emph{C}
#' represents the vector of measured confounders (if any), and
#' \emph{j} corresponds to the number of measured confounders.
#' The number of parameters therefore equals 3 + \emph{j}.
#' @param x1u0_model_coefs The regression coefficients corresponding to the
#' model:
#' \ifelse{html}{\out{log(P(X=1,U=0)/P(X=0,U=0)) = &gamma;<sub>1,0</sub> + &gamma;<sub>1,1</sub>X* + &gamma;<sub>1,2</sub>Y + &gamma;<sub>1,2+j</sub>C<sub>j</sub>, }}{\eqn{log(P(X=1,U=0)/P(X=0,U=0)) = \gamma_{1,0} + \gamma_{1,1} X^* + \gamma_{1,2} Y + \gamma_{1,2+j} C_j, }}
#' where \emph{X} is the binary true exposure, \emph{U} is the binary
#' unmeasured confounder, \emph{X*} is the binary misclassified exposure,
#' \emph{Y} is the outcome, \emph{C}
#' represents the vector of measured confounders (if any), and
#' \emph{j} corresponds to the number of measured confounders.
#' @param x0u1_model_coefs The regression coefficients corresponding to the
#' model:
#' \ifelse{html}{\out{log(P(X=0,U=1)/P(X=0,U=0)) = &gamma;<sub>2,0</sub> + &gamma;<sub>2,1</sub>X* + &gamma;<sub>2,2</sub>Y + &gamma;<sub>2,2+j</sub>C<sub>j</sub>, }}{\eqn{log(P(X=0,U=1)/P(X=0,U=0)) = \gamma_{2,0} + \gamma_{2,1} X^* + \gamma_{2,2} Y + \gamma_{2,2+j} C_j, }}
#' where \emph{X} is the binary true exposure, \emph{U} is the binary
#' unmeasured confounder, \emph{X*} is the binary misclassified exposure,
#' \emph{Y} is the outcome,
#' \emph{C} represents the vector of measured confounders (if any), and
#' \emph{j} corresponds to the number of measured confounders.
#' @param x1u1_model_coefs The regression coefficients corresponding to the
#' model:
#' \ifelse{html}{\out{log(P(X=1,U=1)/P(X=0,U=0)) = &gamma;<sub>3,0</sub> + &gamma;<sub>3,1</sub>X* + &gamma;<sub>3,2</sub>Y + &gamma;<sub>3,2+j</sub>C<sub>j</sub>, }}{\eqn{log(P(X=1,U=1)/P(X=0,U=0)) = \gamma_{3,0} + \gamma_{3,1} X^* + \gamma_{3,2} Y + \gamma_{3,2+j} C_j, }}
#' where \emph{X} is the binary true exposure, \emph{U} is the binary
#' unmeasured confounder, \emph{X*} is the binary misclassified exposure,
#' \emph{Y} is the outcome,
#' \emph{C} represents the vector of measured confounders (if any),
#' and \emph{j} corresponds to the number of measured confounders.
#' @param s_model_coefs The regression coefficients corresponding to the model:
#' \ifelse{html}{\out{logit(P(S=1)) = &beta;<sub>0</sub> + &beta;<sub>1</sub>X* + &beta;<sub>2</sub>Y + &beta;<sub>2+j</sub>C<sub>2+j</sub>, }}{\eqn{logit(P(S=1)) = \beta_0 + \beta_1 X^* + \beta_2 Y + \beta_{2+j} C_j, }}
#' where \emph{S} represents binary selection, \emph{X*} is the
#' binary misclassified exposure, \emph{Y} is the outcome,
#' \emph{C} represents the vector of measured confounders (if any),
#' and \emph{j} corresponds to the number of measured confounders.
#' The number of parameters therefore equals 3 + \emph{j}.
#'
#' @return A list where the first item is the odds ratio estimate of the
#' effect of the exposure on the outcome and the second item is the
#' confidence interval as the vector: (lower bound, upper bound).
#'
#' @examples
#' # Using u_model_coefs, x_model_coefs, s_model_coefs -------------------------
#' adjust_uc_em_sel(
#'   df_uc_em_sel,
#'   exposure = "Xstar",
#'   outcome = "Y",
#'   confounders = c("C1", "C2", "C3"),
#'   u_model_coefs = c(-0.32, 0.59, 0.69),
#'   x_model_coefs = c(-2.44, 1.62, 0.72, 0.32, -0.15, 0.85),
#'   s_model_coefs = c(0.00, 0.26, 0.78, 0.03, -0.02, 0.10)
#' )
#'
#' # Using x1u0_model_coefs, x0u1_model_coefs, x1u1_model_coefs, s_model_coefs
#' adjust_uc_em_sel(
#'   df_uc_em_sel,
#'   exposure = "Xstar",
#'   outcome = "Y",
#'   confounders = c("C1", "C2", "C3"),
#'   x1u0_model_coefs = c(-2.78, 1.62, 0.61, 0.36, -0.27, 0.88),
#'   x0u1_model_coefs = c(-0.17, -0.01, 0.71, -0.08, 0.07, -0.15),
#'   x1u1_model_coefs = c(-2.36, 1.62, 1.29, 0.25, -0.06, 0.74),
#'   s_model_coefs = c(0.00, 0.26, 0.78, 0.03, -0.02, 0.10)
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

adjust_uc_em_sel <- function(
    data,
    exposure,
    outcome,
    confounders = NULL,
    u_model_coefs = NULL,
    x_model_coefs = NULL,
    x1u0_model_coefs = NULL,
    x0u1_model_coefs = NULL,
    x1u1_model_coefs = NULL,
    s_model_coefs,
    level = 0.95) {
  n <- nrow(data)
  len_c <- length(confounders)

  xstar <- data[, exposure]
  y <- data[, outcome]

  if (sum(xstar %in% c(0, 1)) != n) {
    stop("Exposure must be a binary integer.")
  }
  if (sum(y %in% c(0, 1)) == n) {
    y_binary <- TRUE
  } else {
    y_binary <- FALSE
  }

  # check that user correctly specified bias parameters
  if (
    !(
      (is.null(u_model_coefs) && is.null(x_model_coefs)) ||
        ((is.null(x1u0_model_coefs) && is.null(x0u1_model_coefs) &&
            is.null(x1u1_model_coefs)))
    )
  ) {
    stop(
      "In addition to s_model_coefs, bias parameters must be specified for:\n
       1) u_model_coefs and x_model_coefs OR\n
       2) x1u0_model_coefs, x0u1_model_coefs, and x1u1_model_coefs."
    )
  }

  if (!is.null(x_model_coefs)) {
    len_u_coefs <- length(u_model_coefs)
    len_x_coefs <- length(x_model_coefs)
    len_s_coefs <- length(s_model_coefs)

    if (len_u_coefs != 3) {
      stop("Incorrect length of U model coefficients. Length should equal 3.")
    }
    if (len_x_coefs != 3 + len_c) {
      stop(
        paste0(
          "Incorrect length of X model coefficients. ",
          "Length should equal 3 + number of confounders."
        )
      )
    }
    if (len_s_coefs != 3 + len_c) {
      stop(
        paste0(
          "Incorrect length of S model coefficients. ",
          "Length should equal 3 + number of confounders."
        )
      )
    }

    u1_0 <- u_model_coefs[1]
    u1_x <- u_model_coefs[2]
    u1_y <- u_model_coefs[3]

    x1_0 <- x_model_coefs[1]
    x1_xstar <- x_model_coefs[2]
    x1_y <- x_model_coefs[3]

    s1_0 <- s_model_coefs[1]
    s1_xstar <- s_model_coefs[2]
    s1_y <- s_model_coefs[3]

    if (is.null(confounders)) {
      df <- data.frame(Xstar = xstar, Y = y)

      df <- df %>%
        mutate(
          Xpred = rbinom(
            n, 1,
            plogis(x1_0 + x1_xstar * .data$Xstar + x1_y * .data$Y)
          ),
          Upred = rbinom(
            n, 1,
            plogis(u1_0 + u1_x * .data$Xpred + u1_y * .data$Y)
          ),
          pS = plogis(s1_0 + s1_xstar * .data$Xstar + s1_y * .data$Y)
        )

      if (y_binary) {
        suppressWarnings({
          final <- glm(
            Y ~ Xpred + Upred,
            family = binomial(link = "logit"),
            weights = (1 / df$pS),
            data = df
          )
        })
      } else {
        final <- lm(
          Y ~ Xpred + Upred,
          weights = (1 / df$pS),
          data = df
        )
      }
    } else if (len_c == 1) {
      c1 <- data[, confounders]
      df <- data.frame(Xstar = xstar, Y = y, C1 = c1)

      x1_c1 <- x_model_coefs[4]
      s1_c1 <- s_model_coefs[4]

      df <- df %>%
        mutate(
          Xpred = rbinom(
            n, 1, plogis(
              x1_0 + x1_xstar * .data$Xstar +
                x1_y * .data$Y + x1_c1 * .data$C1
            )
          ),
          Upred = rbinom(
            n, 1, plogis(
              u1_0 + u1_x * .data$Xpred +
                u1_y * .data$Y
            )
          ),
          pS = plogis(
            s1_0 + s1_xstar * .data$Xstar + s1_y * .data$Y +
              s1_c1 * .data$C1
          )
        )

      if (y_binary) {
        suppressWarnings({
          final <- glm(
            Y ~ Xpred + C1 + Upred,
            family = binomial(link = "logit"),
            weights = (1 / df$pS),
            data = df
          )
        })
      } else {
        final <- lm(
          Y ~ Xpred + C1 + Upred,
          weights = (1 / df$pS),
          data = df
        )
      }
    } else if (len_c == 2) {
      c1 <- data[, confounders[1]]
      c2 <- data[, confounders[2]]

      df <- data.frame(Xstar = xstar, Y = y, C1 = c1, C2 = c2)

      x1_c1 <- x_model_coefs[4]
      x1_c2 <- x_model_coefs[5]

      s1_c1 <- s_model_coefs[4]
      s1_c2 <- s_model_coefs[5]

      df <- df %>%
        mutate(
          Xpred = rbinom(
            n, 1,
            plogis(
              x1_0 + x1_xstar * .data$Xstar + x1_y * .data$Y +
                x1_c1 * .data$C1 + x1_c2 * .data$C2
            )
          ),
          Upred = rbinom(
            n, 1, plogis(
              u1_0 + u1_x * .data$Xpred +
                u1_y * .data$Y
            )
          ),
          pS = plogis(
            s1_0 + s1_xstar * .data$Xstar + s1_y * .data$Y +
              s1_c1 * .data$C1 + s1_c2 * .data$C2
          )
        )

      if (y_binary) {
        suppressWarnings({
          final <- glm(
            Y ~ Xpred + C1 + C2 + Upred,
            family = binomial(link = "logit"),
            weights = (1 / df$pS),
            data = df
          )
        })
      } else {
        final <- lm(
          Y ~ Xpred + C1 + C2 + Upred,
          weights = (1 / df$pS),
          data = df
        )
      }
    } else if (len_c == 3) {
      c1 <- data[, confounders[1]]
      c2 <- data[, confounders[2]]
      c3 <- data[, confounders[3]]

      df <- data.frame(Xstar = xstar, Y = y, C1 = c1, C2 = c2, C3 = c3)

      x1_c1 <- x_model_coefs[4]
      x1_c2 <- x_model_coefs[5]
      x1_c3 <- x_model_coefs[6]

      s1_c1 <- s_model_coefs[4]
      s1_c2 <- s_model_coefs[5]
      s1_c3 <- s_model_coefs[6]

      df <- df %>%
        mutate(
          Xpred = rbinom(
            n, 1,
            plogis(
              x1_0 + x1_xstar * .data$Xstar + x1_y * .data$Y +
                x1_c1 * .data$C1 + x1_c2 * .data$C2 + x1_c3 * .data$C3
            )
          ),
          Upred = rbinom(
            n, 1, plogis(
              u1_0 + u1_x * .data$Xpred +
                u1_y * .data$Y
            )
          ),
          pS = plogis(
            s1_0 + s1_xstar * .data$Xstar + s1_y * .data$Y +
              s1_c1 * .data$C1 + s1_c2 * .data$C2 + s1_c3 * .data$C3
          )
        )

      if (y_binary) {
        suppressWarnings({
          final <- glm(
            Y ~ Xpred + C1 + C2 + C3 + Upred,
            family = binomial(link = "logit"),
            weights = (1 / df$pS),
            data = df
          )
        })
      } else {
        final <- lm(
          Y ~ Xpred + C1 + C2 + C3 + Upred,
          weights = (1 / df$pS),
          data = df
        )
      }
    } else if (len_c > 3) {
      stop("This function is currently not compatible with >3 confounders.")
    }
  } else if (!is.null(x1u0_model_coefs)) {
    len_x1u0_coefs <- length(x1u0_model_coefs)
    len_x0u1_coefs <- length(x0u1_model_coefs)
    len_x1u1_coefs <- length(x1u1_model_coefs)
    len_s_coefs <- length(s_model_coefs)

    if (len_x1u0_coefs != 3 + len_c) {
      stop(
        paste0(
          "Incorrect length of X1U0 model coefficients. ",
          "Length should equal 3 + number of confounders."
        )
      )
    }
    if (len_x0u1_coefs != 3 + len_c) {
      stop(
        paste0(
          "Incorrect length of X0U1 model coefficients. ",
          "Length should equal 3 + number of confounders."
        )
      )
    }
    if (len_x1u1_coefs != 3 + len_c) {
      stop(
        paste0(
          "Incorrect length of X1U1 model coefficients. ",
          "Length should equal 3 + number of confounders."
        )
      )
    }
    if (len_s_coefs != 3 + len_c) {
      stop(
        paste0(
          "Incorrect length of S model coefficients. ",
          "Length should equal 3 + number of confounders."
        )
      )
    }

    s1_0 <- s_model_coefs[1]
    s1_xstar <- s_model_coefs[2]
    s1_y <- s_model_coefs[3]

    x1u0_0 <- x1u0_model_coefs[1]
    x1u0_xstar <- x1u0_model_coefs[2]
    x1u0_y <- x1u0_model_coefs[3]

    x0u1_0 <- x0u1_model_coefs[1]
    x0u1_xstar <- x0u1_model_coefs[2]
    x0u1_y <- x0u1_model_coefs[3]

    x1u1_0 <- x1u1_model_coefs[1]
    x1u1_xstar <- x1u1_model_coefs[2]
    x1u1_y <- x1u1_model_coefs[3]

    if (is.null(confounders)) {
      df <- data.frame(Xstar = xstar, Y = y)

      p_x1u0 <- exp(x1u0_0 + x1u0_xstar * df$Xstar + x1u0_y * df$Y)
      p_x0u1 <- exp(x0u1_0 + x0u1_xstar * df$Xstar + x0u1_y * df$Y)
      p_x1u1 <- exp(x1u1_0 + x1u1_xstar * df$Xstar + x1u1_y * df$Y)

      denom <- (1 + p_x1u0 + p_x0u1 + p_x1u1)

      x0u0_pred <- 1 / denom
      x1u0_pred <- p_x1u0 / denom
      x0u1_pred <- p_x0u1 / denom
      x1u1_pred <- p_x1u1 / denom

      df_xu_pred <- data.frame(
        X0U0 = x0u0_pred,
        X1U0 = x1u0_pred,
        X0U1 = x0u1_pred,
        X1U1 = x1u1_pred
      )
      df_xu_pred4 <- bind_rows(df_xu_pred, df_xu_pred, df_xu_pred, df_xu_pred)

      combined <- bind_rows(df, df, df, df) %>%
        bind_cols(df_xu_pred4) %>%
        mutate(
          Xbar = rep(c(1, 0, 1, 0), each = n),
          Ubar = rep(c(1, 1, 0, 0), each = n),
          pXU = case_when(
            Xbar == 0 & Ubar == 0 ~ X0U0,
            Xbar == 1 & Ubar == 0 ~ X1U0,
            Xbar == 0 & Ubar == 1 ~ X0U1,
            Xbar == 1 & Ubar == 1 ~ X1U1
          ),
          pS = plogis(s1_0 + s1_xstar * .data$Xstar + s1_y * .data$Y)
        )

      if (y_binary) {
        suppressWarnings({
          final <- glm(
            Y ~ Xbar + Ubar,
            family = binomial(link = "logit"),
            weights = (combined$pXU / combined$pS),
            data = combined
          )
        })
      } else {
        final <- lm(
          Y ~ Xbar + Ubar,
          weights = (combined$pXU / combined$pS),
          data = combined
        )
      }
    } else if (len_c == 1) {
      c1 <- data[, confounders]
      df <- data.frame(Xstar = xstar, Y = y, C1 = c1)

      s1_c1 <- s_model_coefs[4]
      x1u0_c1 <- x1u0_model_coefs[4]
      x0u1_c1 <- x0u1_model_coefs[4]
      x1u1_c1 <- x1u1_model_coefs[4]

      p_x1u0 <- exp(
        x1u0_0 + x1u0_xstar * df$Xstar + x1u0_y * df$Y +
          x1u0_c1 * df$C1
      )
      p_x0u1 <- exp(
        x0u1_0 + x0u1_xstar * df$Xstar + x0u1_y * df$Y +
          x0u1_c1 * df$C1
      )
      p_x1u1 <- exp(
        x1u1_0 + x1u1_xstar * df$Xstar + x1u1_y * df$Y +
          x1u1_c1 * df$C1
      )

      denom <- (1 + p_x1u0 + p_x0u1 + p_x1u1)

      x0u0_pred <- 1 / denom
      x1u0_pred <- p_x1u0 / denom
      x0u1_pred <- p_x0u1 / denom
      x1u1_pred <- p_x1u1 / denom

      df_xu_pred <- data.frame(
        X0U0 = x0u0_pred,
        X1U0 = x1u0_pred,
        X0U1 = x0u1_pred,
        X1U1 = x1u1_pred
      )
      df_xu_pred4 <- bind_rows(df_xu_pred, df_xu_pred, df_xu_pred, df_xu_pred)

      combined <- bind_rows(df, df, df, df) %>%
        bind_cols(df_xu_pred4) %>%
        mutate(
          Xbar = rep(c(1, 0, 1, 0), each = n),
          Ubar = rep(c(1, 1, 0, 0), each = n),
          pXU = case_when(
            Xbar == 0 & Ubar == 0 ~ X0U0,
            Xbar == 1 & Ubar == 0 ~ X1U0,
            Xbar == 0 & Ubar == 1 ~ X0U1,
            Xbar == 1 & Ubar == 1 ~ X1U1
          ),
          pS = plogis(
            s1_0 + s1_xstar * .data$Xstar + s1_y * .data$Y +
              s1_c1 * .data$C1
          )
        )

      if (y_binary) {
        suppressWarnings({
          final <- glm(
            Y ~ Xbar + C1 + Ubar,
            family = binomial(link = "logit"),
            weights = (combined$pXU / combined$pS),
            data = combined
          )
        })
      } else {
        final <- lm(
          Y ~ Xbar + C1 + Ubar,
          weights = (combined$pXU / combined$pS),
          data = combined
        )
      }
    } else if (len_c == 2) {
      c1 <- data[, confounders[1]]
      c2 <- data[, confounders[2]]

      df <- data.frame(Xstar = xstar, Y = y, C1 = c1, C2 = c2)

      s1_c1 <- s_model_coefs[4]
      s1_c2 <- s_model_coefs[5]

      x1u0_c1 <- x1u0_model_coefs[4]
      x1u0_c2 <- x1u0_model_coefs[5]

      x0u1_c1 <- x0u1_model_coefs[4]
      x0u1_c2 <- x0u1_model_coefs[5]

      x1u1_c1 <- x1u1_model_coefs[4]
      x1u1_c2 <- x1u1_model_coefs[5]

      p_x1u0 <- exp(
        x1u0_0 + x1u0_xstar * df$Xstar + x1u0_y * df$Y +
          x1u0_c1 * df$C1 + x1u0_c2 * df$C2
      )
      p_x0u1 <- exp(
        x0u1_0 + x0u1_xstar * df$Xstar + x0u1_y * df$Y +
          x0u1_c1 * df$C1 + x0u1_c2 * df$C2
      )
      p_x1u1 <- exp(
        x1u1_0 + x1u1_xstar * df$Xstar + x1u1_y * df$Y +
          x1u1_c1 * df$C1 + x1u1_c2 * df$C2
      )

      denom <- (1 + p_x1u0 + p_x0u1 + p_x1u1)

      x0u0_pred <- 1 / denom
      x1u0_pred <- p_x1u0 / denom
      x0u1_pred <- p_x0u1 / denom
      x1u1_pred <- p_x1u1 / denom

      df_xu_pred <- data.frame(
        X0U0 = x0u0_pred,
        X1U0 = x1u0_pred,
        X0U1 = x0u1_pred,
        X1U1 = x1u1_pred
      )
      df_xu_pred4 <- bind_rows(df_xu_pred, df_xu_pred, df_xu_pred, df_xu_pred)

      combined <- bind_rows(df, df, df, df) %>%
        bind_cols(df_xu_pred4) %>%
        mutate(
          Xbar = rep(c(1, 0, 1, 0), each = n),
          Ubar = rep(c(1, 1, 0, 0), each = n),
          pXU = case_when(
            Xbar == 0 & Ubar == 0 ~ X0U0,
            Xbar == 1 & Ubar == 0 ~ X1U0,
            Xbar == 0 & Ubar == 1 ~ X0U1,
            Xbar == 1 & Ubar == 1 ~ X1U1
          ),
          pS = plogis(
            s1_0 + s1_xstar * .data$Xstar + s1_y * .data$Y +
              s1_c1 * .data$C1 + s1_c2 * .data$C2
          )
        )

      if (y_binary) {
        suppressWarnings({
          final <- glm(
            Y ~ Xbar + C1 + C2 + Ubar,
            family = binomial(link = "logit"),
            weights = (combined$pXU / combined$pS),
            data = combined
          )
        })
      } else {
        final <- lm(
          Y ~ Xbar + C1 + C2 + Ubar,
          weights = (combined$pXU / combined$pS),
          data = combined
        )
      }
    } else if (len_c == 3) {
      c1 <- data[, confounders[1]]
      c2 <- data[, confounders[2]]
      c3 <- data[, confounders[3]]

      df <- data.frame(Xstar = xstar, Y = y, C1 = c1, C2 = c2, C3 = c3)

      s1_c1 <- s_model_coefs[4]
      s1_c2 <- s_model_coefs[5]
      s1_c3 <- s_model_coefs[6]

      x1u0_c1 <- x1u0_model_coefs[4]
      x1u0_c2 <- x1u0_model_coefs[5]
      x1u0_c3 <- x1u0_model_coefs[6]

      x0u1_c1 <- x0u1_model_coefs[4]
      x0u1_c2 <- x0u1_model_coefs[5]
      x0u1_c3 <- x0u1_model_coefs[6]

      x1u1_c1 <- x1u1_model_coefs[4]
      x1u1_c2 <- x1u1_model_coefs[5]
      x1u1_c3 <- x1u1_model_coefs[6]

      p_x1u0 <- exp(
        x1u0_0 + x1u0_xstar * df$Xstar + x1u0_y * df$Y +
          x1u0_c1 * df$C1 + x1u0_c2 * df$C2 + x1u0_c3 * df$C3
      )
      p_x0u1 <- exp(
        x0u1_0 + x0u1_xstar * df$Xstar + x0u1_y * df$Y +
          x0u1_c1 * df$C1 + x0u1_c2 * df$C2 + x0u1_c3 * df$C3
      )
      p_x1u1 <- exp(
        x1u1_0 + x1u1_xstar * df$Xstar + x1u1_y * df$Y +
          x1u1_c1 * df$C1 + x1u1_c2 * df$C2 + x1u1_c3 * df$C3
      )

      denom <- (1 + p_x1u0 + p_x0u1 + p_x1u1)

      x0u0_pred <- 1 / denom
      x1u0_pred <- p_x1u0 / denom
      x0u1_pred <- p_x0u1 / denom
      x1u1_pred <- p_x1u1 / denom

      df_xu_pred <- data.frame(
        X0U0 = x0u0_pred,
        X1U0 = x1u0_pred,
        X0U1 = x0u1_pred,
        X1U1 = x1u1_pred
      )
      df_xu_pred4 <- bind_rows(df_xu_pred, df_xu_pred, df_xu_pred, df_xu_pred)

      combined <- bind_rows(df, df, df, df) %>%
        bind_cols(df_xu_pred4) %>%
        mutate(
          Xbar = rep(c(1, 0, 1, 0), each = n),
          Ubar = rep(c(1, 1, 0, 0), each = n),
          pXU = case_when(
            Xbar == 0 & Ubar == 0 ~ X0U0,
            Xbar == 1 & Ubar == 0 ~ X1U0,
            Xbar == 0 & Ubar == 1 ~ X0U1,
            Xbar == 1 & Ubar == 1 ~ X1U1
          ),
          pS = plogis(
            s1_0 + s1_xstar * .data$Xstar + s1_y * .data$Y +
              s1_c1 * .data$C1 + s1_c2 * .data$C2 +
              s1_c3 * .data$C3
          )
        )

      if (y_binary) {
        suppressWarnings({
          final <- glm(
            Y ~ Xbar + C1 + C2 + C3 + Ubar,
            family = binomial(link = "logit"),
            weights = (combined$pXU / combined$pS),
            data = combined
          )
        })
      } else {
        final <- lm(
          Y ~ Xbar + C1 + C2 + C3 + Ubar,
          weights = (combined$pXU / combined$pS),
          data = combined
        )
      }
    } else if (len_c > 3) {
      stop("This function is currently not compatible with >3 confounders.")
    }
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
