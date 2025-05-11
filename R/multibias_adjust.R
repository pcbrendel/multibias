#' Simultaneously adjust for multiple biases
#'
#' `multibias_adjust` returns the exposure-outcome odds ratio and confidence
#' interval, adjusted for one or more biases.
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
#' the validation data used to adjust for bias in the observed data. The
#' validation data should have data for the same variables as in
#' `data_observed`, plus data for the missing variables leading to bias.
#' @param bias_params Object of class 'bias_params' corresponding to the
#' bias parameters used to adjust for bias in the observed data. There must
#' be parameters corresponding to the bias or biases specified in
#' `data_observed`.
#' @param bootstrap Boolean for whether to perform bootstrapping to obtain
#' the estimate and confidence interval.
#' @param bootstrap_reps Integer number of bootstrap samples to run in
#' bootstrapping.
#' @param level Value from 0-1 representing the full range of the confidence
#' interval. Default is 0.95.
#'
#' @return A list including: the bias-adjusted effect estimate of the exposure
#' on the outcome, the standard error, and the confidence interval as the
#' vector: (lower bound, upper bound).
#'
#' @examples
#' # Adjust for exposure misclassification -------------------------------------
#' df_observed <- data_observed(
#'   data = df_em,
#'   bias = "em",
#'   exposure = "Xstar",
#'   outcome = "Y",
#'   confounders = "C1"
#' )
#'
#' # Using validation data
#' df_validation <- data_validation(
#'   data = df_em_source,
#'   true_exposure = "X",
#'   true_outcome = "Y",
#'   confounders = "C1",
#'   misclassified_exposure = "Xstar"
#' )
#'
#' multibias_adjust(
#'   data_observed = df_observed,
#'   data_validation = df_validation
#' )
#'
#' # Using bias_params
#' bp <- bias_params(coef_list = list(x = c(-2.10, 1.62, 0.63, 0.35)))
#'
#' multibias_adjust(
#'   data_observed = df_observed,
#'   bias_params = bp
#' )
#'
#' # Adjust for three biases ---------------------------------------------------
#' df_observed <- data_observed(
#'   data = df_uc_om_sel,
#'   bias = c("uc", "om", "sel"),
#'   exposure = "X",
#'   outcome = "Ystar",
#'   confounders = c("C1", "C2", "C3")
#' )
#'
#' # Using validation data
#' df_validation <- data_validation(
#'   data = df_uc_om_sel_source,
#'   true_exposure = "X",
#'   true_outcome = "Y",
#'   confounders = c("C1", "C2", "C3", "U"),
#'   misclassified_outcome = "Ystar",
#'   selection = "S"
#' )
#'
#' multibias_adjust(
#'   data_observed = df_observed,
#'   data_validation = df_validation
#' )
#'
#' # Using bias_params
#' bp1 <- bias_params(
#'   coef_list = list(
#'     u = c(-0.32, 0.59, 0.69),
#'     y = c(-2.85, 0.71, 1.63, 0.40, -0.85, 0.22),
#'     s = c(0.00, 0.74, 0.19, 0.02, -0.06, 0.02)
#'   )
#' )
#'
#' multibias_adjust(
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
#' # with bootstrapping
#' multibias_adjust(
#'   data_observed = df_observed,
#'   bias_params = bp2,
#'   bootstrap = TRUE,
#'   bootstrap_reps = 10
#' )
#'
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom stats as.formula binomial coef lm median glm plogis rbinom
#' @importFrom stats qnorm quantile
#' @importFrom rlang .data
#'
#' @export

multibias_adjust <- function(
    data_observed,
    data_validation = NULL,
    bias_params = NULL,
    bootstrap = FALSE,
    bootstrap_reps = 100,
    level = 0.95) {
  check_inputs2(data_validation, bias_params)

  bias_dict <- list(
    "uc" = adjust_uc,
    "em" = adjust_em,
    "om" = adjust_om,
    "sel" = adjust_sel,
    "em_uc" = adjust_uc_em,
    "om_uc" = adjust_uc_om,
    "sel_uc" = adjust_uc_sel,
    "em_om" = adjust_em_om,
    "em_sel" = adjust_em_sel,
    "om_sel" = adjust_om_sel,
    "em_sel_uc" = adjust_uc_em_sel,
    "om_sel_uc" = adjust_uc_om_sel
  )

  bias_key <- paste(sort(data_observed$bias), collapse = "_")
  adjust_fn <- bias_dict[[bias_key]]

  if (bootstrap == FALSE) {
    output <- adjust_fn(data_observed, data_validation, bias_params, level)
  } else if (bootstrap == TRUE) {
    df <- data_observed$data
    n <- nrow(df)
    est <- numeric(bootstrap_reps)

    indices <- matrix(
      sample.int(n, n * bootstrap_reps, replace = TRUE),
      nrow = n,
      ncol = bootstrap_reps
    )

    for (i in seq_len(bootstrap_reps)) {
      data_observed$data <- df[indices[, i], ]
      final <- adjust_fn(
        data_observed,
        data_validation,
        bias_params
      )
      est[i] <- final$estimate
    }

    alpha <- 1 - level
    output <- list(
      estimate = median(est),
      std.error = sd(est),
      ci = as.vector(quantile(est, c(alpha / 2, 1 - alpha / 2)))
    )
  }

  return(output)
}
