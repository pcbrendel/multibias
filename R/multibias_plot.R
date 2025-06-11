#' Create a Forest Plot comparing observed and adjusted effect estimates
#'
#' This function generates a forest plot comparing the observed effect estimate
#' with adjusted effect estimates from sensitivity analyses. The plot includes
#' point estimates and confidence intervals for each analysis.
#'
#' @param data_observed Object of class `data_observed` representing the
#' observed causal data and effect of interest.
#' @param multibias_result_list A named list of sensitivity analysis results.
#' Each element should be a result from [multibias_adjust()].
#' @param log_scale Boolean indicating whether to display the x-axis on the
#' log scale. Default is FALSE.
#'
#' @return A ggplot object showing a forest plot with:
#'   \itemize{
#'     \item Point estimates (blue dots)
#'     \item Confidence intervals (gray horizontal lines)
#'     \item A vertical reference line at x=1 (dashed)
#'     \item Appropriate labels and title
#'   }
#'
#' @examples
#' df_observed <- data_observed(
#'   data = df_em,
#'   bias = "em",
#'   exposure = "Xstar",
#'   outcome = "Y",
#'   confounders = "C1"
#' )
#'
#' bp1 <- bias_params(coef_list = list(x = c(-2.10, 1.62, 0.63, 0.35)))
#' bp2 <- bias_params(coef_list = list(x = c(-2.10 * 2, 1.62 * 2, 0.63 * 2, 0.35 * 2)))
#'
#' result1 <- multibias_adjust(
#'   data_observed = df_observed,
#'   bias_params = bp1
#' )
#' result2 <- multibias_adjust(
#'   data_observed = df_observed,
#'   bias_params = bp2
#' )
#'
#' multibias_plot(
#'   data_observed = df_observed,
#'   multibias_result_list = list(
#'     "Adjusted with bias params" = result1,
#'     "Adjusted with bias params doubled" = result2
#'   )
#' )
#'
#' @export

multibias_plot <- function(
    data_observed,
    multibias_result_list,
    log_scale = FALSE) {
  summary_quietly <- purrr::quietly(summary)
  observed_summary <- summary_quietly(data_observed)$result

  observed_est <- observed_summary$estimate[2]
  observed_ci_low <- observed_summary$conf.low[2]
  observed_ci_high <- observed_summary$conf.high[2]

  df_observed <- tibble(
    type = "Observed",
    est = observed_est,
    ci_low = observed_ci_low,
    ci_high = observed_ci_high
  )

  list_names <- names(multibias_result_list)

  list_adjusted <- purrr::map(
    seq_along(multibias_result_list),
    function(list_index) {
      list_item <- multibias_result_list[[list_index]]
      name <- list_names[list_index]

      is_valid <- is.list(list_item) &&
        all(c("estimate", "ci") %in% names(list_item)) &&
        length(list_item$ci) == 2 &&
        is.numeric(list_item$estimate) &&
        is.numeric(list_item$ci)

      if (!is_valid) {
        warning(
          paste0(
            "Item ",
            name,
            " in multibias_result_list is not structured as expected. ",
            "Skipping this item."
          ),
          call. = FALSE
        )
        return(NULL)
      }
      tibble(
        type = name,
        est = list_item$estimate,
        ci_low = list_item$ci[1],
        ci_high = list_item$ci[2]
      )
    }
  )

  df_adjusted <- purrr::list_rbind(list_adjusted)
  df <- rbind(df_observed, as.data.frame(df_adjusted))

  # Set factor levels to control y-axis order
  df$type <- factor(df$type, levels = rev(c("Observed", list_names)))

  final <- ggplot2::ggplot(
    data = df, ggplot2::aes(x = .data$est, y = .data$type)
  ) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = 1), linetype = "dashed") +
    ggplot2::geom_point(color = "blue") +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = .data$ci_low, xmax = .data$ci_high),
      size = .5,
      height = .2,
      color = "gray50"
    ) +
    ggplot2::labs(
      title = "Observed vs. Adjusted Effect Estimates",
      subtitle = paste0(
        "Sensitivity analysis of the ",
        data_observed$exposure,
        "-",
        data_observed$outcome,
        " relationship"
      ),
      x = "Effect Estimate",
      y = ""
    ) +
    ggplot2::theme_bw()

  if (log_scale) {
    if (any(df$est <= 0 | df$ci_low <= 0, na.rm = TRUE)) {
      warning(
        paste0(
          "Some estimates or lower CI bounds are zero or negative. ",
          "Log scale may produce NaNs, errors, or missing data points. ",
          "Review data or consider if log scale is appropriate."
        ),
        call. = FALSE
      )
    }
    final <- final +
      ggplot2::scale_x_log10()
  }

  return(final)
}
