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

  # adjusted_est <- multibias_result$estimate
  # adjusted_ci_low <- multibias_result$ci[1]
  # adjusted_ci_high <- multibias_result$ci[2]

  # print(observed_est)
  # print(observed_se)
  # print(adjusted_est)
  # print(adjusted_se)

  # xmin <- min(observed_est - observed_se * 5, adjusted_est - adjusted_se * 5)
  # xmax <- max(observed_est + observed_se * 5, adjusted_est + adjusted_se * 5)

  # df <- data.frame(
  #   type = c("Observed", "Adjusted"),
  #   est = c(observed_est, adjusted_est),
  #   ci_low = c(observed_ci_low, adjusted_ci_low),
  #   ci_high = c(observed_ci_high, adjusted_ci_high)
  # )

  list_names <- names(multibias_result_list)

  df_adjusted <- purrr::map_dfr(seq_along(multibias_result_list), function(list_index) {
    list_item <- multibias_result_list[[list_index]]
    name <- list_names[list_index]

    is_valid <- is.list(list_item) &&
      all(c("estimate", "ci") %in% names(list_item)) &&
      length(list_item$ci) == 2 &&
      is.numeric(list_item$estimate) &&
      is.numeric(list_item$ci)

    if (!is_valid) {
      warning(paste("Item", name, "in multibias_result_list is not structured as expected. Skipping this item."), call. = FALSE)
      return(NULL)
    }
    tibble(
      type = name,
      est = list_item$estimate,
      ci_low  = list_item$ci[1],
      ci_high = list_item$ci[2]
    )
  })

  df <- rbind(df_observed, as.data.frame(df_adjusted))

  final <- ggplot(data = df, aes(x = est, y = type)) +
    geom_vline(aes(xintercept = 1), linetype = "dashed") +
    geom_point(color = "blue") +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), size = .5, height = .2, color = "gray50") +
    labs(
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
    theme_bw()

  if (log_scale) {
    if(any(df$est <= 0 | df$ci_low <= 0, na.rm = TRUE)) {
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
      scale_x_log10()
  }

  return(final)
}
