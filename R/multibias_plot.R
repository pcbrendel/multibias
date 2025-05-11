multibias_plot <- function(
    data_observed,
    multibias_result) {
  summary_quietly <- quietly(summary)
  observed_summary <- summary_quietly(data_observed)$result

  observed_est <- observed_summary[[2]][2]
  observed_se <- observed_summary[[3]][2]
  adjusted_est <- multibias_result$estimate
  adjusted_se <- multibias_result$std.error

  xmin <- min(observed_est - observed_se * 4, adjusted_est - adjusted_se * 4)
  xmax <- max(observed_est + observed_se * 4, adjusted_est + adjusted_se * 4)

  final <- ggplot() +
    xlim(xmin, xmax) +
    stat_function(
      fun = dnorm,
      args = list(mean = observed_est, sd = observed_se),
      aes(color = "Observed")
    ) +
    stat_function(
      fun = dnorm,
      args = list(mean = adjusted_est, sd = adjusted_se),
      aes(color = "Adjusted")
    ) +
    labs(
      title = "Observed vs. Adjusted Estimates",
      x = "Effect Estimate",
      y = "Density",
    ) +
    theme_bw()

  return(final)
}
