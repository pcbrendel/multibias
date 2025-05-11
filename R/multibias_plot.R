multibias_plot <- function(
  data_observed,
  multibias_result
) {
  observed_est <- unname(summary(df_observed)[2, 2])
  observed_se <- unname(summary(df_observed)[2, 3])
  adjusted_est <- multibias_result$estimate
  adjusted_se <- multibias_result$std.error

  final <- x

  return(final)
}