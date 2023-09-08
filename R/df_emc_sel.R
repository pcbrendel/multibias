#' Simulated data with exposure misclassification and selection bias
#'
#' Data containing two sources of bias, a known confounder, and
#'  100,000 observations. This data is obtained by sampling with replacement
#'  with probability = S from \code{df_emc_sel_source} then removing the
#'  columns X and S. The resulting data corresponds to what a researcher
#'  would see in the real-world: a misclassified exposure, Xstar,
#'  and missing data for those not selected into the study
#'  (S=0). As seen in \code{df_emc_sel_source}, the true, unbiased
#'  exposure-outcome odds ratio = 2.
#'
#' @format A dataframe with 100,000 rows and 3 columns:
#' \describe{
#'     \item{Xstar}{misclassified exposure, 1 = present and 0 = absent}
#'     \item{Y}{outcome, 1 = present and 0 = absent}
#'     \item{C1}{1st confounder, 1 = present and 0 = absent}
#' }
"df_emc_sel"
