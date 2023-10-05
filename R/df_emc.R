#' Simulated data with exposure misclassification
#'
#' Data containing one source of bias, a known confounder, and
#'  100,000 observations. This data is obtained from \code{df_emc_source}
#'  by removing the column X. The resulting data corresponds to
#'  what a researcher would see in the real-world: a misclassified exposure,
#'  Xstar, and no data on the true exposure. As seen in
#'  \code{df_emc_source}, the true, unbiased exposure-outcome odds ratio = 2.
#'
#' @format A dataframe with 100,000 rows and 3 columns:
#' \describe{
#'     \item{Xstar}{misclassified exposure, 1 = present and 0 = absent}
#'     \item{Y}{outcome, 1 = present and 0 = absent}
#'     \item{C1}{1st confounder, 1 = present and 0 = absent}
#' }
"df_emc"