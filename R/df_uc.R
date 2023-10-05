#' Simulated data with uncontrolled confounding
#'
#' Data containing one source of bias, a known confounder, and
#'  100,000 observations. This data is obtained from \code{df_uc_source}
#'  by removing the column U. The resulting data corresponds to
#'  what a researcher would see in the real-world: a known confounder,
#'  C1, but no data on the other confounder, U. As seen in
#'  \code{df_uc_source}, the true, unbiased exposure-outcome odds ratio = 2.
#'
#' @format A dataframe with 100,000 rows and 3 columns:
#' \describe{
#'     \item{X}{exposure, 1 = present and 0 = absent}
#'     \item{Y}{outcome, 1 = present and 0 = absent}
#'     \item{C1}{1st confounder, 1 = present and 0 = absent}
#' }
"df_uc"