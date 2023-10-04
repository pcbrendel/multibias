#' Simulated data with outcome misclassification
#'
#' Data containing one source of bias, a known confounder, and
#'  100,000 observations. This data is obtained from \code{df_omc_source}
#'  by removing the column Y. The resulting data corresponds to
#'  what a researcher would see in the real-world: a misclassified outcome,
#'  Ystar, and no data on the true outcome. As seen in
#'  \code{df_omc_source}, the true, unbiased exposure-outcome odds ratio = 2.
#'
#' @format A dataframe with 100,000 rows and 3 columns:
#' \describe{
#'     \item{X}{exposure, 1 = present and 0 = absent}
#'     \item{Ystar}{misclassified outcome, 1 = present and 0 = absent}
#'     \item{C1}{1st confounder, 1 = present and 0 = absent}
#' }
"df_omc"