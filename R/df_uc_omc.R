#' Simulated data with uncontrolled confounding and outcome misclassification
#'
#' Data containing two sources of bias, a known confounder, and
#'  100,000 observations. This data is obtained from \code{df_uc_omc_source}
#'  by removing the columns Y and U. The resulting data corresponds to
#'  what a researcher would see in the real-world: a misclassified outcome,
#'  Ystar, and missing data on a confounder U. As seen in
#'  \code{df_uc_omc_source}, the true, unbiased exposure-outcome odds ratio = 2.
#'
#' @format A dataframe with 100,000 rows and 3 columns:
#' \describe{
#'     \item{X}{exposure, 1 = present and 0 = absent}
#'     \item{Ystar}{misclassified outcome, 1 = present and 0 = absent}
#'     \item{C1}{confounder, 1 = present and 0 = absent}
#' }
"df_uc_omc"