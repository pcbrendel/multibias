#' Simulated data with uncontrolled confounding and exposure misclassification
#'
#' Data containing two sources of bias, a known confounders, and
#'  100,000 observations. This data is obtained from \code{df_uc_emc_source}
#'  then removing the columns X and U. The resulting data corresponds to
#'  what a researcher would see in the real-world: a misclassified exposure,
#'  Xstar, and missing data on a confounder U. As seen in
#'  \code{df_uc_emc_source}, the true, unbiased exposure-outcome odds ratio = 2.
#'
#' @format A dataframe with 100,000 rows and 3 columns:
#' \describe{
#'     \item{Xstar}{misclassified exposure, 1 = present and 0 = absent}
#'     \item{Y}{outcome, 1 = present and 0 = absent}
#'     \item{C1}{confounder, 1 = present and 0 = absent}
#' }
"df_uc_emc"