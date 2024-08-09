#' Simulated data with uncontrolled confounding and exposure misclassification
#'
#' Data containing two sources of bias, three known confounders, and
#' 100,000 observations. This data is obtained from \code{df_uc_emc_source}
#' by removing the columns \emph{X} and \emph{U}. The resulting data
#' corresponds to what a researcher would see in the real-world: a
#' misclassified exposure, \emph{Xstar}, and missing data on a confounder
#' \emph{U}. As seen in \code{df_uc_emc_source}, the true, unbiased
#' exposure-outcome odds ratio = 2.
#'
#' @format A dataframe with 100,000 rows and 5 columns:
#' \describe{
#'     \item{Xstar}{misclassified exposure, 1 = present and 0 = absent}
#'     \item{Y}{outcome, 1 = present and 0 = absent}
#'     \item{C1}{1st confounder, 1 = present and 0 = absent}
#'     \item{C2}{2nd confounder, 1 = present and 0 = absent}
#'     \item{C3}{3rd confounder, 1 = present and 0 = absent}
#' }
"df_uc_emc"