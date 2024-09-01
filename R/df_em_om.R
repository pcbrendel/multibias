#' Simulated data with exposure misclassification and outcome misclassification
#'
#' Data containing two sources of bias, three known confounders, and
#' 100,000 observations. This data is obtained from \code{df_emc_omc_source}
#' by removing the columns \emph{X} and \emph{Y}. The resulting data corresponds
#' to what a researcher would see in the real-world: a misclassified exposure,
#' \emph{Xstar}, and a misclassified outcome, \emph{Ystar}. As seen in
#' \code{df_em_om_source}, the true, unbiased exposure-outcome
#' odds ratio = 2.
#'
#' @format A dataframe with 100,000 rows and 5 columns:
#' \describe{
#'     \item{Xstar}{misclassified exposure, 1 = present and 0 = absent}
#'     \item{Ystar}{misclassified outcome, 1 = present and 0 = absent}
#'     \item{C1}{1st confounder, 1 = present and 0 = absent}
#'     \item{C2}{2nd confounder, 1 = present and 0 = absent}
#'     \item{C3}{3rd confounder, 1 = present and 0 = absent}
#' }
"df_em_om"