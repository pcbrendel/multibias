#' Simulated data with uncontrolled confounding, exposure misclassification,
#' and selection bias
#'
#' Data containing three sources of bias, three known confounders, and
#' 100,000 observations. This data is obtained by sampling with replacement
#' with probability = \emph{S} from \code{df_uc_emc_sel_source} then removing
#' the columns \emph{X}, \emph{U}, and \emph{S}. The resulting data corresponds
#' to what a researcher would see in the real-world: a misclassified exposure,
#' \emph{Xstar}; missing data on a confounder \emph{U}; and missing data for
#' those not selected into the study (\emph{S}=0). As seen in
#' \code{df_uc_emc_sel_source}, the true, unbiased exposure-outcome
#' odds ratio = 2.
#'
#' @format A dataframe with 100,000 rows and 5 columns:
#' \describe{
#'     \item{Xstar}{misclassified exposure, 1 = present and 0 = absent}
#'     \item{Y}{outcome, 1 = present and 0 = absent}
#'     \item{C1}{1st confounder, 1 = present and 0 = absent}
#'     \item{C2}{2nd confounder, 1 = present and 0 = absent}
#'     \item{C3}{3rd confounder, 1 = present and 0 = absent}
#' }
"df_uc_emc_sel"