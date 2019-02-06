#' Simulated dataset with uncontrolled confounding, exposure misclassification, 
#' and selection bias.
#'
#' A dataset with three sources of bias, three known confounders, and 100,000 observations.
#' The true, unbiased exposure-outcome odds ratio = 2.
#'
#' @format A dataframe with 500,000 rows and 5 variables:
#' \describe{
#'     \item{Xstar}{misclassified exposure, 1 = present and 0 = absent}
#'     \item{Y}{outcome, 1 = present and 0 = absent}
#'     \item{C1}{1st confounder, 1 = present and 0 = absent}
#'     \item{C2}{2nd confounder, 1 = present and 0 = absent}
#'     \item{C3}{3rd confounder, 1 = present and 0 = absent}
#' }
"df_uc_mc_sel"