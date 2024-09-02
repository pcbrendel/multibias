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

#' Data source for \code{df_em_om}
#'
#' Data with complete information on the two sources of bias, three known
#' confounders, and 100,000 observations. This data is used to derive
#' \code{df_em_om} and can be used to obtain bias parameters for purposes
#' of validating the simultaneous multi-bias adjustment method with
#' \code{df_em_om}. With this source data, the fitted regression
#' \ifelse{html}{\out{logit(P(Y=1)) = &alpha;<sub>0</sub> + &alpha;<sub>1</sub>X + &alpha;<sub>2</sub>C1 + &alpha;<sub>3</sub>C2 + &alpha;<sub>4</sub>C3}}{\eqn{logit(P(Y=1)) = \alpha_0 + \alpha_1 X + \alpha_2 C1 + \alpha_3 C2 + \alpha_4 C3}}
#' shows that the true, unbiased exposure-outcome odds ratio = 2.
#'
#' @format A dataframe with 100,000 rows and 7 columns:
#' \describe{
#'     \item{X}{true exposure, 1 = present and 0 = absent}
#'     \item{Y}{outcome, 1 = present and 0 = absent}
#'     \item{C1}{1st confounder, 1 = present and 0 = absent}
#'     \item{C2}{2nd confounder, 1 = present and 0 = absent}
#'     \item{C3}{3rd confounder, 1 = present and 0 = absent}
#'     \item{Xstar}{misclassified exposure, 1 = present and 0 = absent}
#'     \item{Ystar}{misclassified outcome, 1 = present and 0 = absent}
#' }
"df_em_om_source"