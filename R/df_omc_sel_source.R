#' Data source for df_omc_sel
#'
#' Data with complete information on the two sources of bias, a known
#'  confounder, and 100,000 observations. This data is used to derive
#'  \code{df_omc_sel} and can be used to obtain bias parameters for purposes
#'  of validating the simultaneous multi-bias adjustment method with
#'  \code{df_omc_sel}. With this source data, the fitted regression
#'  \ifelse{html}{\out{logit(P(Y=1)) = &alpha;<sub>0</sub> + &alpha;<sub>1</sub>X + &alpha;<sub>2</sub>C1}}{\eqn{logit(P(Y=1)) = \alpha_0 + \alpha_1 X + \alpha_2 C1}}
#'  shows that the true, unbiased exposure-outcome odds ratio = 2.
#'
#' @format A dataframe with 100,000 rows and 5 columns:
#' \describe{
#'     \item{X}{exposure, 1 = present and 0 = absent}
#'     \item{Y}{true outcome, 1 = present and 0 = absent}
#'     \item{C1}{1st confounder, 1 = present and 0 = absent}
#'     \item{Ystar}{misclassified outcome, 1 = present and 0 = absent}
#'     \item{S}{selection, 1 = selected into the study and 0 = not selected into the study}
#' }
"df_omc_sel_source"