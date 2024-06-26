#' Data source for df_uc_omc_sel
#'
#' Data with complete information on the three sources of bias, three known
#'  confounders, and 100,000 observations. This data is used to derive
#'  \code{df_uc_omc_sel} and can be used to obtain bias parameters for purposes
#'  of validating the simultaneous multi-bias adjustment method with
#'  \code{df_uc_omc_sel}. With this source data, the fitted regression
#'  \ifelse{html}{\out{logit(P(Y=1)) = &alpha;<sub>0</sub> + &alpha;<sub>1</sub>X + &alpha;<sub>2</sub>C1 + &alpha;<sub>3</sub>C2 + &alpha;<sub>4</sub>C3 + &alpha;<sub>5</sub>U}}{\eqn{logit(P(Y=1)) = \alpha_0 + \alpha_1 X + \alpha_2 C1 + \alpha_3 C2 + \alpha_4 C3 + \alpha_5 U}}
#'  shows that the true, unbiased exposure-outcome odds ratio = 2.
#'
#' @format A dataframe with 100,000 rows and 8 columns:
#' \describe{
#'     \item{X}{exposure, 1 = present and 0 = absent}
#'     \item{Y}{true outcome, 1 = present and 0 = absent}
#'     \item{C1}{1st confounder, 1 = present and 0 = absent}
#'     \item{C2}{2nd confounder, 1 = present and 0 = absent}
#'     \item{C3}{3rd confounder, 1 = present and 0 = absent}
#'     \item{U}{unmeasured confounder, 1 = present and 0 = absent}
#'     \item{Ystar}{misclassified outcome, 1 = present and 0 = absent}
#'     \item{S}{selection, 1 = selected into the study and 0 = not selected into the study}
#' }
"df_uc_omc_sel_source"