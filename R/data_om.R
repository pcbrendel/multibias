#' Simulated data with outcome misclassification
#'
#' Data containing one source of bias, three known confounders, and
#' 100,000 observations. This data is obtained from \code{df_om_source}
#' by removing the column \emph{Y}. The resulting data corresponds to
#' what a researcher would see in the real-world: a misclassified outcome,
#' \emph{Ystar}, and no data on the true outcome. As seen in
#' \code{df_om_source}, the true, unbiased exposure-outcome odds ratio = 2.
#'
#' @format A dataframe with 100,000 rows and 5 columns:
#' \describe{
#'     \item{X}{exposure, 1 = present and 0 = absent}
#'     \item{Ystar}{misclassified outcome, 1 = present and 0 = absent}
#'     \item{C1}{1st confounder, 1 = present and 0 = absent}
#'     \item{C2}{2nd confounder, 1 = present and 0 = absent}
#'     \item{C3}{3rd confounder, 1 = present and 0 = absent}
#' }
"df_om"

#' Data source for \code{df_om}
#'
#' Data with complete information on one sources of bias, three known
#' confounders, and 100,000 observations. This data is used to derive
#' \code{df_om} and can be used to obtain bias parameters for purposes
#' of validating the simultaneous multi-bias adjustment method with
#' \code{df_om}. With this source data, the fitted regression
#' \ifelse{html}{\out{logit(P(Y=1)) = &alpha;<sub>0</sub> + &alpha;<sub>1</sub>X + &alpha;<sub>2</sub>C1 + &alpha;<sub>3</sub>C2 + &alpha;<sub>4</sub>C3}}{\eqn{logit(P(Y=1)) = \alpha_0 + \alpha_1 X + \alpha_2 C1 + \alpha_3 C2 + \alpha_4 C3}}
#' shows that the true, unbiased exposure-outcome odds ratio = 2.
#'
#' @format A dataframe with 100,000 rows and 6 columns:
#' \describe{
#'     \item{X}{exposure, 1 = present and 0 = absent}
#'     \item{Y}{true outcome, 1 = present and 0 = absent}
#'     \item{C1}{1st confounder, 1 = present and 0 = absent}
#'     \item{C2}{2nd confounder, 1 = present and 0 = absent}
#'     \item{C3}{3rd confounder, 1 = present and 0 = absent}
#'     \item{Ystar}{misclassified outcome, 1 = present and 0 = absent}
#' }
"df_om_source"