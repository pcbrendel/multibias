#' Data source for df_uc_omc
#'
#' Data with complete information on the two sources of bias, a known
#'  confounder, and 100,000 observations. This data is used to derive
#'  \code{df_uc_omc} and can be used to obtain bias parameters for purposes
#'  of validating the simultaneous multi-bias adjustment method with
#'  \code{df_uc_omc}. The regression \ifelse{html}{\out{logit(P(Y=1)) =
#'  &alpha;<sub>0</sub> + &alpha;<sub>1</sub>X + &alpha;<sub>2</sub>C1 +
#'  &alpha;<sub>5</sub>U} shows that the true, unbiased exposure-outcome odds
#'  ratio = 2.}{\eqn{logit(P(Y=1)) =}}
#'
#' @format A dataframe with 100,000 rows and 5 columns:
#' \describe{
#'     \item{X}{true exposure, 1 = present and 0 = absent}
#'     \item{Y}{outcome, 1 = present and 0 = absent}
#'     \item{C1}{1st confounder, 1 = present and 0 = absent}
#'     \item{U}{unmeasured confounder, 1 = present and 0 = absent}
#'     \item{Ystar}{misclassified outcome, 1 = present and 0 = absent}
#' }
"df_uc_omc_source"