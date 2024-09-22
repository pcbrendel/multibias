#' Simulated data with uncontrolled confounding
#'
#' Data containing one source of bias, three known confounders, and
#' 100,000 observations. This data is obtained from `df_uc_source`
#' by removing the column *U*. The resulting data corresponds to
#' what a researcher would see in the real-world: information on known
#' confounders (*C1*, *C2*, and *C3*), but not for
#' confounder *U*.
#' As seen in `df_uc_source`, the true, unbiased exposure-outcome
#' effect estimate = 2.
#'
#' @format A dataframe with 100,000 rows and 7 columns:
#' \describe{
#'     \item{X_bi}{binary exposure, 1 = present and 0 = absent}
#'     \item{X_cont}{continuous exposure}
#'     \item{Y_bi}{binary outcome corresponding to exposure *X_bi*, 1 = present and 0 = absent}
#'     \item{Y_cont}{continuous outcome corresponding to exposure *X_cont*}
#'     \item{C1}{1st confounder, 1 = present and 0 = absent}
#'     \item{C2}{2nd confounder, 1 = present and 0 = absent}
#'     \item{C3}{3rd confounder, 1 = present and 0 = absent}
#' }
"df_uc"

#' Data source for `df_uc`
#'
#' Data with complete information on one source of bias, three known
#' confounders, and 100,000 observations. This data is used to derive
#' `df_uc` and can be used to obtain bias parameters for purposes
#' of validating the simultaneous multi-bias adjustment method with
#' `df_uc`. With this source data, the fitted regression
#' \ifelse{html}{\out{g(P(Y=1)) = &alpha;<sub>0</sub> + &alpha;<sub>1</sub>X + &alpha;<sub>2</sub>C1 + &alpha;<sub>3</sub>C2 + &alpha;<sub>4</sub>C3 + &alpha;<sub>5</sub>U}}{\eqn{logit(P(Y=1)) = \alpha_0 + \alpha_1 X + \alpha_2 C1 + \alpha_3 C2 + \alpha_4 C3 + \alpha_5 U}}
#' shows that the true, unbiased exposure-outcome effect estimate = 2 when:
#' \enumerate{
#'   \item g = logit, Y = *Y_bi*, and X = *X_bi* or
#'   \item g = identity, Y = *Y_cont*, X = *X_cont*.
#' }
#'
#' @format A dataframe with 100,000 rows and 8 columns:
#' \describe{
#'     \item{X_bi}{binary exposure, 1 = present and 0 = absent}
#'     \item{X_cont}{continuous exposure}
#'     \item{Y_bi}{binary outcome corresponding to exposure *X_bi*, 1 = present and 0 = absent}
#'     \item{Y_cont}{continuous outcome corresponding to exposure *X_cont*}
#'     \item{C1}{1st confounder, 1 = present and 0 = absent}
#'     \item{C2}{2nd confounder, 1 = present and 0 = absent}
#'     \item{C3}{3rd confounder, 1 = present and 0 = absent}
#'     \item{U}{uncontrolled confounder, 1 = present and 0 = absent}
#' }
"df_uc_source"
