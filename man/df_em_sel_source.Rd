% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_em_sel.R
\docType{data}
\name{df_em_sel_source}
\alias{df_em_sel_source}
\title{Data source for \code{df_em_sel}}
\format{
A dataframe with 100,000 rows and 7 columns:
\describe{
\item{X}{true exposure, 1 = present and 0 = absent}
\item{Y}{outcome, 1 = present and 0 = absent}
\item{C1}{1st confounder, 1 = present and 0 = absent}
\item{C2}{2nd confounder, 1 = present and 0 = absent}
\item{C3}{3rd confounder, 1 = present and 0 = absent}
\item{Xstar}{misclassified exposure, 1 = present and 0 = absent}
\item{S}{selection, 1 = selected into the study and 0 = not selected into the study}
}
}
\usage{
df_em_sel_source
}
\description{
Data with complete information on the two sources of bias, three known
confounders, and 100,000 observations. This data is used to derive
\code{df_em_sel} and can be used to obtain bias parameters for purposes
of validating the simultaneous multi-bias adjustment method with
\code{df_em_sel}. With this source data, the fitted regression
\ifelse{html}{\out{logit(P(Y=1)) = &alpha;<sub>0</sub> + &alpha;<sub>1</sub>X + &alpha;<sub>2</sub>C1 + &alpha;<sub>3</sub>C2 + &alpha;<sub>4</sub>C3}}{\eqn{logit(P(Y=1)) = \alpha_0 + \alpha_1 X + \alpha_2 C1 + \alpha_3 C2 + \alpha_4 C3}}
shows that the true, unbiased exposure-outcome odds ratio = 2.
}
\keyword{datasets}
