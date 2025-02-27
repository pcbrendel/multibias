% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adust_om_sel.R
\name{adjust_om_sel}
\alias{adjust_om_sel}
\title{Adust for outcome misclassification and selection bias.}
\usage{
adjust_om_sel(
  data_observed,
  data_validation = NULL,
  y_model_coefs = NULL,
  s_model_coefs = NULL,
  level = 0.95
)
}
\arguments{
\item{data_observed}{Object of class \code{data_observed} corresponding to the
data to perform bias analysis on.}

\item{data_validation}{Object of class \code{data_validation} corresponding to
the validation data used to adjust for bias in the observed data. Here, the
validation data should have data for the same variables as in the observed
data, plus data for the true and misclassified outcome,
corresponding to the observed outcome in \code{data_observed}. There should also
be a selection indicator representing whether the observation in
\code{data_validation} was selected in \code{data_observed}.}

\item{y_model_coefs}{The regression coefficients corresponding to the model:
\ifelse{html}{\out{logit(P(Y=1)) = &delta;<sub>0</sub> + &delta;<sub>1</sub>X + &delta;<sub>2</sub>Y* + &delta;<sub>2+j</sub>C<sub>j</sub>, }}{\eqn{logit(P(Y=1)) = \delta_0 + \delta_1 X + \delta_2 Y^* + \delta_{2+j} C_j, }}
where \emph{Y} represents the binary true outcome, \emph{X} is the exposure,
\emph{Y}* is the binary misclassified outcome, \emph{C} represents
the vector of measured confounders (if any), and \emph{j} corresponds
to the number of measured confounders. The number of parameters is
therefore 3 + \emph{j}.}

\item{s_model_coefs}{The regression coefficients corresponding to the model:
\ifelse{html}{\out{logit(P(S=1)) = &beta;<sub>0</sub> + &beta;<sub>1</sub>X + &beta;<sub>2</sub>Y* + &beta;<sub>2+j</sub>C<sub>j</sub>, }}{\eqn{logit(P(S=1)) = \beta_0 + \beta_1 X + \beta_2 Y^* + \beta_{2+j} C_j, }}
where \emph{S} represents binary selection,
\emph{X} is the exposure, \emph{Y}* is the binary misclassified outcome,
\emph{C} represents the vector of measured confounders (if any),
and \emph{j} corresponds to the number of measured confounders.
The number of parameters is therefore 3 + \emph{j}.}

\item{level}{Value from 0-1 representing the full range of the confidence
interval. Default is 0.95.}
}
\value{
A list where the first item is the odds ratio estimate of the
effect of the exposure on the outcome and the second item is the
confidence interval as the vector: (lower bound, upper bound).
}
\description{
\code{adjust_om_sel} returns the exposure-outcome odds ratio and confidence
interval, adjusted for outcome misclassification and selection bias.
}
\details{
Bias adjustment can be performed by inputting either a validation dataset or
the necessary bias parameters. Values for the bias parameters
can be applied as fixed values or as single draws from a probability
distribution (ex: \code{rnorm(1, mean = 2, sd = 1)}). The latter has
the advantage of allowing the researcher to capture the uncertainty
in the bias parameter estimates. To incorporate this uncertainty in the
estimate and confidence interval, this function should be run in loop across
bootstrap samples of the dataframe for analysis. The estimate and
confidence interval would then be obtained from the median and quantiles
of the distribution of odds ratio estimates.
}
\examples{
df_observed <- data_observed(
  data = df_om_sel,
  exposure = "X",
  outcome = "Ystar",
  confounders = "C1"
)

# Using validation data -----------------------------------------------------
df_validation <- data_validation(
  data = df_om_sel_source,
  true_exposure = "X",
  true_outcome = "Y",
  confounders = "C1",
  misclassified_outcome = "Ystar",
  selection = "S"
)

adjust_om_sel(
  data_observed = df_observed,
  data_validation = df_validation
)

# Using y_model_coefs and s_model_coefs -------------------------------------
adjust_om_sel(
  data_observed = df_observed,
  y_model_coefs = c(-3.24, 0.58, 1.59, 0.45),
  s_model_coefs = c(0.03, 0.92, 0.12, 0.05)
)

}
