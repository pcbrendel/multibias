% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/causal_data.R
\name{data_observed}
\alias{data_observed}
\title{Represent observed causal data}
\usage{
data_observed(data, bias, exposure, outcome, confounders = NULL)
}
\arguments{
\item{data}{Dataframe for bias analysis.}

\item{bias}{String type(s) of bias distorting the effect of the exposure
on the outcome. Can choose from a subset of the following: "uc", "em", "om",
"sel". These correspond to uncontrolled confounding, exposure
misclassification, outcome misclassification, and selection bias,
respectively.}

\item{exposure}{String name of the column in \code{data} corresponding to the
exposure variable.}

\item{outcome}{String name of the column in \code{data} corresponding to the
outcome variable.}

\item{confounders}{String name(s) of the column(s) in \code{data} corresponding
to the confounding variable(s).}
}
\value{
An object of class \code{data_observed} containing:
\item{data}{A dataframe with the selected columns}
\item{bias}{The type(s) of bias present}
\item{exposure}{The name of the exposure variable}
\item{outcome}{The name of the outcome variable}
\item{confounders}{The name(s) of the confounder variable(s)}
}
\description{
\code{data_observed} combines the observed dataframe with specific identification
of the columns corresponding to the exposure, outcome, and confounders. It is
an essential input of the \code{\link[=multibias_adjust]{multibias_adjust()}} function.
}
\examples{
df <- data_observed(
  data = df_sel,
  bias = "uc",
  exposure = "X",
  outcome = "Y",
  confounders = c("C1", "C2", "C3")
)

}
