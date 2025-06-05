#' Represent observed causal data
#'
#' @description
#' `data_observed` combines the observed dataframe with specific identification
#' of the columns corresponding to the exposure, outcome, and confounders. It is
#' an essential input of the [multibias_adjust()] function.
#'
#' @param data Dataframe for bias analysis.
#' @param bias String type(s) of bias distorting the effect of the exposure
#' on the outcome. Can choose from a subset of the following: "uc", "em", "om",
#' "sel". These correspond to uncontrolled confounding, exposure
#' misclassification, outcome misclassification, and selection bias,
#' respectively.
#' @param exposure String name of the column in `data` corresponding to the
#' exposure variable.
#' @param outcome String name of the column in `data` corresponding to the
#' outcome variable.
#' @param confounders String name(s) of the column(s) in `data` corresponding
#' to the confounding variable(s).
#'
#' @return An object of class `data_observed` containing:
#'   \item{data}{A dataframe with the selected columns}
#'   \item{bias}{The type(s) of bias present}
#'   \item{exposure}{The name of the exposure variable}
#'   \item{outcome}{The name of the outcome variable}
#'   \item{confounders}{The name(s) of the confounder variable(s)}
#'
#' @examples
#' df <- data_observed(
#'   data = df_sel,
#'   bias = "uc",
#'   exposure = "X",
#'   outcome = "Y",
#'   confounders = c("C1", "C2", "C3")
#' )
#'
#' @export

data_observed <- function(
    data,
    bias,
    exposure,
    outcome,
    confounders = NULL) {
  stopifnot(
    is.data.frame(data),
    is.character(bias),
    is.character(exposure) & length(exposure) == 1,
    is.character(outcome) & length(outcome) == 1,
    is.character(confounders) | is.null(confounders)
  )

  required_vars <- c(exposure, outcome, confounders)
  required_vars <- required_vars[!is.null(required_vars)]

  if (!all(required_vars %in% names(data))) {
    missing_vars <- required_vars[!required_vars %in% names(data)]
    stop(
      paste0(
        "Input variable(s) missing in data: ",
        paste(missing_vars, collapse = ", ")
      )
    )
  }

  df <- data %>%
    select(all_of(c(exposure, outcome, confounders)))

  acceptable_biases <- c("uc", "em", "om", "sel")
  if (!all(bias %in% acceptable_biases)) {
    stop(
      paste0(
        "Unacceptable bias input. ",
        "Biases must include a subset of the following: ",
        "uc, em, om, sel"
      )
    )
  }

  obj <- structure(
    list(
      data = df,
      bias = bias,
      exposure = exposure,
      outcome = outcome,
      confounders = confounders
    ),
    class = "data_observed"
  )

  return(obj)
}


#' Print method for data_observed objects
#'
#' @description
#' Prints a formatted summary of a `data_observed` object, including:
#' - The types of biases present
#' - The exposure, outcome, and confounder variables
#' - A preview of the first 5 rows of data
#'
#' @param x A `data_observed` object
#' @param ... Additional arguments passed to print
#'
#' @return The input object invisibly, allowing for method chaining
#'
#' @method print data_observed
#' @keywords internal
#' @export

print.data_observed <- function(x, ...) {
  bias_map <- c(
    "uc" = "Uncontrolled Confounding",
    "em" = "Exposure Misclassification",
    "om" = "Outcome Misclassification",
    "sel" = "Selection Bias"
  )
  bias_description <- bias_map[x$bias]

  cat("Observed Data\n")
  cat("---------------------------------\n")
  cat("The following biases are present:", "\n")
  cat(paste(bias_description, collapse = "\n"), "\n")
  cat("---------------------------------\n")
  cat("Exposure:", x$exposure, "\n")
  cat("Outcome:", x$outcome, "\n")
  if (!is.null(x$confounders)) {
    cat("Confounders:", paste(x$confounders, collapse = ", "), "\n")
  }
  cat("---------------------------------\n")
  cat("Data head: \n")
  print(x$data[seq_len(min(5, nrow(x$data))), ])
  invisible(x)
}


#' Summary method for data_observed objects
#'
#' @description
#' Provides a statistical summary of the observed data by fitting either:
#' - A logistic regression model for binary outcomes
#' - A linear regression model for continuous outcomes
#'
#' The model includes the exposure and all confounders as predictors.
#' For binary outcomes, estimates are exponentiated to show odds ratios.
#'
#' @param object A `data_observed` object
#' @param ... Additional arguments passed to summary
#'
#' @return A data frame containing model coefficients, standard errors,
#'         confidence intervals, and p-values. For binary outcomes,
#'         coefficients are exponentiated to show odds ratios.
#'
#' @method summary data_observed
#' @importFrom broom tidy
#' @keywords internal
#' @export

summary.data_observed <- function(object, ...) {
  df <- data.frame(
    exposure = object$data[, object$exposure],
    Y = object$data[, object$outcome]
  )
  names(df)[1] <- object$exposure

  df <- bind_cols(
    df,
    object$data %>%
      select(all_of(object$confounders))
  )

  if (all(df$Y %in% 0:1)) {
    y_binary <- TRUE
  } else {
    y_binary <- FALSE
  }

  if (y_binary) {
    mod <- glm(
      Y ~ .,
      family = binomial(link = "logit"),
      data = df
    )
    final <- broom::tidy(mod, conf.int = TRUE, exponentiate = TRUE)
    cat(
      "Note: Estimates are exponentiated (odds ratios) for binary outcomes\n\n"
    )
  } else {
    mod <- lm(
      Y ~ .,
      data = df
    )
    final <- broom::tidy(mod, conf.int = TRUE, exponentiate = FALSE)
    cat("Note: Estimates are not exponentiated for continuous outcomes\n\n")
  }

  return(final)
}


#' Represent validation causal data
#'
#' @description
#' `data_validation` is one of two different options to represent bias
#' assumptions for bias adjustment. It combines the validation dataframe
#' with specific identification of the appropriate columns for bias adjustment,
#' including: true exposure, true outcome, confounders, misclassified exposure,
#' misclassified outcome, and selection. The purpose of validation data is to
#' use an external data source to transport the necessary causal relationships
#' that are missing in the observed data.
#'
#' @param data Dataframe of validation data
#' @param true_exposure String name of the column in `data` corresponding to
#' the true exposure.
#' @param true_outcome String name of the column in `data` corresponding to
#' the true outcome.
#' @param confounders String name(s) of the column(s) in `data` corresponding
#' to the confounding variable(s).
#' @param misclassified_exposure String name of the column in `data`
#' corresponding to the misclassified exposure.
#' @param misclassified_outcome String name of the column in `data`
#' corresponding to the misclassified outcome.
#' @param selection String name of the column in `data` corresponding to the
#' selection indicator.
#'
#' @return An object of class `data_validation` containing:
#'   \item{data}{A dataframe with the selected columns}
#'   \item{true_exposure}{The name of the true exposure variable}
#'   \item{true_outcome}{The name of the true outcome variable}
#'   \item{confounders}{The name(s) of the confounder variable(s)}
#'   \item{misclassified_exposure}{The name of the misclassified exposure variable}
#'   \item{misclassified_outcome}{The name of the misclassified outcome variable}
#'   \item{selection}{The name of the selection indicator variable}
#'
#' @examples
#' df <- data_validation(
#'   data = df_sel_source,
#'   true_exposure = "X",
#'   true_outcome = "Y",
#'   confounders = c("C1", "C2", "C3"),
#'   selection = "S"
#' )
#'
#' @export

data_validation <- function(
    data,
    true_exposure,
    true_outcome,
    confounders = NULL,
    misclassified_exposure = NULL,
    misclassified_outcome = NULL,
    selection = NULL) {
  stopifnot(
    is.data.frame(data),
    is.character(true_exposure) & length(true_exposure) == 1,
    is.character(true_outcome) & length(true_outcome) == 1,
    is.character(confounders) | is.null(confounders),
    (is.character(misclassified_exposure) &
      length(misclassified_exposure) == 1) |
      is.null(misclassified_exposure),
    (is.character(misclassified_outcome) &
      length(misclassified_outcome) == 1) |
      is.null(misclassified_outcome),
    (is.character(selection) & length(selection) == 1) | is.null(selection)
  )

  required_vars <- c(
    true_exposure, true_outcome,
    confounders, misclassified_exposure,
    misclassified_outcome, selection
  )
  required_vars <- required_vars[!is.null(required_vars)]

  if (!all(required_vars %in% names(data))) {
    missing_vars <- required_vars[!required_vars %in% names(data)]
    stop(
      paste0(
        "Input variable(s) missing in data: ",
        paste(missing_vars, collapse = ", ")
      )
    )
  }

  cols <- c(
    true_exposure,
    true_outcome,
    confounders,
    misclassified_exposure,
    misclassified_outcome,
    selection
  )

  df <- data %>%
    select(all_of(cols))

  obj <- list(
    data = df,
    true_exposure = true_exposure,
    true_outcome = true_outcome,
    confounders = confounders,
    misclassified_exposure = misclassified_exposure,
    misclassified_outcome = misclassified_outcome,
    selection = selection
  )

  class(obj) <- "data_validation"
  return(obj)
}


#' Print method for data_validation objects
#'
#' @description
#' Prints a formatted summary of a `data_validation` object, including:
#' - The true exposure and outcome variables
#' - Any confounders, misclassified variables, or selection indicators
#' - A preview of the first 5 rows of data
#'
#' @param x A `data_validation` object
#' @param ... Additional arguments passed to print
#'
#' @return The input object invisibly
#'
#' @keywords internal
#' @export

print.data_validation <- function(x, ...) {
  cat("Validation Data\n")
  cat("------------------\n")
  cat("True exposure:", x$true_exposure, "\n")
  cat("True outcome:", x$true_outcome, "\n")
  if (!is.null(x$confounders)) {
    cat("Confounders:", paste(x$confounders, collapse = ", "), "\n")
  }
  if (!is.null(x$misclassified_exposure)) {
    cat("Misclassified exposure:", x$misclassified_exposure, "\n")
  }
  if (!is.null(x$misclassified_outcome)) {
    cat("Misclassified outcome:", x$misclassified_outcome, "\n")
  }
  if (!is.null(x$selection)) {
    cat("Selection:", x$selection, "\n")
  }
  cat("Data head: \n")
  print(x$data[seq_len(min(5, nrow(x$data))), ])
  invisible(x)
}
