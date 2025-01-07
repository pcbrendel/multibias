#' Represent observed causal data
#'
#' @description
#' `data_observed` combines the observed dataframe with specific identification
#' of the columns corresponding to the exposure, outcome, and confounders. It is
#' an essential input of all `adjust` functions.
#'
#' @param data Dataframe for bias analysis.
#' @param exposure String name of the column in `data` corresponding to the
#' exposure variable.
#' @param outcome String name of the column in `data` corresponding to the
#' outcome variable.
#' @param confounders String name(s) of the column(s) in `data` corresponding
#' to the confounding variable(s).
#'
#' @examples
#' df <- data_observed(
#'   data = df_sel,
#'   exposure = "X",
#'   outcome = "Y",
#'   confounders = c("C1", "C2", "C3")
#' )
#'
#' @export

data_observed <- function(
    data,
    exposure,
    outcome,
    confounders = NULL) {
  stopifnot(
    is.data.frame(data),
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

  obj <- structure(
    list(
      data = df,
      exposure = exposure,
      outcome = outcome,
      confounders = confounders
    ),
    class = "data_observed"
  )

  return(obj)
}

#' @export

print.data_observed <- function(x, ...) {
  cat("Observed Data\n")
  cat("------------------\n")
  cat("Exposure:", x$exposure, "\n")
  cat("Outcome:", x$outcome, "\n")
  if (!is.null(x$confounders)) {
    cat("Confounders:", paste(x$confounders, collapse = ", "), "\n")
  }
  cat("Data head: \n")
  print(x$data[seq_len(min(5, nrow(x$data))), ])
  invisible(x)
}


#' Represent validation causal data
#'
#' @description
#' `data_validation` combines the validation dataframe with specific
#' identification of the appropriate columns for bias adjustment, including:
#' true exposure, true outcome, confounders, misclassified exposure,
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
