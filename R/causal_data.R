#' Represent observed causal data
#'
#' @description
#' `data_observed` combines the observed dataframe with specific identification
#' of the columns corresponding to the exposure, outcome, and confounders. It is
#' an essential input of all `adjust` functions.
#'
#' @param data Dataframe for bias analysis.
#' @param exposure String name of the exposure variable.
#' @param outcome String name of the outcome variable.
#' @param confounders String name(s) of the confounder(s).
#'
#' @examples
#' df <- data_observed(
#'   data = df_uc,
#'   exposure = "X",
#'   outcome = "Ystar",
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

  obj <- list(
    data = data,
    exposure = exposure,
    outcome = outcome,
    confounders = confounders
  )

  class(obj) <- "data_observed"
  return(obj)
}

print.data_observed <- function(x, ...) {
  cat("Observed Data\n")
  cat("------------------\n")
  cat("Exposure:", x$exposure, "\n")
  cat("Outcome:", x$outcome, "\n")
  if (!is.null(x$confounders)) {
    cat("Confounders:", paste(x$confounders, collapse = ", "), "\n")
  }
  cat("Data head: \n")
  print(x$data[1:5, ])
}


# Class for Validation Data

# data_validation <- function(
#     data,
#     true_exposure,
#     true_outcome,
#     confounders = NULL,
#     misclassified_exposure = NULL,
#     misclassified_outcome = NULL,
#     selection = NULL,
#     remove = NULL) {
#   stopifnot(
#     is.data.frame(data),
#     is.character(true_exposure) & length(true_exposure) == 1,
#     is.character(true_outcome) & length(true_outcome) == 1,
#     is.character(confounders) | is.null(confounders),
#     (is.character(misclassified_exposure) &
#       length(misclassified_exposure) == 1) |
#       is.null(misclassified_exposure),
#     (is.character(misclassified_outcome) &
#       length(misclassified_outcome) == 1) |
#       is.null(misclassified_outcome),
#     (is.character(selection) & length(selection) == 1) | is.null(selection)
#   )

#   obj <- list(
#     data = data[, !colnames(data) %in% remove],
#     true_exposure = true_exposure,
#     true_outcome = true_outcome,
#     confounders = confounders,
#     misclassified_exposure = misclassified_exposure,
#     misclassified_outcome = misclassified_outcome,
#     selection = selection
#   )

#   class(obj) <- "data_validation"
#   return(obj)
# }

# print.data_validation <- function(x, ...) {
#   cat("Validation Data\n")
#   cat("------------------\n")
#   cat("True exposure:", x$true_exposure, "\n")
#   cat("True outcome:", x$true_outcome, "\n")
#   if (!is.null(x$confounders)) {
#     cat("Confounders:", paste(x$confounders, collapse = ", "), "\n")
#   }
#   if (!is.null(x$misclassified_exposure)) {
#     cat("Misclassified exposure:", x$misclassified_exposure, "\n")
#   }
#   if (!is.null(x$misclassified_outcome)) {
#     cat("Misclassified outcome:", x$misclassified_outcome, "\n")
#   }
#   if (!is.null(x$selection)) {
#     cat("Selection:", x$selection, "\n")
#   }
#   cat("Data head: \n")
#   print(x$data[1:5, ])
# }
