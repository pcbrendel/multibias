#' Represent bias parameters
#'
#' @description
#' `bias_params` ...
#'
#' @param coef_list List ...
#'
#' @examples
#' list_for_uc <- list(
#'   u = c(-0.19, 0.61, 0.70, -0.09, 0.10, -0.15)
#' )
#'
#' bp <- bias_params(coef_list = list_for_uc)
#'
#' @export

bias_params <- function(coef_list) {
  # Check if the input is a list (R's equivalent of a dictionary)
  if (!is.list(coef_list)) {
    stop("Input must be a list (dictionary).")
  }

  # Check if all keys are strings
  if (!all(sapply(names(coef_list), is.character))) {
    stop("All keys in the dictionary must be strings.")
  }

  # Check if all values are vectors
  if (!all(sapply(coef_list, is.vector))) {
    stop("All values in the dictionary must be vectors.")
  }

  obj <- list(data = coef_list)
  class(obj) <- "bias_params"

  return(obj)
}

#' @export

print.bias_params <- function(x, ...) {
  cat("Bias Parameters\n")
  cat("---------------------------------\n")
  for (key in names(x$data)) {
    cat("Model", key, "\n")
    cat("  Coefs: ", paste(x$data[[key]], collapse = ", "), "\n")
  }
  invisible(x)
}