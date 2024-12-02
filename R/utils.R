is_binary <- function(x) {
  if (is.logical(x)) {
    return(TRUE)
  }
  if (!is.numeric(x)) {
    return(FALSE)
  }
  unique_vals <- unique(x)
  length(unique_vals) == 2 &&
    all(near(unique_vals, 0) | near(unique_vals, 1))
}

is_continuous <- function(x) {
  is.numeric(x) && length(unique(x)) > 2
}

force_binary <- function(col, message) {
  if (!all(is_binary(col))) {
    stop(message)
  }
}

force_len <- function(len_observed, len_required, message) {
  if (len_observed != len_required) {
    stop(message)
  }
}

force_match <- function(col1, col2, message) {
  col1_type <- if_else(
    is_binary(col1),
    "binary",
    if_else(
      is_continuous(col1),
      "continuous",
      "other"
    )
  )
  col2_type <- if_else(
    is_binary(col2),
    "binary",
    if_else(
      is_continuous(col2),
      "continuous",
      "other"
    )
  )
  if (col1_type != col2_type) {
    stop(message)
  }
}