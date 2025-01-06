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

is_any_null <- function(x) {
  if (is.list(x) && inherits(x, "data_validation") == FALSE) {
    any(sapply(x, is.null))
  } else {
    is.null(x)
  }
}

is_all_null <- function(x) {
  if (is.list(x) && inherits(x, "data_validation") == FALSE) {
    all(sapply(x, is.null))
  } else {
    is.null(x)
  }
}

check_inputs2 <- function(input1, input2) {
  if (!is_any_null(input1)) {
    if (!is_all_null(input2)) {
      stop(
        paste0(
          "No other bias-adjusting parameters should be specified when ",
          deparse(substitute(input1)),
          " is used."
        )
      )
    }
  } else if (!is_any_null(input2)) {
    if (!is_all_null(input1)) {
      stop(
        paste0(
          "No other bias-adjusting parameters should be specified when ",
          deparse(substitute(input2)),
          " is used."
        )
      )
    }
  } else if (all(is_any_null(input1), is_any_null(input2))) {
    stop(
      paste0(
        "One of:\n",
        "1. ", deparse(substitute(input1)), "\n",
        "2. ", deparse(substitute(input2)), "\n",
        "must be non-null."
      )
    )
  }
}

check_inputs3 <- function(input1, input2, input3) {
  if (!is_any_null(input1)) {
    if (!all(is_all_null(input2), is_all_null(input3))) {
      stop(
        paste0(
          "No other bias-adjusting parameters should be specified when ",
          deparse(substitute(input1)),
          " is used."
        )
      )
    }
  } else if (!is_any_null(input2)) {
    if (!all(is_all_null(input1), is_all_null(input3))) {
      stop(
        paste0(
          "No other bias-adjusting parameters should be specified when ",
          deparse(substitute(input2)),
          " is used."
        )
      )
    }
  } else if (!is_any_null(input3)) {
    if (!all(is_all_null(input1), is_all_null(input2))) {
      stop(
        paste0(
          "No other bias-adjusting parameters should be specified when ",
          deparse(substitute(input3)),
          " is used."
        )
      )
    }
  } else if (
    all(
      is_any_null(input1),
      is_any_null(input2),
      is_any_null(input3)
    )
  ) {
    stop(
      paste0(
        "One of:\n",
        "1. ", deparse(substitute(input1)), "\n",
        "2. ", deparse(substitute(input2)), "\n",
        "3. ", deparse(substitute(input3)), "\n",
        "must be non-null."
      )
    )
  }
}