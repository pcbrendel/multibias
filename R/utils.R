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
    stop(message, call. = FALSE)
  }
}

force_len <- function(len_observed, len_required, message) {
  if (len_observed != len_required) {
    stop(message, call. = FALSE)
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
    stop(message, call. = FALSE)
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
  input1_name <- deparse(substitute(input1))
  input2_name <- deparse(substitute(input2))

  if (!is_any_null(input1)) {
    if (!is_all_null(input2)) {
      stop(
        paste0(
          "No other bias-adjusting parameters should be specified when ",
          input1_name,
          " is used."
        ),
        call. = FALSE
      )
    }
  } else if (!is_any_null(input2)) {
    if (!is_all_null(input1)) {
      stop(
        paste0(
          "No other bias-adjusting parameters should be specified when ",
          input2_name,
          " is used."
        ),
        call. = FALSE
      )
    }
  } else if (all(is_any_null(input1), is_any_null(input2))) {
    stop(
      paste0(
        "One of:\n",
        "1. ", input1_name, "\n",
        "2. ", input2_name, "\n",
        "must be non-null."
      ),
      call. = FALSE
    )
  }
}

check_inputs3 <- function(input1, input2, input3, ignore = NULL) {
  input1_name <- deparse(substitute(input1))
  input2_name <- deparse(substitute(input2))
  input3_name <- deparse(substitute(input3))

  remove_ignore <- function(x, ignore) {
    if (is.list(x) && !inherits(x, "data_validation")) {
      x[sapply(x, function(y) !(identical(y, ignore)))]
    } else {
      x
    }
  }

  if (!is.null(ignore)) {
    input1 <- remove_ignore(input1, ignore)
    input2 <- remove_ignore(input2, ignore)
    input3 <- remove_ignore(input3, ignore)
  }

  if (!is_any_null(input1)) {
    if (!all(is_all_null(input2), is_all_null(input3))) {
      stop(
        paste0(
          "No other bias-adjusting parameters should be specified when ",
          input1_name,
          " is used."
        ),
        call. = FALSE
      )
    }
  } else if (!is_any_null(input2)) {
    if (!all(is_all_null(input1), is_all_null(input3))) {
      stop(
        paste0(
          "No other bias-adjusting parameters should be specified when ",
          input2_name,
          " is used."
        ),
        call. = FALSE
      )
    }
  } else if (!is_any_null(input3)) {
    if (!all(is_all_null(input1), is_all_null(input2))) {
      stop(
        paste0(
          "No other bias-adjusting parameters should be specified when ",
          input3_name,
          " is used."
        ),
        call. = FALSE
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
        "1. ", input1_name, "\n",
        "2. ", input2_name, "\n",
        "3. ", input3_name, "\n",
        "must be non-null."
      ),
      call. = FALSE
    )
  }
}