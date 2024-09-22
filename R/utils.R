force_binary <- function(col, message) {
  if (!all(col %in% 0:1)) {
    stop(message)
  }
}

force_len <- function(len_observed, len_required, message) {
  if (len_observed != len_required) {
    stop(message)
  }
}