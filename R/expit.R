# Inverse logit (expit) function

expit <- function(p) {
  exp(p) / (1 + exp(p))
}
