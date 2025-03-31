set.seed(1234)
n <- 50000
nreps <- 10

# SEPARATE BIAS PARAMETERS

# 0 confounders

nobias_model <- glm(
  Y ~ X,
  family = binomial(link = "logit"),
  data = df_em_om_source
)

x_model <- glm(
  X ~ Xstar + Ystar,
  family = binomial(link = "logit"),
  data = df_em_om_source
)

y_model <- glm(
  Y ~ X + Ystar,
  family = binomial(link = "logit"),
  data = df_em_om_source
)

df_observed <- data_observed(
  df_em_om,
  bias = c("em", "om"),
  exposure = "Xstar",
  outcome = "Ystar",
  confounders = NULL
)
list_for_em_om <- list(
  x = as.vector(coef(x_model)),
  y = as.vector(coef(y_model))
)
bp_em_om <- bias_params(coef_list = list_for_em_om)

single_run <- adjust_em_om(
  df_observed,
  bias_params = bp_em_om
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_em_om[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    bias = c("em", "om"),
    exposure = "Xstar",
    outcome = "Ystar",
    confounders = NULL
  )
  results <- adjust_em_om(
    df_observed,
    bias_params = bp_em_om
  )
  est[i] <- results$estimate
}

or_true <- exp(summary(nobias_model)$coef[2, 1])
or_adjusted <- median(est)

test_that("odds ratio and confidence interval output", {
  expect_gt(or_adjusted, or_true - 0.1)
  expect_lt(or_adjusted, or_true + 0.1)
  expect_vector(
    single_run$ci,
    ptype = double(),
    size = 2
  )
})

# 3 confounders

nobias_model <- glm(
  Y ~ X + C1 + C2 + C3,
  family = binomial(link = "logit"),
  data = df_em_om_source
)

x_model <- glm(
  X ~ Xstar + Ystar + C1 + C2 + C3,
  family = binomial(link = "logit"),
  data = df_em_om_source
)

y_model <- glm(
  Y ~ X + Ystar + C1 + C2 + C3,
  family = binomial(link = "logit"),
  data = df_em_om_source
)

df_observed <- data_observed(
  df_em_om,
  bias = c("em", "om"),
  exposure = "Xstar",
  outcome = "Ystar",
  confounders = c("C1", "C2", "C3")
)
list_for_em_om <- list(
  x = as.vector(coef(x_model)),
  y = as.vector(coef(y_model))
)
bp_em_om <- bias_params(coef_list = list_for_em_om)

single_run <- adjust_em_om(
  df_observed,
  bias_params = bp_em_om
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_em_om[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    bias = c("em", "om"),
    exposure = "Xstar",
    outcome = "Ystar",
    confounders = c("C1", "C2", "C3")
  )
  results <- adjust_em_om(
    df_observed,
    bias_params = bp_em_om
  )
  est[i] <- results$estimate
}

or_true <- exp(summary(nobias_model)$coef[2, 1])
or_adjusted <- median(est)

test_that("odds ratio and confidence interval output", {
  expect_gt(or_adjusted, or_true - 0.1)
  expect_lt(or_adjusted, or_true + 0.1)
  expect_vector(
    single_run$ci,
    ptype = double(),
    size = 2
  )
})

# MULTINOMIAL BIAS PARAMETERS

# library(nnet)

set.seed(1234)
n <- 10000
nreps <- 10

# 0 confounders

nobias_model <- glm(
  Y ~ X,
  family = binomial(link = "logit"),
  data = df_em_om_source
)

# xy_model <- multinom(
#   paste(X, Y) ~ Xstar + Ystar,
#   data = df_em_om_source
# )
# summary(xy_model)

df_observed <- data_observed(
  df_em_om,
  bias = c("em", "om"),
  exposure = "Xstar",
  outcome = "Ystar",
  confounders = NULL
)
list_for_em_om <- list(
  x1y0 = c(-1.99, 1.62, 0.23),
  x0y1 = c(-2.96, 0.22, 1.60),
  x1y1 = c(-4.36, 1.82, 1.82)

)
bp_em_om <- bias_params(coef_list = list_for_em_om)

single_run <- adjust_em_om(
  df_observed,
  bias_params = bp_em_om
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_em_om[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    bias = c("em", "om"),
    exposure = "Xstar",
    outcome = "Ystar",
    confounders = NULL
  )
  results <- adjust_em_om(
    df_observed,
    bias_params = bp_em_om
  )
  est[i] <- results$estimate
}

or_true <- exp(summary(nobias_model)$coef[2, 1])
or_adjusted <- median(est)

test_that("0 confounders: odds ratio and confidence interval output", {
  expect_gt(or_adjusted, or_true - 0.1)
  expect_lt(or_adjusted, or_true + 0.1)
  expect_vector(
    single_run$ci,
    ptype = double(),
    size = 2
  )
})

# 3 confounders

nobias_model <- glm(
  Y ~ X + C1 + C2 + C3,
  family = binomial(link = "logit"),
  data = df_em_om_source
)

# xy_model <- multinom(
#   paste(X, Y) ~ Xstar + Ystar + C1 + C2 + C3,
#   data = df_em_om_source
# )
# summary(xy_model)

df_observed <- data_observed(
  df_em_om,
  bias = c("em", "om"),
  exposure = "Xstar",
  outcome = "Ystar",
  confounders = c("C1", "C2", "C3")
)
list_for_em_om <- list(
  x1y0 = c(-2.86, 1.63, 0.23, 0.37, -0.22, 0.87),
  x0y1 = c(-3.26, 0.22, 1.60, 0.41, -0.93, 0.28),
  x1y1 = c(-5.62, 1.83, 1.83, 0.74, -1.15, 1.19)

)
bp_em_om <- bias_params(coef_list = list_for_em_om)

single_run <- adjust_em_om(
  df_observed,
  bias_params = bp_em_om
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_em_om[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    bias = c("em", "om"),
    exposure = "Xstar",
    outcome = "Ystar",
    confounders = c("C1", "C2", "C3")
  )
  results <- adjust_em_om(
    df_observed,
    bias_params = bp_em_om
  )
  est[i] <- results$estimate
}

or_true <- exp(summary(nobias_model)$coef[2, 1])
or_adjusted <- median(est)

test_that("3 confounders: odds ratio and confidence interval output", {
  expect_gt(or_adjusted, or_true - 0.1)
  expect_lt(or_adjusted, or_true + 0.1)
  expect_vector(
    single_run$ci,
    ptype = double(),
    size = 2
  )
})

# adjust with validation data

or_val <- adjust_em_om(
  data_observed = data_observed(
    df_em_om,
    bias = c("em", "om"),
    exposure = "Xstar",
    outcome = "Ystar",
    confounders = c("C1", "C2", "C3")
  ),
  data_validation = data_validation(
    df_em_om_source,
    true_exposure = "X",
    true_outcome = "Y",
    confounders = c("C1", "C2", "C3"),
    misclassified_exposure = "Xstar",
    misclassified_outcome = "Ystar"
  )
)

test_that("adjust_em_om, validation data", {
  expect_gt(or_val$estimate, or_true - 0.1)
  expect_lt(or_val$estimate, or_true + 0.1)
})
