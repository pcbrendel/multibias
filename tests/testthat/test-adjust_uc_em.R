set.seed(1234)
n <- 50000
nreps <- 10

# SEPARATE BIAS PARAMETERS

# cont Y just for testing that function runs
df_uc_em$Y_cont <- plogis(df_uc_em$Y) +
  rnorm(nrow(df_uc_em), mean = 0, sd = 0.1)

# 0 confounders

nobias_model <- glm(
  Y ~ X + U,
  family = binomial(link = "logit"),
  data = df_uc_em_source
)

u_model <- glm(
  U ~ X + Y,
  family = binomial(link = "logit"),
  data = df_uc_em_source
)
x_model <- glm(
  X ~ Xstar + Y,
  family = binomial(link = "logit"),
  data = df_uc_em_source
)

df_observed <- data_observed(
  df_uc_em,
  bias = c("uc", "em"),
  exposure = "Xstar",
  outcome = "Y",
  confounders = NULL
)
list_for_uc_em <- list(
  u = as.vector(coef(u_model)),
  x = as.vector(coef(x_model))
)
bp_uc_em <- bias_params(coef_list = list_for_uc_em)

single_run <- adjust_uc_em(
  df_observed,
  bias_params = bp_uc_em
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_em[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    bias = c("uc", "em"),
    exposure = "Xstar",
    outcome = "Y",
    confounders = NULL
  )
  results <- adjust_uc_em(
    df_observed,
    bias_params = bp_uc_em
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
  Y ~ X + C1 + C2 + C3 + U,
  family = binomial(link = "logit"),
  data = df_uc_em_source
)

u_model <- glm(
  U ~ X + Y,
  family = binomial(link = "logit"),
  data = df_uc_em_source
)
x_model <- glm(
  X ~ Xstar + Y + C1 + C2 + C3,
  family = binomial(link = "logit"),
  data = df_uc_em_source
)

df_observed <- data_observed(
  df_uc_em,
  bias = c("uc", "em"),
  exposure = "Xstar",
  outcome = "Y",
  confounders = c("C1", "C2", "C3")
)
list_for_uc_em <- list(
  u = as.vector(coef(u_model)),
  x = as.vector(coef(x_model))
)
bp_uc_em <- bias_params(coef_list = list_for_uc_em)

single_run <- adjust_uc_em(
  df_observed,
  bias_params = bp_uc_em
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_em[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    bias = c("uc", "em"),
    exposure = "Xstar",
    outcome = "Y",
    confounders = c("C1", "C2", "C3")
  )
  results <- adjust_uc_em(
    df_observed,
    bias_params = bp_uc_em
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
n <- 20000
nreps <- 10

# 0 confounders

nobias_model <- glm(
  Y ~ X + U,
  family = binomial(link = "logit"),
  data = df_uc_em_source
)

# xu_model <- multinom(
#   paste0(X, U) ~ Xstar + Y,
#   data = df_uc_em_source
# )
# summary(xu_model)

df_observed <- data_observed(
  df_uc_em,
  bias = c("uc", "em"),
  exposure = "Xstar",
  outcome = "Y",
  confounders = NULL
)
list_for_uc_em <- list(
  x1u0 = c(-1.91, 1.63, 0.71),
  x0u1 = c(-0.23, 0.02, 0.01),
  x1u1 = c(-1.47, 1.63, 0.71)

)
bp_uc_em <- bias_params(coef_list = list_for_uc_em)

single_run <- adjust_uc_em(
  df_observed,
  bias_params = bp_uc_em
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_em[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    bias = c("uc", "em"),
    exposure = "Xstar",
    outcome = "Y",
    confounders = NULL
  )
  results <- adjust_uc_em(
    df_observed,
    bias_params = bp_uc_em
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
  Y ~ X + C1 + C2 + C3 + U,
  family = binomial(link = "logit"),
  data = df_uc_em_source
)

# xu_model <- multinom(
#   paste0(X, U) ~ Xstar + Y + C1 + C2 + C3,
#   data = df_uc_em_source
# )
# summary(xu_model)

df_observed <- data_observed(
  df_uc_em,
  bias = c("uc", "em"),
  exposure = "Xstar",
  outcome = "Y",
  confounders = c("C1", "C2", "C3")
)
list_for_uc_em <- list(
  x1u0 = c(-2.79, 1.63, 0.63, 0.35, -0.21, 0.89),
  x0u1 = c(-0.12, 0.02, 0.02, -0.04, 0.02, -0.12),
  x1u1 = c(-2.20, 1.63, 0.64, 0.29, -0.15, 0.74)

)
bp_uc_em <- bias_params(coef_list = list_for_uc_em)

single_run <- adjust_uc_em(
  df_observed,
  bias_params = bp_uc_em
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_em[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    bias = c("uc", "em"),
    exposure = "Xstar",
    outcome = "Y",
    confounders = c("C1", "C2", "C3")
  )
  results <- adjust_uc_em(
    df_observed,
    bias_params = bp_uc_em
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

# adjust with validation data

or_val <- adjust_uc_em(
  data_observed = data_observed(
    df_uc_em,
    bias = c("uc", "em"),
    exposure = "Xstar",
    outcome = "Y",
    confounders = c("C1", "C2", "C3")
  ),
  data_validation = data_validation(
    df_uc_em_source,
    true_exposure = "X",
    true_outcome = "Y",
    confounders = c("C1", "C2", "C3", "U"),
    misclassified_exposure = "Xstar"
  )
)

test_that("adjust_uc_em, validation data", {
  expect_gt(or_val$estimate, or_true - 0.1)
  expect_lt(or_val$estimate, or_true + 0.1)
})
