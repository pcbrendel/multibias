set.seed(1234)
n <- 50000
nreps <- 10

# SEPARATE BIAS PARAMETERS

# cont Y just for testing that function runs
df_uc_em_sel$Y_cont <- plogis(df_uc_em_sel$Y) +
  rnorm(nrow(df_uc_em_sel), mean = 0, sd = 0.1)

# 0 confounders

nobias_model <- glm(
  Y ~ X + U,
  family = binomial(link = "logit"),
  data = df_uc_em_sel_source
)

u_model <- glm(
  U ~ X + Y,
  family = binomial(link = "logit"),
  data = df_uc_em_sel_source
)
x_model <- glm(
  X ~ Xstar + Y,
  family = binomial(link = "logit"),
  data = df_uc_em_sel_source
)
s_model <- glm(
  S ~ Xstar + Y,
  family = binomial(link = "logit"),
  data = df_uc_em_sel_source
)

df_observed <- data_observed(
  df_uc_em_sel,
  bias = c("uc", "em", "sel"),
  exposure = "Xstar",
  outcome = "Y_cont",
  confounders = NULL
)
list_for_uc_em_sel <- list(
  u = as.vector(coef(u_model)),
  x = as.vector(coef(x_model)),
  s = as.vector(coef(s_model))
)
bp_uc_em_sel <- bias_params(coef_list = list_for_uc_em_sel)

single_run <- adjust_uc_em_sel(
  df_observed,
  bias_params = bp_uc_em_sel
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_em_sel[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    bias = c("uc", "em", "sel"),
    exposure = "Xstar",
    outcome = "Y",
    confounders = NULL
  )
  results <- adjust_uc_em_sel(
    df_observed,
    bias_params = bp_uc_em_sel
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
  Y ~ X + C1 + C2 + C3 + U,
  family = binomial(link = "logit"),
  data = df_uc_em_sel_source
)

u_model <- glm(
  U ~ X + Y,
  family = binomial(link = "logit"),
  data = df_uc_em_sel_source
)
x_model <- glm(
  X ~ Xstar + Y + C1 + C2 + C3,
  family = binomial(link = "logit"),
  data = df_uc_em_sel_source
)
s_model <- glm(
  S ~ Xstar + Y + C1 + C2 + C3,
  family = binomial(link = "logit"),
  data = df_uc_em_sel_source
)

df_observed <- data_observed(
  df_uc_em_sel,
  bias = c("uc", "em", "sel"),
  exposure = "Xstar",
  outcome = "Y_cont",
  confounders = c("C1", "C2", "C3")
)
list_for_uc_em_sel <- list(
  u = as.vector(coef(u_model)),
  x = as.vector(coef(x_model)),
  s = as.vector(coef(s_model))
)
bp_uc_em_sel <- bias_params(coef_list = list_for_uc_em_sel)

single_run <- adjust_uc_em_sel(
  df_observed,
  bias_params = bp_uc_em_sel
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_em_sel[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    bias = c("uc", "em", "sel"),
    exposure = "Xstar",
    outcome = "Y",
    confounders = c("C1", "C2", "C3")
  )
  results <- adjust_uc_em_sel(
    df_observed,
    bias_params = bp_uc_em_sel
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

# MULTINOMIAL BIAS PARAMETERS

# library(nnet)

set.seed(1234)
n <- 10000
nreps <- 10

# cont Y just for testing that function runs
df_uc_em_sel$Y_cont <- plogis(df_uc_em_sel$Y) +
  rnorm(nrow(df_uc_em_sel), mean = 0, sd = 0.1)

# 0 confounders

nobias_model <- glm(
  Y ~ X + U,
  family = binomial(link = "logit"),
  data = df_uc_em_sel_source
)

# xu_model <- multinom(
#   paste(X, U) ~ Xstar + Y,
#   data = df_uc_em_sel_source
# )

s_model <- glm(
  S ~ Xstar + Y,
  family = binomial(link = "logit"),
  data = df_uc_em_sel_source
)

df_observed <- data_observed(
  df_uc_em_sel,
  bias = c("uc", "em", "sel"),
  exposure = "Xstar",
  outcome = "Y_cont",
  confounders = NULL
)
list_for_uc_em_sel <- list(
  x1u0 = c(-1.92, 1.62, 0.70),
  x0u1 = c(-0.30, -0.01, 0.69),
  x1u1 = c(-1.64, 1.62, 1.35),
  s = as.vector(coef(s_model))
)
bp_uc_em_sel <- bias_params(coef_list = list_for_uc_em_sel)

single_run <- adjust_uc_em_sel(
  df_observed,
  bias_params = bp_uc_em_sel
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_em_sel[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    bias = c("uc", "em", "sel"),
    exposure = "Xstar",
    outcome = "Y",
    confounders = NULL
  )
  results <- adjust_uc_em_sel(
    df_observed,
    bias_params = bp_uc_em_sel
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
  Y ~ X + C1 + C2 + C3 + U,
  family = binomial(link = "logit"),
  data = df_uc_em_sel_source
)

# xu_model <- multinom(
#   paste(X, U) ~ Xstar + Y + C1 + C2 + C3,
#   data = df_uc_em_sel_source
# )

s_model <- glm(
  S ~ Xstar + Y + C1 + C2 + C3,
  family = binomial(link = "logit"),
  data = df_uc_em_sel_source
)

df_observed <- data_observed(
  df_uc_em_sel,
  bias = c("uc", "em", "sel"),
  exposure = "Xstar",
  outcome = "Y_cont",
  confounders = c("C1", "C2", "C3")
)
list_for_uc_em_sel <- list(
  x1u0 = c(-2.78, 1.62, 0.61, 0.36, -0.27, 0.88),
  x0u1 = c(-0.17, -0.01, 0.71, -0.08, 0.07, -0.15),
  x1u1 = c(-2.36, 1.62, 1.29, 0.25, -0.06, 0.74),
  s = as.vector(coef(s_model))
)
bp_uc_em_sel <- bias_params(coef_list = list_for_uc_em_sel)

single_run <- adjust_uc_em_sel(
  df_observed,
  bias_params = bp_uc_em_sel
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_em_sel[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    bias = c("uc", "em", "sel"),
    exposure = "Xstar",
    outcome = "Y",
    confounders = c("C1", "C2", "C3")
  )
  results <- adjust_uc_em_sel(
    df_observed,
    bias_params = bp_uc_em_sel
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

or_val <- adjust_uc_em_sel(
  data_observed = data_observed(
    df_uc_em_sel,
    bias = c("uc", "em", "sel"),
    exposure = "Xstar",
    outcome = "Y",
    confounders = c("C1", "C2", "C3")
  ),
  data_validation = data_validation(
    df_uc_em_sel_source,
    true_exposure = "X",
    true_outcome = "Y",
    confounders = c("C1", "C2", "C3", "U"),
    misclassified_exposure = "Xstar",
    selection = "S"
  )
)

test_that("adjust_uc_em_sel, validation data", {
  expect_gt(or_val$estimate, or_true - 0.1)
  expect_lt(or_val$estimate, or_true + 0.1)
})
