set.seed(1234)
n <- 10000
nreps <- 10

# SEPARATE BIAS PARAMETERS

# cont X just for testing that function runs
df_uc_om_sel$X_cont <- plogis(df_uc_om_sel$X) +
  rnorm(nrow(df_uc_om_sel), mean = 0, sd = 0.1)

# 0 confounders

nobias_model <- glm(
  Y ~ X + U,
  family = binomial(link = "logit"),
  data = df_uc_om_sel_source
)

u_model <- glm(
  U ~ X + Y,
  family = binomial(link = "logit"),
  data = df_uc_om_sel_source
)
y_model <- glm(
  Y ~ X + Ystar,
  family = binomial(link = "logit"),
  data = df_uc_om_sel_source
)
s_model <- glm(
  S ~ X + Ystar,
  family = binomial(link = "logit"),
  data = df_uc_om_sel_source
)

df_observed <- data_observed(
  df_uc_om_sel,
  bias = c("uc", "om", "sel"),
  exposure = "X_cont",
  outcome = "Ystar",
  confounders = NULL
)
list_for_uc_om_sel <- list(
  u = as.vector(coef(u_model)),
  y = as.vector(coef(y_model)),
  s = as.vector(coef(s_model))
)
bp_uc_om_sel <- bias_params(coef_list = list_for_uc_om_sel)

single_run <- adjust_uc_om_sel(
  df_observed,
  bias_params = bp_uc_om_sel
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_om_sel[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    bias = c("uc", "om", "sel"),
    exposure = "X",
    outcome = "Ystar",
    confounders = NULL
  )
  results <- adjust_uc_om_sel(
    df_observed,
    bias_params = bp_uc_om_sel
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
  data = df_uc_om_sel_source
)

u_model <- glm(
  U ~ X + Y,
  family = binomial(link = "logit"),
  data = df_uc_om_sel_source
)
y_model <- glm(
  Y ~ X + Ystar + C1 + C2 + C3,
  family = binomial(link = "logit"),
  data = df_uc_om_sel_source
)
s_model <- glm(
  S ~ X + Ystar + C1 + C2 + C3,
  family = binomial(link = "logit"),
  data = df_uc_om_sel_source
)

df_observed <- data_observed(
  df_uc_om_sel,
  bias = c("uc", "om", "sel"),
  exposure = "X_cont",
  outcome = "Ystar",
  confounders = c("C1", "C2", "C3")
)
list_for_uc_om_sel <- list(
  u = as.vector(coef(u_model)),
  y = as.vector(coef(y_model)),
  s = as.vector(coef(s_model))
)
bp_uc_om_sel <- bias_params(coef_list = list_for_uc_om_sel)

single_run <- adjust_uc_om_sel(
  df_observed,
  bias_params = bp_uc_om_sel
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_om_sel[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    bias = c("uc", "om", "sel"),
    exposure = "X",
    outcome = "Ystar",
    confounders = c("C1", "C2", "C3")
  )
  results <- adjust_uc_om_sel(
    df_observed,
    bias_params = bp_uc_om_sel
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

# cont X just for testing that function runs
df_uc_om_sel$X_cont <- plogis(df_uc_om_sel$X) +
  rnorm(nrow(df_uc_om_sel), mean = 0, sd = 0.1)

# 0 confounders

nobias_model <- glm(
  Y ~ X + U,
  family = binomial(link = "logit"),
  data = df_uc_om_sel_source
)

# uy_model <- multinom(
#   paste(U, Y) ~ X + Ystar,
#   data = df_uc_om_sel_source
# )

s_model <- glm(
  S ~ X + Ystar,
  family = binomial(link = "logit"),
  data = df_uc_om_sel_source
)

df_observed <- data_observed(
  df_uc_om_sel,
  bias = c("uc", "om", "sel"),
  exposure = "X",
  outcome = "Ystar",
  confounders = NULL
)
list_for_uc_om_sel <- list(
  u1y0 = c(-0.32, 0.59, 0.01),
  u0y1 = c(-2.98, 0.71, 1.65),
  u1y1 = c(-2.59, 1.27, 1.64),
  s = as.vector(coef(s_model))
)
bp_uc_om_sel <- bias_params(coef_list = list_for_uc_om_sel)

single_run <- adjust_uc_om_sel(
  df_observed,
  bias_params = bp_uc_om_sel
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_om_sel[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    bias = c("uc", "om", "sel"),
    exposure = "X",
    outcome = "Ystar",
    confounders = NULL
  )
  results <- adjust_uc_om_sel(
    df_observed,
    bias_params = bp_uc_om_sel
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
  data = df_uc_om_sel_source
)

# uy_model <- multinom(
#   paste(U, Y) ~ X + Ystar + C1 + C2 + C3,
#   data = df_uc_om_sel_source
# )

s_model <- glm(
  S ~ X + Ystar + C1 + C2 + C3,
  family = binomial(link = "logit"),
  data = df_uc_om_sel_source
)

df_observed <- data_observed(
  df_uc_om_sel,
  bias = c("uc", "om", "sel"),
  exposure = "X",
  outcome = "Ystar",
  confounders = c("C1", "C2", "C3")
)
list_for_uc_om_sel <- list(
  u1y0 = c(-0.20, 0.62, 0.01, -0.08, 0.10, -0.15),
  u0y1 = c(-3.28, 0.63, 1.65, 0.42, -0.85, 0.26),
  u1y1 = c(-2.70, 1.22, 1.64, 0.32, -0.77, 0.09),
  s = as.vector(coef(s_model))
)
bp_uc_om_sel <- bias_params(coef_list = list_for_uc_om_sel)

single_run <- adjust_uc_om_sel(
  df_observed,
  bias_params = bp_uc_om_sel
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_om_sel[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    bias = c("uc", "om", "sel"),
    exposure = "X",
    outcome = "Ystar",
    confounders = c("C1", "C2", "C3")
  )
  results <- adjust_uc_om_sel(
    df_observed,
    bias_params = bp_uc_om_sel
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

or_val <- adjust_uc_om_sel(
  data_observed = data_observed(
    df_uc_om_sel,
    bias = c("uc", "om", "sel"),
    exposure = "X",
    outcome = "Ystar",
    confounders = c("C1", "C2", "C3")
  ),
  data_validation = data_validation(
    df_uc_om_sel_source,
    true_exposure = "X",
    true_outcome = "Y",
    confounders = c("C1", "C2", "C3", "U"),
    misclassified_outcome = "Ystar",
    selection = "S"
  )
)

test_that("adjust_uc_om_sel, validation data", {
  expect_gt(or_val$estimate, or_true - 0.1)
  expect_lt(or_val$estimate, or_true + 0.1)
})