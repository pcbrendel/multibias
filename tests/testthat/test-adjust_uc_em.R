set.seed(1234)
nreps <- 10

# SEPARATE BIAS PARAMETERS
# 0 confounders

nobias_model <- glm(
  Y ~ X + U,
  family = binomial(link = "logit"),
  data = df_uc_em_source
)
or_true <- exp(summary(nobias_model)$coef[2, 1])

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

single_run <- multibias_adjust(
  df_observed,
  bias_params = bp_uc_em
)

bs_run <- multibias_adjust(
  df_observed,
  bias_params = bp_uc_em,
  bootstrap = TRUE,
  bootstrap_reps = nreps
)

test_that("odds ratio and confidence interval output", {
  expect_gt(bs_run$estimate, or_true - 0.1)
  expect_lt(bs_run$estimate, or_true + 0.1)
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
or_true <- exp(summary(nobias_model)$coef[2, 1])

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

single_run <- multibias_adjust(
  df_observed,
  bias_params = bp_uc_em
)

bs_run <- multibias_adjust(
  df_observed,
  bias_params = bp_uc_em,
  bootstrap = TRUE,
  bootstrap_reps = nreps
)

test_that("odds ratio and confidence interval output", {
  expect_gt(bs_run$estimate, or_true - 0.1)
  expect_lt(bs_run$estimate, or_true + 0.1)
  expect_vector(
    single_run$ci,
    ptype = double(),
    size = 2
  )
})

# MULTINOMIAL BIAS PARAMETERS
# 0 confounders

nobias_model <- glm(
  Y ~ X + U,
  family = binomial(link = "logit"),
  data = df_uc_em_source
)
or_true <- exp(summary(nobias_model)$coef[2, 1])

# library(nnet)
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

single_run <- multibias_adjust(
  df_observed,
  bias_params = bp_uc_em
)

bs_run <- multibias_adjust(
  df_observed,
  bias_params = bp_uc_em,
  bootstrap = TRUE,
  bootstrap_reps = nreps
)

test_that("odds ratio and confidence interval output", {
  expect_gt(bs_run$estimate, or_true - 0.1)
  expect_lt(bs_run$estimate, or_true + 0.1)
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
or_true <- exp(summary(nobias_model)$coef[2, 1])

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

single_run <- multibias_adjust(
  df_observed,
  bias_params = bp_uc_em
)

bs_run <- multibias_adjust(
  df_observed,
  bias_params = bp_uc_em,
  bootstrap = TRUE,
  bootstrap_reps = nreps
)

test_that("odds ratio and confidence interval output", {
  expect_gt(bs_run$estimate, or_true - 0.1)
  expect_lt(bs_run$estimate, or_true + 0.1)
  expect_vector(
    single_run$ci,
    ptype = double(),
    size = 2
  )
})

# adjust with validation data

val_run <- multibias_adjust(
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
  ),
  bootstrap = TRUE,
  bootstrap_reps = nreps
)

test_that("adjust_uc_em, validation data", {
  expect_gt(val_run$estimate, or_true - 0.1)
  expect_lt(val_run$estimate, or_true + 0.1)
})
