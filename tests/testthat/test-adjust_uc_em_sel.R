set.seed(1234)
nreps <- 10

# SEPARATE BIAS PARAMETERS
# 0 confounders

nobias_model <- glm(
  Y ~ X + U,
  family = binomial(link = "logit"),
  data = df_uc_em_sel_source
)
or_true <- exp(summary(nobias_model)$coef[2, 1])

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
  outcome = "Y",
  confounders = NULL
)
list_for_uc_em_sel <- list(
  u = as.vector(coef(u_model)),
  x = as.vector(coef(x_model)),
  s = as.vector(coef(s_model))
)
bp_uc_em_sel <- bias_params(coef_list = list_for_uc_em_sel)

single_run <- multibias_adjust(
  df_observed,
  bias_params = bp_uc_em_sel
)

bs_run <- multibias_adjust(
  df_observed,
  bias_params = bp_uc_em_sel,
  bootstrap = TRUE,
  bootstrap_reps = nreps
)

test_that("0 confounders: odds ratio and confidence interval output", {
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
  data = df_uc_em_sel_source
)
or_true <- exp(summary(nobias_model)$coef[2, 1])

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
  outcome = "Y",
  confounders = c("C1", "C2", "C3")
)
list_for_uc_em_sel <- list(
  u = as.vector(coef(u_model)),
  x = as.vector(coef(x_model)),
  s = as.vector(coef(s_model))
)
bp_uc_em_sel <- bias_params(coef_list = list_for_uc_em_sel)

single_run <- multibias_adjust(
  df_observed,
  bias_params = bp_uc_em_sel
)

bs_run <- multibias_adjust(
  df_observed,
  bias_params = bp_uc_em_sel,
  bootstrap = TRUE,
  bootstrap_reps = nreps
)

test_that("3 confounders: odds ratio and confidence interval output", {
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
  data = df_uc_em_sel_source
)
or_true <- exp(summary(nobias_model)$coef[2, 1])

# library(nnet)
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
  outcome = "Y",
  confounders = NULL
)
list_for_uc_em_sel <- list(
  x1u0 = c(-1.92, 1.62, 0.70),
  x0u1 = c(-0.30, -0.01, 0.69),
  x1u1 = c(-1.64, 1.62, 1.35),
  s = as.vector(coef(s_model))
)
bp_uc_em_sel <- bias_params(coef_list = list_for_uc_em_sel)

single_run <- multibias_adjust(
  df_observed,
  bias_params = bp_uc_em_sel
)

bs_run <- multibias_adjust(
  df_observed,
  bias_params = bp_uc_em_sel,
  bootstrap = TRUE,
  bootstrap_reps = nreps
)

test_that("0 confounders: odds ratio and confidence interval output", {
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
  data = df_uc_em_sel_source
)
or_true <- exp(summary(nobias_model)$coef[2, 1])

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
  outcome = "Y",
  confounders = c("C1", "C2", "C3")
)
list_for_uc_em_sel <- list(
  x1u0 = c(-2.78, 1.62, 0.61, 0.36, -0.27, 0.88),
  x0u1 = c(-0.17, -0.01, 0.71, -0.08, 0.07, -0.15),
  x1u1 = c(-2.36, 1.62, 1.29, 0.25, -0.06, 0.74),
  s = as.vector(coef(s_model))
)
bp_uc_em_sel <- bias_params(coef_list = list_for_uc_em_sel)

single_run <- multibias_adjust(
  df_observed,
  bias_params = bp_uc_em_sel
)

bs_run <- multibias_adjust(
  df_observed,
  bias_params = bp_uc_em_sel,
  bootstrap = TRUE,
  bootstrap_reps = nreps
)

test_that("3 confounders: odds ratio and confidence interval output", {
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
  ),
  bootstrap = TRUE,
  bootstrap_reps = nreps
)

test_that("adjust_uc_em_sel, validation data", {
  expect_gt(val_run$estimate, or_true - 0.1)
  expect_lt(val_run$estimate, or_true + 0.1)
})
