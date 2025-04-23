set.seed(1234)
nreps <- 10

# adjust with coefs
# 0 confounders

nobias_model <- glm(
  Y_bi ~ X_bi + U,
  family = binomial(link = "logit"),
  data = df_uc_source
)
or_true <- exp(summary(nobias_model)$coef[2, 1])

u_model <- glm(
  U ~ X_bi + Y_bi,
  family = binomial(link = "logit"),
  data = df_uc_source
)

df_observed <- data_observed(
  df_uc,
  bias = "uc",
  exposure = "X_bi",
  outcome = "Y_bi",
  confounders = NULL
)
list_for_uc <- list(u = as.vector(coef(u_model)))
bp_uc <- bias_params(coef_list = list_for_uc)

single_run <- multibias_adjust(
  df_observed,
  bias_params = bp_uc
)

bs_run <- multibias_adjust(
  df_observed,
  bias_params = bp_uc,
  bootstrap = TRUE,
  bootstrap_reps = nreps
)

test_that("adjust_uc, 0 confounders: OR and CI output", {
  expect_gt(bs_run$estimate, or_true - 0.1)
  expect_lt(bs_run$estimate, or_true + 0.1)
  expect_vector(
    single_run$ci,
    ptype = double(),
    size = 2
  )
})

# 3 confounders

nobias_model <- lm(Y_cont ~ X_cont + C1 + C2 + C3 + U,
  data = df_uc_source
)
or_true <- summary(nobias_model)$coef[2, 1]

u_model <- glm(U ~ X_cont + Y_cont + C1 + C2 + C3,
  family = binomial(link = "logit"),
  data = df_uc_source
)

df_observed <- data_observed(
  df_uc,
  bias = "uc",
  exposure = "X_cont",
  outcome = "Y_cont",
  confounders = c("C1", "C2", "C3")
)
list_for_uc <- list(u = as.vector(coef(u_model)))
bp_uc <- bias_params(coef_list = list_for_uc)

single_run <- adjust_uc(
  df_observed,
  bias_params = bp_uc
)

bs_run <- multibias_adjust(
  df_observed,
  bias_params = bp_uc,
  bootstrap = TRUE,
  bootstrap_reps = nreps
)

test_that("adjust_uc, 3 confounders: OR and CI output", {
  expect_gt(bs_run$estimate, or_true - 0.1)
  expect_lt(bs_run$estimate, or_true + 0.1)
  expect_vector(
    single_run$ci,
    ptype = double(),
    size = 2
  )
})

# adjust with validation data

or_val <- adjust_uc(
  data_observed = data_observed(
    df_uc,
    bias = "uc",
    exposure = "X_cont",
    outcome = "Y_cont",
    confounders = c("C1", "C2", "C3")
  ),
  data_validation = data_validation(
    df_uc_source,
    true_exposure = "X_cont",
    true_outcome = "Y_cont",
    confounders = c("C1", "C2", "C3", "U")
  )
)

test_that("adjust_uc, validation data", {
  expect_gt(or_val$estimate, or_true - 0.1)
  expect_lt(or_val$estimate, or_true + 0.1)
})