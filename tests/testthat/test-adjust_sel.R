set.seed(1234)
nreps <- 10

s_model <- glm(S ~ X + Y,
  family = binomial(link = "logit"),
  data = df_sel_source
)

# 0 confounders

# nobias_model <- glm(
#   Y ~ X,
#   family = binomial(link = "logit"),
#   data = df_sel_source
# )
# or_true <- exp(summary(nobias_model)$coef[2, 1])

# df_observed <- data_observed(
#   df_sel,
#   bias = "sel",
#   exposure = "X",
#   outcome = "Y",
#   confounders = NULL
# )
# list_for_sel <- list(s = as.vector(coef(s_model)))
# bp_sel <- bias_params(coef_list = list_for_sel)

# single_run <- multibias_adjust(
#   df_observed,
#   bias_params = bp_sel
# )

# bs_run <- multibias_adjust(
#   df_observed,
#   bias_params = bp_sel,
#   bootstrap = TRUE,
#   bootstrap_reps = nreps
# )

# test_that("odds ratio and confidence interval output", {
#   expect_gt(bs_run$estimate, or_true - 0.1)
#   expect_lt(bs_run$estimate, or_true + 0.1)
#   expect_vector(
#     single_run$ci,
#     ptype = double(),
#     size = 2
#   )
# })

# 3 confounders

nobias_model <- glm(Y ~ X + C1 + C2 + C3,
  family = binomial(link = "logit"),
  data = df_sel_source
)
or_true <- exp(summary(nobias_model)$coef[2, 1])

df_observed <- data_observed(
  df_sel,
  bias = "sel",
  exposure = "X",
  outcome = "Y",
  confounders = c("C1", "C2", "C3")
)
list_for_sel <- list(s = as.vector(coef(s_model)))
bp_sel <- bias_params(coef_list = list_for_sel)

single_run <- multibias_adjust(
  df_observed,
  bias_params = bp_sel
)

bs_run <- multibias_adjust(
  df_observed,
  bias_params = bp_sel,
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
    df_sel,
    bias = "sel",
    exposure = "X",
    outcome = "Y",
    confounders = c("C1", "C2", "C3")
  ),
  data_validation = data_validation(
    df_sel_source,
    true_exposure = "X",
    true_outcome = "Y",
    confounders = c("C1", "C2", "C3"),
    selection = "S"
  ),
  bootstrap = TRUE,
  bootstrap_reps = nreps
)

test_that("adjust_sel, validation data", {
  expect_gt(val_run$estimate, or_true - 0.1)
  expect_lt(val_run$estimate, or_true + 0.1)
})
