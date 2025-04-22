set.seed(1234)
nreps <- 10

# 0 confounders

nobias_model <- glm(
  Y ~ X,
  family = binomial(link = "logit"),
  data = df_om_source
)
or_true <- exp(summary(nobias_model)$coef[2, 1])

y_model <- glm(
  Y ~ X + Ystar,
  family = binomial(link = "logit"),
  data = df_om_source
)

df_observed <- data_observed(
  df_om,
  bias = "om",
  exposure = "X",
  outcome = "Ystar",
  confounders = NULL
)
list_for_om <- list(y = as.vector(coef(y_model)))
bp_om <- bias_params(coef_list = list_for_om)

single_run <- adjust_om(
  df_observed,
  bias_params = bp_om
)

bs_run <- multibias_adjust(
  df_observed,
  bias_params = bp_om,
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

nobias_model <- glm(Y ~ X + C1 + C2 + C3,
  family = binomial(link = "logit"),
  data = df_om_source
)
or_true <- exp(summary(nobias_model)$coef[2, 1])

y_model <- glm(Y ~ X + Ystar + C1 + C2 + C3,
  family = binomial(link = "logit"),
  data = df_om_source
)

df_observed <- data_observed(
  df_om,
  bias = "om",
  exposure = "X",
  outcome = "Ystar",
  confounders = c("C1", "C2", "C3")
)
list_for_om <- list(y = as.vector(coef(y_model)))
bp_om <- bias_params(coef_list = list_for_om)

single_run <- adjust_om(
  df_observed,
  bias_params = bp_om
)

bs_run <- multibias_adjust(
  df_observed,
  bias_params = bp_om,
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

or_val <- adjust_om(
  data_observed = data_observed(
    df_om,
    bias = "om",
    exposure = "X",
    outcome = "Ystar",
    confounders = c("C1", "C2", "C3")
  ),
  data_validation = data_validation(
    df_om_source,
    true_exposure = "X",
    true_outcome = "Y",
    confounders = c("C1", "C2", "C3"),
    misclassified_outcome = "Ystar"
  )
)

test_that("adjust_om, validation data", {
  expect_gt(or_val$estimate, or_true - 0.1)
  expect_lt(or_val$estimate, or_true + 0.1)
})
