set.seed(1234)
n <- 10000
nreps <- 10

# cont X just for testing that function runs
df_om_sel$X_cont <- plogis(df_om_sel$X) +
  rnorm(nrow(df_om_sel), mean = 0, sd = 0.1)

# 0 confounders

nobias_model <- glm(
  Y ~ X,
  family = binomial(link = "logit"),
  data = df_om_sel_source
)

y_model <- glm(
  Y ~ X + Ystar,
  family = binomial(link = "logit"),
  data = df_om_sel_source
)

s_model <- glm(
  S ~ X + Ystar,
  family = binomial(link = "logit"),
  data = df_om_sel_source
)

df_observed <- data_observed(
  df_om_sel,
  bias = c("om", "sel"),
  exposure = "X_cont",
  outcome = "Ystar",
  confounders = NULL
)
list_for_om_sel <- list(
  y = as.vector(coef(y_model)),
  s = as.vector(coef(s_model))
)
bp_om_sel <- bias_params(coef_list = list_for_om_sel)

single_run <- adjust_om_sel(
  df_observed,
  bias_params = bp_om_sel
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_om_sel[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    bias = c("om", "sel"),
    exposure = "X",
    outcome = "Ystar",
    confounders = NULL
  )
  results <- adjust_om_sel(
    df_observed,
    bias_params = bp_om_sel
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
  data = df_om_sel_source
)

y_model <- glm(
  Y ~ X + Ystar + C1 + C2 + C3,
  family = binomial(link = "logit"),
  data = df_om_sel_source
)

s_model <- glm(
  S ~ X + Ystar + C1 + C2 + C3,
  family = binomial(link = "logit"),
  data = df_om_sel_source
)

df_observed <- data_observed(
  df_om_sel,
  bias = c("om", "sel"),
  exposure = "X_cont",
  outcome = "Ystar",
  confounders = c("C1", "C2", "C3")
)
list_for_om_sel <- list(
  y = as.vector(coef(y_model)),
  s = as.vector(coef(s_model))
)
bp_om_sel <- bias_params(coef_list = list_for_om_sel)

single_run <- adjust_om_sel(
  df_observed,
  bias_params = bp_om_sel
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_om_sel[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    bias = c("om", "sel"),
    exposure = "X",
    outcome = "Ystar",
    confounders = c("C1", "C2", "C3")
  )
  results <- adjust_om_sel(
    df_observed,
    bias_params = bp_om_sel
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

or_val <- adjust_om_sel(
  data_observed = data_observed(
    df_om_sel,
    bias = c("om", "sel"),
    exposure = "X",
    outcome = "Ystar",
    confounders = c("C1", "C2", "C3")
  ),
  data_validation = data_validation(
    df_om_sel_source,
    true_exposure = "X",
    true_outcome = "Y",
    confounders = c("C1", "C2", "C3"),
    misclassified_outcome = "Ystar",
    selection = "S"
  )
)

test_that("adjust_om_sel, validation data", {
  expect_gt(or_val$estimate, or_true - 0.1)
  expect_lt(or_val$estimate, or_true + 0.1)
})