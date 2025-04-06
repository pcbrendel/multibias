set.seed(1234)
n <- 50000
nreps <- 10

# cont Y just for testing that function runs
df_em_sel$Y_cont <- plogis(df_em_sel$Y) +
  rnorm(nrow(df_em_sel), mean = 0, sd = 0.1)

# 0 confounders

nobias_model <- glm(
  Y ~ X,
  family = binomial(link = "logit"),
  data = df_em_sel_source
)

x_model <- glm(
  X ~ Xstar + Y,
  family = binomial(link = "logit"),
  data = df_em_sel_source
)

s_model <- glm(
  S ~ Xstar + Y,
  family = binomial(link = "logit"),
  data = df_em_sel_source
)

df_observed <- data_observed(
  df_em_sel,
  bias = c("em", "sel"),
  exposure = "Xstar",
  outcome = "Y_cont",
  confounders = NULL
)
list_for_em_sel <- list(
  x = as.vector(coef(x_model)),
  s = as.vector(coef(s_model))
)
bp_em_sel <- bias_params(coef_list = list_for_em_sel)


single_run <- adjust_em_sel(
  df_observed,
  bias_params = bp_em_sel
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_em_sel[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    bias = c("em", "sel"),
    exposure = "Xstar",
    outcome = "Y",
    confounders = NULL
  )
  results <- adjust_em_sel(
    df_observed,
    bias_params = bp_em_sel
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

nobias_model <- glm(Y ~ X + C1 + C2 + C3,
  family = binomial(link = "logit"),
  data = df_em_sel_source
)

x_model <- glm(X ~ Xstar + Y + C1 + C2 + C3,
  family = binomial(link = "logit"),
  data = df_em_sel_source
)

s_model <- glm(S ~ Xstar + Y + C1 + C2 + C3,
  family = binomial(link = "logit"),
  data = df_em_sel_source
)

df_observed <- data_observed(
  df_em_sel,
  bias = c("em", "sel"),
  exposure = "Xstar",
  outcome = "Y_cont",
  confounders = c("C1", "C2", "C3")
)
list_for_em_sel <- list(
  x = as.vector(coef(x_model)),
  s = as.vector(coef(s_model))
)
bp_em_sel <- bias_params(coef_list = list_for_em_sel)

single_run <- adjust_em_sel(
  df_observed,
  bias_params = bp_em_sel
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_em_sel[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    bias = c("em", "sel"),
    exposure = "Xstar",
    outcome = "Y",
    confounders = c("C1", "C2", "C3")
  )
  results <- adjust_em_sel(
    df_observed,
    bias_params = bp_em_sel
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

or_val <- adjust_em_sel(
  data_observed = data_observed(
    df_em_sel,
    bias = c("em", "sel"),
    exposure = "Xstar",
    outcome = "Y",
    confounders = c("C1", "C2", "C3")
  ),
  data_validation = data_validation(
    df_em_sel_source,
    true_exposure = "X",
    true_outcome = "Y",
    confounders = c("C1", "C2", "C3"),
    misclassified_exposure = "Xstar",
    selection = "S"
  )
)

test_that("adjust_em_sel, validation data", {
  expect_gt(or_val$estimate, or_true - 0.1)
  expect_lt(or_val$estimate, or_true + 0.1)
})
