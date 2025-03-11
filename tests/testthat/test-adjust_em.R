set.seed(1234)
n <- 20000
nreps <- 10

# cont Y just for testing that function runs
df_em$Y_cont <- plogis(df_em$Y) + rnorm(nrow(df_em), mean = 0, sd = 0.1)

# 0 confounders

nobias_model <- glm(
  Y ~ X,
  family = binomial(link = "logit"),
  data = df_em_source
)

x_model <- glm(
  X ~ Xstar + Y,
  family = binomial(link = "logit"),
  data = df_em_source
)

df_observed <- data_observed(
  df_em,
  bias = "em",
  exposure = "Xstar",
  outcome = "Y_cont",
  confounders = NULL
)
list_for_em <- list(x = as.vector(coef(x_model)))
bp_em <- bias_params(coef_list = list_for_em)

single_run <- adjust_em(
  df_observed,
  bias_params = bp_em
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_em[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    bias = "em",
    exposure = "Xstar",
    outcome = "Y",
    confounders = NULL
  )
  results <- adjust_em(
    df_observed,
    bias_params = bp_em
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
  data = df_em_source
)

x_model <- glm(
  X ~ Xstar + Y + C1 + C2 + C3,
  family = binomial(link = "logit"),
  data = df_em_source
)

df_observed <- data_observed(
  df_em,
  bias = "em",
  exposure = "Xstar",
  outcome = "Y_cont",
  confounders = c("C1", "C2", "C3")
)
list_for_em <- list(x = as.vector(coef(x_model)))
bp_em <- bias_params(coef_list = list_for_em)

single_run <- adjust_em(
  df_observed,
  bias_params = bp_em
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_em[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    bias = "em",
    exposure = "Xstar",
    outcome = "Y",
    confounders = c("C1", "C2", "C3")
  )
  results <- adjust_em(
    df_observed,
    bias_params = bp_em
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

or_val <- adjust_em(
  data_observed = data_observed(
    df_em,
    bias = "em",
    exposure = "Xstar",
    outcome = "Y",
    confounders = c("C1", "C2", "C3")
  ),
  data_validation = data_validation(
    df_em_source,
    true_exposure = "X",
    true_outcome = "Y",
    confounders = c("C1", "C2", "C3"),
    misclassified_exposure = "Xstar"
  )
)

test_that("adjust_em, validation data", {
  expect_gt(or_val$estimate, or_true - 0.1)
  expect_lt(or_val$estimate, or_true + 0.1)
})
