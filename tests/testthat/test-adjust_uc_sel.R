set.seed(1234)
n <- 50000
nreps <- 10

# cont Y just for testing that function runs
df_uc_sel$Y_cont <- plogis(df_uc_sel$Y) + rnorm(nrow(df_uc_sel), mean = 0, sd = 0.1)

# 0 confounders

nobias_model <- glm(
  Y ~ X + U,
  family = binomial(link = "logit"),
  data = df_uc_sel_source
)

u_model <- glm(
  U ~ X + Y,
  family = binomial(link = "logit"),
  data = df_uc_sel_source
)

s_model <- glm(
  S ~ X + Y,
  family = binomial(link = "logit"),
  data = df_uc_sel_source
)

df_observed <- data_observed(
  df_uc_sel,
  bias = c("uc", "sel"),
  exposure = "X",
  outcome = "Y",
  confounders = NULL
)
list_for_uc_sel <- list(
  u = as.vector(coef(u_model)),
  s = as.vector(coef(s_model))
)
bp_uc_sel <- bias_params(coef_list = list_for_uc_sel)

single_run <- adjust_uc_sel(
  df_observed,
  bias_params = bp_uc_sel
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_sel[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    bias = c("uc", "sel"),
    exposure = "X",
    outcome = "Y",
    confounders = NULL
  )
  results <- adjust_uc_sel(
    df_observed,
    bias_params = bp_uc_sel
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
  Y ~ X + C1 + C2 + C3 + U,
  family = binomial(link = "logit"),
  data = df_uc_sel_source
)

u_model <- glm(
  U ~ X + Y + C1 + C2 + C3,
  family = binomial(link = "logit"),
  data = df_uc_sel_source
)

s_model <- glm(
  S ~ X + Y,
  family = binomial(link = "logit"),
  data = df_uc_sel_source
)

df_observed <- data_observed(
  df_uc_sel,
  bias = c("uc", "sel"),
  exposure = "X",
  outcome = "Y",
  confounders = c("C1", "C2", "C3")
)
list_for_uc_sel <- list(
  u = as.vector(coef(u_model)),
  s = as.vector(coef(s_model))
)
bp_uc_sel <- bias_params(coef_list = list_for_uc_sel)

single_run <- adjust_uc_sel(
  df_observed,
  bias_params = bp_uc_sel
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_sel[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    bias = c("uc", "sel"),
    exposure = "X",
    outcome = "Y",
    confounders = c("C1", "C2", "C3")
  )
  results <- adjust_uc_sel(
    df_observed,
    bias_params = bp_uc_sel
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

or_val <- adjust_uc_sel(
  data_observed = data_observed(
    df_uc_sel,
    bias = c("uc", "sel"),
    exposure = "X",
    outcome = "Y",
    confounders = c("C1", "C2", "C3")
  ),
  data_validation = data_validation(
    df_uc_sel_source,
    true_exposure = "X",
    true_outcome = "Y",
    confounders = c("C1", "C2", "C3", "U"),
    selection = "S"
  )
)

test_that("adjust_uc_sel, validation data", {
  expect_gt(or_val$estimate, or_true - 0.1)
  expect_lt(or_val$estimate, or_true + 0.1)
})