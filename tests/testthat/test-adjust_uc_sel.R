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
  exposure = "X",
  outcome = "Y_cont",
  confounders = NULL
)

single_run <- adjust_uc_sel(
  df_observed,
  u_model_coefs = c(
    u_model$coef[1],
    u_model$coef[2],
    u_model$coef[3]
  ),
  s_model_coefs = c(
    s_model$coef[1],
    s_model$coef[2],
    s_model$coef[3]
  )
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_sel[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    exposure = "X",
    outcome = "Y",
    confounders = NULL
  )
  results <- adjust_uc_sel(
    df_observed,
    u_model_coefs = c(
      u_model$coef[1],
      u_model$coef[2],
      u_model$coef[3]
    ),
    s_model_coefs = c(
      s_model$coef[1],
      s_model$coef[2],
      s_model$coef[3]
    )
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

# 1 confounder

nobias_model <- glm(
  Y ~ X + C1 + U,
  family = binomial(link = "logit"),
  data = df_uc_sel_source
)

u_model <- glm(
  U ~ X + Y + C1,
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
  exposure = "X",
  outcome = "Y_cont",
  confounders = "C1"
)

single_run <- adjust_uc_sel(
  df_observed,
  u_model_coefs = c(
    u_model$coef[1],
    u_model$coef[2],
    u_model$coef[3],
    u_model$coef[4]
  ),
  s_model_coefs = c(
    s_model$coef[1],
    s_model$coef[2],
    s_model$coef[3]
  )
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_sel[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    exposure = "X",
    outcome = "Y",
    confounders = "C1"
  )
  results <- adjust_uc_sel(
    df_observed,
    u_model_coefs = c(
      u_model$coef[1],
      u_model$coef[2],
      u_model$coef[3],
      u_model$coef[4]
    ),
    s_model_coefs = c(
      s_model$coef[1],
      s_model$coef[2],
      s_model$coef[3]
    )
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

# 2 confounders

nobias_model <- glm(
  Y ~ X + C1 + C2 + U,
  family = binomial(link = "logit"),
  data = df_uc_sel_source
)

u_model <- glm(
  U ~ X + Y + C1 + C2,
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
  exposure = "X",
  outcome = "Y_cont",
  confounders = c("C1", "C2")
)

single_run <- adjust_uc_sel(
  df_observed,
  u_model_coefs = c(
    u_model$coef[1],
    u_model$coef[2],
    u_model$coef[3],
    u_model$coef[4],
    u_model$coef[5]
  ),
  s_model_coefs = c(
    s_model$coef[1],
    s_model$coef[2],
    s_model$coef[3]
  )
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_sel[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    exposure = "X",
    outcome = "Y",
    confounders = c("C1", "C2")
  )
  results <- adjust_uc_sel(
    df_observed,
    u_model_coefs = c(
      u_model$coef[1],
      u_model$coef[2],
      u_model$coef[3],
      u_model$coef[4],
      u_model$coef[5]
    ),
    s_model_coefs = c(
      s_model$coef[1],
      s_model$coef[2],
      s_model$coef[3]
    )
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
  exposure = "X",
  outcome = "Y_cont",
  confounders = c("C1", "C2", "C3")
)

single_run <- adjust_uc_sel(
  df_observed,
  u_model_coefs = c(
    u_model$coef[1],
    u_model$coef[2],
    u_model$coef[3],
    u_model$coef[4],
    u_model$coef[5],
    u_model$coef[6]
  ),
  s_model_coefs = c(
    s_model$coef[1],
    s_model$coef[2],
    s_model$coef[3]
  )
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_sel[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    exposure = "X",
    outcome = "Y",
    confounders = c("C1", "C2", "C3")
  )
  results <- adjust_uc_sel(
    df_observed,
    u_model_coefs = c(
      u_model$coef[1],
      u_model$coef[2],
      u_model$coef[3],
      u_model$coef[4],
      u_model$coef[5],
      u_model$coef[6]
    ),
    s_model_coefs = c(
      s_model$coef[1],
      s_model$coef[2],
      s_model$coef[3]
    )
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