set.seed(1234)
n <- 20000
nreps <- 10

# 0 confounders

nobias_model <- glm(
  Y_bi ~ X_bi + U,
  family = binomial(link = "logit"),
  data = df_uc_source
)

u_model <- glm(
  U ~ X_bi + Y_bi,
  family = binomial(link = "logit"),
  data = df_uc_source
)

df_observed <- data_observed(
  df_uc,
  exposure = "X_bi",
  outcome = "Y_bi",
  confounders = NULL
)

single_run <- adjust_uc(
  df_observed,
  u_model_coefs = c(
    u_model$coef[1],
    u_model$coef[2],
    u_model$coef[3]
  )
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    exposure = "X_bi",
    outcome = "Y_bi",
    confounders = NULL
  )
  results <- adjust_uc(
    df_observed,
    u_model_coefs = c(
      u_model$coef[1],
      u_model$coef[2],
      u_model$coef[3]
    )
  )
  est[i] <- results$estimate
}

or_true <- exp(summary(nobias_model)$coef[2, 1])
or_adjusted <- median(est)

test_that("adjust_uc, 0 confounders: OR and CI output", {
  expect_gt(or_adjusted, or_true - 0.1)
  expect_lt(or_adjusted, or_true + 0.1)
  expect_vector(
    single_run$ci,
    ptype = double(),
    size = 2
  )
})

# 1 confounder

nobias_model <- glm(Y_bi ~ X_bi + C1 + U,
  family = binomial(link = "logit"),
  data = df_uc_source
)

u_model <- glm(U ~ X_bi + Y_bi + C1,
  family = binomial(link = "logit"),
  data = df_uc_source
)

df_observed <- data_observed(
  df_uc,
  exposure = "X_bi",
  outcome = "Y_bi",
  confounders = "C1"
)

single_run <- adjust_uc(
  df_observed,
  u_model_coefs = c(
    u_model$coef[1],
    u_model$coef[2],
    u_model$coef[3],
    u_model$coef[4]
  )
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    exposure = "X_bi",
    outcome = "Y_bi",
    confounders = "C1"
  )
  results <- adjust_uc(
    df_observed,
    u_model_coefs = c(
      u_model$coef[1],
      u_model$coef[2],
      u_model$coef[3],
      u_model$coef[4]
    )
  )
  est[i] <- results$estimate
}

or_true <- exp(summary(nobias_model)$coef[2, 1])
or_adjusted <- median(est)

test_that("adjust_uc, 1 confounder: OR and CI output", {
  expect_gt(or_adjusted, or_true - 0.1)
  expect_lt(or_adjusted, or_true + 0.1)
  expect_vector(
    single_run$ci,
    ptype = double(),
    size = 2
  )
})

# 2 confounders

nobias_model <- lm(Y_cont ~ X_cont + C1 + C2 + U, data = df_uc_source)

u_model <- glm(U ~ X_cont + Y_cont + C1 + C2,
  family = binomial(link = "logit"),
  data = df_uc_source
)

df_observed <- data_observed(
  df_uc,
  exposure = "X_cont",
  outcome = "Y_cont",
  confounders = c("C1", "C2")
)

single_run <- adjust_uc(
  df_observed,
  u_model_coefs = c(
    u_model$coef[1],
    u_model$coef[2],
    u_model$coef[3],
    u_model$coef[4],
    u_model$coef[5]
  )
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    exposure = "X_cont",
    outcome = "Y_cont",
    confounders = c("C1", "C2")
  )
  results <- adjust_uc(
    df_observed,
    u_model_coefs = c(
      u_model$coef[1],
      u_model$coef[2],
      u_model$coef[3],
      u_model$coef[4],
      u_model$coef[5]
    )
  )
  est[i] <- results$estimate
}

or_true <- summary(nobias_model)$coef[2, 1]
or_adjusted <- median(est)

test_that("adjust_uc, 2 confounders: OR and CI output", {
  expect_gt(or_adjusted, or_true - 0.1)
  expect_lt(or_adjusted, or_true + 0.1)
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

u_model <- glm(U ~ X_cont + Y_cont + C1 + C2 + C3,
  family = binomial(link = "logit"),
  data = df_uc_source
)

df_observed <- data_observed(
  df_uc,
  exposure = "X_cont",
  outcome = "Y_cont",
  confounders = c("C1", "C2", "C3")
)

single_run <- adjust_uc(
  df_observed,
  u_model_coefs = c(
    u_model$coef[1],
    u_model$coef[2],
    u_model$coef[3],
    u_model$coef[4],
    u_model$coef[5],
    u_model$coef[6]
  )
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    exposure = "X_cont",
    outcome = "Y_cont",
    confounders = c("C1", "C2", "C3")
  )
  results <- adjust_uc(
    df_observed,
    u_model_coefs = c(
      u_model$coef[1],
      u_model$coef[2],
      u_model$coef[3],
      u_model$coef[4],
      u_model$coef[5],
      u_model$coef[6]
    )
  )
  est[i] <- results$estimate
}

or_true <- summary(nobias_model)$coef[2, 1]
or_adjusted <- median(est)

test_that("adjust_uc, 3 confounders: OR and CI output", {
  expect_gt(or_adjusted, or_true - 0.1)
  expect_lt(or_adjusted, or_true + 0.1)
  expect_vector(
    single_run$ci,
    ptype = double(),
    size = 2
  )
})
