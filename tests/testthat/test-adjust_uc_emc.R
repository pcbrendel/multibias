set.seed(1234)
n <- 20000
nreps <- 10

# SEPARATE BIAS PARAMETERS

# cont Y just for testing that function runs
df_uc_emc$Y_cont <- plogis(df_uc_emc$Y) +
  rnorm(nrow(df_uc_emc), mean = 0, sd = 0.1)

# 0 confounders

nobias_model <- glm(Y ~ X + U,
                    family = binomial(link = "logit"),
                    data = df_uc_emc_source)

u_model <- glm(U ~ X + Y,
               family = binomial(link = "logit"),
               data = df_uc_emc_source)
x_model <- glm(X ~ Xstar + Y,
               family = binomial(link = "logit"),
               data = df_uc_emc_source)

single_run <- adjust_uc_emc(
  df_uc_emc,
  exposure = "Xstar",
  outcome = "Y_cont",
  u_model_coefs = c(
    u_model$coef[1],
    u_model$coef[2],
    u_model$coef[3]
  ),
  x_model_coefs = c(
    x_model$coef[1],
    x_model$coef[2],
    x_model$coef[3]
  )
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_emc[sample(seq_len(n), n, replace = TRUE), ]
  results <- adjust_uc_emc(
    bdf,
    exposure = "Xstar",
    outcome = "Y",
    u_model_coefs = c(
      u_model$coef[1],
      u_model$coef[2],
      u_model$coef[3]
    ),
    x_model_coefs = c(
      x_model$coef[1],
      x_model$coef[2],
      x_model$coef[3]
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

nobias_model <- glm(Y ~ X + C1 + U,
                    family = binomial(link = "logit"),
                    data = df_uc_emc_source)

u_model <- glm(U ~ X + Y,
               family = binomial(link = "logit"),
               data = df_uc_emc_source)
x_model <- glm(X ~ Xstar + Y + C1,
               family = binomial(link = "logit"),
               data = df_uc_emc_source)

single_run <- adjust_uc_emc(
  df_uc_emc,
  exposure = "Xstar",
  outcome = "Y_cont",
  confounders = "C1",
  u_model_coefs = c(
    u_model$coef[1],
    u_model$coef[2],
    u_model$coef[3]
  ),
  x_model_coefs = c(
    x_model$coef[1],
    x_model$coef[2],
    x_model$coef[3],
    x_model$coef[4]
  )
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_emc[sample(seq_len(n), n, replace = TRUE), ]
  results <- adjust_uc_emc(
    bdf,
    exposure = "Xstar",
    outcome = "Y",
    confounders = "C1",
    u_model_coefs = c(
      u_model$coef[1],
      u_model$coef[2],
      u_model$coef[3]
    ),
    x_model_coefs = c(
      x_model$coef[1],
      x_model$coef[2],
      x_model$coef[3],
      x_model$coef[4]
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

nobias_model <- glm(Y ~ X + C1 + C2 + U,
                    family = binomial(link = "logit"),
                    data = df_uc_emc_source)

u_model <- glm(U ~ X + Y,
               family = binomial(link = "logit"),
               data = df_uc_emc_source)
x_model <- glm(X ~ Xstar + Y + C1 + C2,
               family = binomial(link = "logit"),
               data = df_uc_emc_source)

single_run <- adjust_uc_emc(
  df_uc_emc,
  exposure = "Xstar",
  outcome = "Y_cont",
  confounders = c("C1", "C2"),
  u_model_coefs = c(
    u_model$coef[1],
    u_model$coef[2],
    u_model$coef[3]
  ),
  x_model_coefs = c(
    x_model$coef[1],
    x_model$coef[2],
    x_model$coef[3],
    x_model$coef[4],
    x_model$coef[5]
  )
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_emc[sample(seq_len(n), n, replace = TRUE), ]
  results <- adjust_uc_emc(
    bdf,
    exposure = "Xstar",
    outcome = "Y",
    confounders = c("C1", "C2"),
    u_model_coefs = c(
      u_model$coef[1],
      u_model$coef[2],
      u_model$coef[3]
    ),
    x_model_coefs = c(
      x_model$coef[1],
      x_model$coef[2],
      x_model$coef[3],
      x_model$coef[4],
      x_model$coef[5]
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

nobias_model <- glm(Y ~ X + C1 + C2 + C3 + U,
                    family = binomial(link = "logit"),
                    data = df_uc_emc_source)

u_model <- glm(U ~ X + Y,
               family = binomial(link = "logit"),
               data = df_uc_emc_source)
x_model <- glm(X ~ Xstar + Y + C1 + C2 + C3,
               family = binomial(link = "logit"),
               data = df_uc_emc_source)

single_run <- adjust_uc_emc(
  df_uc_emc,
  exposure = "Xstar",
  outcome = "Y_cont",
  confounders = c("C1", "C2", "C3"),
  u_model_coefs = c(
    u_model$coef[1],
    u_model$coef[2],
    u_model$coef[3]
  ),
  x_model_coefs = c(
    x_model$coef[1],
    x_model$coef[2],
    x_model$coef[3],
    x_model$coef[4],
    x_model$coef[5],
    x_model$coef[6]
  )
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_emc[sample(seq_len(n), n, replace = TRUE), ]
  results <- adjust_uc_emc(
    bdf,
    exposure = "Xstar",
    outcome = "Y",
    confounders = c("C1", "C2", "C3"),
    u_model_coefs = c(
      u_model$coef[1],
      u_model$coef[2],
      u_model$coef[3]
    ),
    x_model_coefs = c(
      x_model$coef[1],
      x_model$coef[2],
      x_model$coef[3],
      x_model$coef[4],
      x_model$coef[5],
      x_model$coef[6]
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

# MULTINOMIAL BIAS PARAMETERS

# library(nnet)

set.seed(1234)
n <- 20000
nreps <- 10

# 0 confounders

nobias_model <- glm(Y ~ X + U,
                    family = binomial(link = "logit"),
                    data = df_uc_emc_source)

# xu_model <- multinom(
#   paste0(X, U) ~ Xstar + Y,
#   data = df_uc_emc_source
# )
# summary(xu_model)

single_run <- adjust_uc_emc(
  df_uc_emc,
  exposure = "Xstar",
  outcome = "Y",
  x1u0_model_coefs = c(-1.91, 1.63, 0.71),
  x0u1_model_coefs = c(-0.23, 0.02, 0.01),
  x1u1_model_coefs = c(-1.47, 1.63, 0.71)
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_emc[sample(seq_len(n), n, replace = TRUE), ]
  results <- adjust_uc_emc(
    bdf,
    exposure = "Xstar",
    outcome = "Y",
    x1u0_model_coefs = c(-1.91, 1.63, 0.71),
    x0u1_model_coefs = c(-0.23, 0.02, 0.01),
    x1u1_model_coefs = c(-1.47, 1.63, 0.71)
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

nobias_model <- glm(Y ~ X + C1 + U,
                    family = binomial(link = "logit"),
                    data = df_uc_emc_source)

# xu_model <- multinom(
#   paste0(X, U) ~ Xstar + Y + C1,
#   data = df_uc_emc_source
# )
# summary(xu_model)

single_run <- adjust_uc_emc(
  df_uc_emc,
  exposure = "Xstar",
  outcome = "Y",
  confounders = "C1",
  x1u0_model_coefs = c(-2.82, 1.62, 0.68, -0.06),
  x0u1_model_coefs = c(-0.20, 0.00, 0.68, -0.05),
  x1u1_model_coefs = c(-2.36, 1.62, 1.29, 0.27)
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_emc[sample(seq_len(n), n, replace = TRUE), ]
  results <- adjust_uc_emc(
    bdf,
    exposure = "Xstar",
    outcome = "Y",
    confounders = "C1",
    x1u0_model_coefs = c(-2.82, 1.62, 0.68, -0.06),
    x0u1_model_coefs = c(-0.20, 0.00, 0.68, -0.05),
    x1u1_model_coefs = c(-2.36, 1.62, 1.29, 0.27)
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

nobias_model <- glm(Y ~ X + C1 + C2 + U,
                    family = binomial(link = "logit"),
                    data = df_uc_emc_source)

# xu_model <- multinom(
#   paste0(X, U) ~ Xstar + Y + C1 + C2,
#   data = df_uc_emc_source
# )
# summary(xu_model)

single_run <- adjust_uc_emc(
  df_uc_emc,
  exposure = "Xstar",
  outcome = "Y",
  confounders = c("C1", "C2"),
  x1u0_model_coefs = c(-2.05, 1.63, 0.66, 0.34, -0.20),
  x0u1_model_coefs = c(-0.22, 0.02, 0.02, -0.04, 0.02),
  x1u1_model_coefs = c(-1.58, 1.63, 0.67, 0.28, -0.14)
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_emc[sample(seq_len(n), n, replace = TRUE), ]
  results <- adjust_uc_emc(
    bdf,
    exposure = "Xstar",
    outcome = "Y",
    confounders = c("C1", "C2"),
    x1u0_model_coefs = c(-2.05, 1.63, 0.66, 0.34, -0.20),
    x0u1_model_coefs = c(-0.22, 0.02, 0.02, -0.04, 0.02),
    x1u1_model_coefs = c(-1.58, 1.63, 0.67, 0.28, -0.14)
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

nobias_model <- glm(Y ~ X + C1 + C2 + C3 + U,
                    family = binomial(link = "logit"),
                    data = df_uc_emc_source)

# xu_model <- multinom(
#   paste0(X, U) ~ Xstar + Y + C1 + C2 + C3,
#   data = df_uc_emc_source
# )
# summary(xu_model)

single_run <- adjust_uc_emc(
  df_uc_emc,
  exposure = "Xstar",
  outcome = "Y",
  confounders = c("C1", "C2", "C3"),
  x1u0_model_coefs = c(-2.79, 1.63, 0.63, 0.35, -0.21, 0.89),
  x0u1_model_coefs = c(-0.12, 0.02, 0.02, -0.04, 0.02, -0.12),
  x1u1_model_coefs = c(-2.20, 1.63, 0.64, 0.29, -0.15, 0.74)
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_emc[sample(seq_len(n), n, replace = TRUE), ]
  results <- adjust_uc_emc(
    bdf,
    exposure = "Xstar",
    outcome = "Y",
    confounders = c("C1", "C2", "C3"),
    x1u0_model_coefs = c(-2.79, 1.63, 0.63, 0.35, -0.21, 0.89),
    x0u1_model_coefs = c(-0.12, 0.02, 0.02, -0.04, 0.02, -0.12),
    x1u1_model_coefs = c(-2.20, 1.63, 0.64, 0.29, -0.15, 0.74)
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
