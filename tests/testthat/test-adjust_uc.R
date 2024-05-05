set.seed(1234)
n <- 20000
nreps <- 10

# 0 confounders

nobias_model <- glm(Y ~ X + U,
                    family = binomial(link = "logit"),
                    data = df_uc_source)

u_model <- glm(U ~ X + Y,
               family = binomial(link = "logit"),
               data = df_uc_source)

single_run <- adjust_uc(
  df_uc,
  exposure = "X",
  outcome = "Y",
  u_model_coefs = c(
    u_model$coef[1],
    u_model$coef[2],
    u_model$coef[3]
  )
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc[sample(seq_len(n), n, replace = TRUE), ]
  results <- adjust_uc(
    bdf,
    exposure = "X",
    outcome = "Y",
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
                    data = df_uc_source)

u_model <- glm(U ~ X + Y + C1,
               family = binomial(link = "logit"),
               data = df_uc_source)

single_run <- adjust_uc(
  df_uc,
  exposure = "X",
  outcome = "Y",
  confounders = "C1",
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
  results <- adjust_uc(
    bdf,
    exposure = "X",
    outcome = "Y",
    confounders = "C1",
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
                    data = df_uc_source)

u_model <- glm(U ~ X + Y + C1 + C2,
               family = binomial(link = "logit"),
               data = df_uc_source)

single_run <- adjust_uc(
  df_uc,
  exposure = "X",
  outcome = "Y",
  confounders = c("C1", "C2"),
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
  results <- adjust_uc(
    bdf,
    exposure = "X",
    outcome = "Y",
    confounders = c("C1", "C2"),
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
                    data = df_uc_source)

u_model <- glm(U ~ X + Y + C1 + C2 + C3,
               family = binomial(link = "logit"),
               data = df_uc_source)

single_run <- adjust_uc(
  df_uc,
  exposure = "X",
  outcome = "Y",
  confounders = c("C1", "C2", "C3"),
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
  results <- adjust_uc(
    bdf,
    exposure = "X",
    outcome = "Y",
    confounders = c("C1", "C2", "C3"),
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
