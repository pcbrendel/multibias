set.seed(1234)
n <- 10000
nreps <- 10

# cont X just for testing that function runs
df_om$X_cont <- plogis(df_om$X) + rnorm(nrow(df_om), mean = 0, sd = 0.1)

# 0 confounders

nobias_model <- glm(Y ~ X,
                    family = binomial(link = "logit"),
                    data = df_om_source)

y_model <- glm(Y ~ X + Ystar,
               family = binomial(link = "logit"),
               data = df_om_source)

single_run <- adjust_om(
  df_om,
  exposure = "X_cont",
  outcome = "Ystar",
  y_model_coefs = c(
    y_model$coef[1],
    y_model$coef[2],
    y_model$coef[3]
  )
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_om[sample(seq_len(n), n, replace = TRUE), ]
  results <- adjust_om(
    bdf,
    exposure = "X",
    outcome = "Ystar",
    y_model_coefs = c(
      y_model$coef[1],
      y_model$coef[2],
      y_model$coef[3]
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

nobias_model <- glm(Y ~ X + C1,
                    family = binomial(link = "logit"),
                    data = df_om_source)

y_model <- glm(Y ~ X + Ystar + C1,
               family = binomial(link = "logit"),
               data = df_om_source)

single_run <- adjust_om(
  df_om,
  exposure = "X_cont",
  outcome = "Ystar",
  confounders = "C1",
  y_model_coefs = c(
    y_model$coef[1],
    y_model$coef[2],
    y_model$coef[3],
    y_model$coef[4]
  )
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_om[sample(seq_len(n), n, replace = TRUE), ]
  results <- adjust_om(
    bdf,
    exposure = "X",
    outcome = "Ystar",
    confounders = "C1",
    y_model_coefs = c(
      y_model$coef[1],
      y_model$coef[2],
      y_model$coef[3],
      y_model$coef[4]
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

nobias_model <- glm(Y ~ X + C1 + C2,
                    family = binomial(link = "logit"),
                    data = df_om_source)

y_model <- glm(Y ~ X + Ystar + C1 + C2,
               family = binomial(link = "logit"),
               data = df_om_source)

single_run <- adjust_om(
  df_om,
  exposure = "X_cont",
  outcome = "Ystar",
  confounders = c("C1", "C2"),
  y_model_coefs = c(
    y_model$coef[1],
    y_model$coef[2],
    y_model$coef[3],
    y_model$coef[4],
    y_model$coef[5]
  )
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_om[sample(seq_len(n), n, replace = TRUE), ]
  results <- adjust_om(
    bdf,
    exposure = "X",
    outcome = "Ystar",
    confounders = c("C1", "C2"),
    y_model_coefs = c(
      y_model$coef[1],
      y_model$coef[2],
      y_model$coef[3],
      y_model$coef[4],
      y_model$coef[5]
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

nobias_model <- glm(Y ~ X + C1 + C2 + C3,
                    family = binomial(link = "logit"),
                    data = df_om_source)

y_model <- glm(Y ~ X + Ystar + C1 + C2 + C3,
               family = binomial(link = "logit"),
               data = df_om_source)

single_run <- adjust_om(
  df_om,
  exposure = "X_cont",
  outcome = "Ystar",
  confounders = c("C1", "C2", "C3"),
  y_model_coefs = c(
    y_model$coef[1],
    y_model$coef[2],
    y_model$coef[3],
    y_model$coef[4],
    y_model$coef[5],
    y_model$coef[6]
  )
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_om[sample(seq_len(n), n, replace = TRUE), ]
  results <- adjust_om(
    bdf,
    exposure = "X",
    outcome = "Ystar",
    confounders = c("C1", "C2", "C3"),
    y_model_coefs = c(
      y_model$coef[1],
      y_model$coef[2],
      y_model$coef[3],
      y_model$coef[4],
      y_model$coef[5],
      y_model$coef[6]
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