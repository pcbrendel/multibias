# library(nnet)

set.seed(1234)
n <- 10000
nreps <- 10

# 0 confounders

nobias_model <- glm(Y ~ X,
                    family = binomial(link = "logit"),
                    data = df_emc_omc_source)

# xy_model <- multinom(
#   paste(X, Y) ~ Xstar + Ystar,
#   data = df_emc_omc_source
# )
# summary(xy_model)

single_run <- adjust_multinom_emc_omc(
  df_emc_omc,
  exposure = "Xstar",
  outcome = "Ystar",
  x1y0_model_coefs = c(-1.99, 1.62, 0.23),
  x0y1_model_coefs = c(-2.96, 0.22, 1.60),
  x1y1_model_coefs = c(-4.36, 1.82, 1.82)
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_emc_omc[sample(seq_len(n), n, replace = TRUE), ]
  results <- adjust_multinom_emc_omc(
    bdf,
    exposure = "Xstar",
    outcome = "Ystar",
    x1y0_model_coefs = c(-1.99, 1.62, 0.23),
    x0y1_model_coefs = c(-2.96, 0.22, 1.60),
    x1y1_model_coefs = c(-4.36, 1.82, 1.82)
  )
  est[i] <- results$estimate
}

or_true <- exp(summary(nobias_model)$coef[2, 1])
or_adjusted <- median(est)

test_that("0 confounders: odds ratio and confidence interval output", {
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
                    data = df_emc_omc_source)

# xy_model <- multinom(
#   paste(X, Y) ~ Xstar + Ystar + C1,
#   data = df_emc_omc_source
# )
# summary(xy_model)

single_run <- adjust_multinom_emc_omc(
  df_emc_omc,
  exposure = "Xstar",
  outcome = "Ystar",
  confounders = "C1",
  x1y0_model_coefs = c(-2.18, 1.63, 0.23, 0.36),
  x0y1_model_coefs = c(-3.17, 0.22, 1.60, 0.40),
  x1y1_model_coefs = c(-4.76, 1.82, 1.83, 0.72)
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_emc_omc[sample(seq_len(n), n, replace = TRUE), ]
  results <- adjust_multinom_emc_omc(
    bdf,
    exposure = "Xstar",
    outcome = "Ystar",
    confounders = "C1",
    x1y0_model_coefs = c(-2.18, 1.63, 0.23, 0.36),
    x0y1_model_coefs = c(-3.17, 0.22, 1.60, 0.40),
    x1y1_model_coefs = c(-4.76, 1.82, 1.83, 0.72)
  )
  est[i] <- results$estimate
}

or_true <- exp(summary(nobias_model)$coef[2, 1])
or_adjusted <- median(est)

test_that("1 confounder: odds ratio and confidence interval output", {
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
                    data = df_emc_omc_source)

# xy_model <- multinom(
#   paste(X, Y) ~ Xstar + Ystar + C1 + C2,
#   data = df_emc_omc_source
# )
# summary(xy_model)

single_run <- adjust_multinom_emc_omc(
  df_emc_omc,
  exposure = "Xstar",
  outcome = "Ystar",
  confounders = c("C1", "C2"),
  x1y0_model_coefs = c(-2.14, 1.63, 0.23, 0.36, -0.21),
  x0y1_model_coefs = c(-3.03, 0.22, 1.60, 0.41, -0.92),
  x1y1_model_coefs = c(-4.60, 1.83, 1.83, 0.72, -1.14)
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_emc_omc[sample(seq_len(n), n, replace = TRUE), ]
  results <- adjust_multinom_emc_omc(
    bdf,
    exposure = "Xstar",
    outcome = "Ystar",
    confounders = c("C1", "C2"),
    x1y0_model_coefs = c(-2.14, 1.63, 0.23, 0.36, -0.21),
    x0y1_model_coefs = c(-3.03, 0.22, 1.60, 0.41, -0.92),
    x1y1_model_coefs = c(-4.60, 1.83, 1.83, 0.72, -1.14)
  )
  est[i] <- results$estimate
}

or_true <- exp(summary(nobias_model)$coef[2, 1])
or_adjusted <- median(est)

test_that("2 confounders: odds ratio and confidence interval output", {
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
                    data = df_emc_omc_source)

# xy_model <- multinom(
#   paste(X, Y) ~ Xstar + Ystar + C1 + C2 + C3,
#   data = df_emc_omc_source
# )
# summary(xy_model)

single_run <- adjust_multinom_emc_omc(
  df_emc_omc,
  exposure = "Xstar",
  outcome = "Ystar",
  confounders = c("C1", "C2", "C3"),
  x1y0_model_coefs = c(-2.86, 1.63, 0.23, 0.37, -0.22, 0.87),
  x0y1_model_coefs = c(-3.26, 0.22, 1.60, 0.41, -0.93, 0.28),
  x1y1_model_coefs = c(-5.62, 1.83, 1.83, 0.74, -1.15, 1.19)
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_emc_omc[sample(seq_len(n), n, replace = TRUE), ]
  results <- adjust_multinom_emc_omc(
    bdf,
    exposure = "Xstar",
    outcome = "Ystar",
    confounders = c("C1", "C2", "C3"),
    x1y0_model_coefs = c(-2.86, 1.63, 0.23, 0.37, -0.22, 0.87),
    x0y1_model_coefs = c(-3.26, 0.22, 1.60, 0.41, -0.93, 0.28),
    x1y1_model_coefs = c(-5.62, 1.83, 1.83, 0.74, -1.15, 1.19)
  )
  est[i] <- results$estimate
}

or_true <- exp(summary(nobias_model)$coef[2, 1])
or_adjusted <- median(est)

test_that("3 confounders: odds ratio and confidence interval output", {
  expect_gt(or_adjusted, or_true - 0.1)
  expect_lt(or_adjusted, or_true + 0.1)
  expect_vector(
    single_run$ci,
    ptype = double(),
    size = 2
  )
})
