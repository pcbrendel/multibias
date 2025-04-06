set.seed(1234)
n <- 10000
nreps <- 10

# SEPARATE BIAS PARAMETERS

# cont X just for testing that function runs
df_uc_om$X_cont <- plogis(df_uc_om$X) +
  rnorm(nrow(df_uc_om), mean = 0, sd = 0.1)

# 0 confounders

nobias_model <- glm(Y ~ X + U,
                    family = binomial(link = "logit"),
                    data = df_uc_om_source)

u_model <- glm(U ~ X + Y,
               family = binomial(link = "logit"),
               data = df_uc_om_source)
y_model <- glm(Y ~ X + Ystar,
               family = binomial(link = "logit"),
               data = df_uc_om_source)

df_observed <- data_observed(
  df_uc_om,
  bias = c("uc", "om"),
  exposure = "X_cont",
  outcome = "Ystar",
  confounders = NULL
)
list_for_uc_om <- list(
  u = as.vector(coef(u_model)),
  y = as.vector(coef(y_model))
)
bp_uc_om <- bias_params(coef_list = list_for_uc_om)

single_run <- adjust_uc_om(
  df_observed,
  bias_params = bp_uc_om
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_om[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    bias = c("uc", "om"),
    exposure = "X",
    outcome = "Ystar",
    confounders = NULL
  )
  results <- adjust_uc_om(
    df_observed,
    bias_params = bp_uc_om
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
                    data = df_uc_om_source)

u_model <- glm(U ~ X + Y,
               family = binomial(link = "logit"),
               data = df_uc_om_source)
y_model <- glm(Y ~ X + Ystar + C1 + C2 + C3,
               family = binomial(link = "logit"),
               data = df_uc_om_source)

df_observed <- data_observed(
  df_uc_om,
  bias = c("uc", "om"),
  exposure = "X_cont",
  outcome = "Ystar",
  confounders = c("C1", "C2", "C3")
)
list_for_uc_om <- list(
  u = as.vector(coef(u_model)),
  y = as.vector(coef(y_model))
)
bp_uc_om <- bias_params(coef_list = list_for_uc_om)

single_run <- adjust_uc_om(
  df_observed,
  bias_params = bp_uc_om
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_om[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    bias = c("uc", "om"),
    exposure = "X",
    outcome = "Ystar",
    confounders = c("C1", "C2", "C3")
  )
  results <- adjust_uc_om(
    df_observed,
    bias_params = bp_uc_om
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
n <- 10000
nreps <- 10

# 0 confounders

nobias_model <- glm(Y ~ X + U,
                    family = binomial(link = "logit"),
                    data = df_uc_om_source)

# uy_model <- multinom(
#   paste0(U, Y) ~ X + Ystar,
#   data = df_uc_om_source
# )
# summary(uy_model, digits = 2)

df_observed <- data_observed(
  df_uc_om,
  bias = c("uc", "om"),
  exposure = "X_cont",
  outcome = "Ystar",
  confounders = NULL
)
list_for_uc_om <- list(
  u1y0 = c(-0.32, 0.59, 0.01),
  u0y1 = c(-2.98, 0.71, 1.65),
  u1y1 = c(-2.59, 1.27, 1.64)

)
bp_uc_om <- bias_params(coef_list = list_for_uc_om)

single_run <- adjust_uc_om(
  df_observed,
  bias_params = bp_uc_om
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_om[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    bias = c("uc", "om"),
    exposure = "X",
    outcome = "Ystar",
    confounders = NULL
  )
  results <- adjust_uc_om(
    df_observed,
    bias_params = bp_uc_om
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
                    data = df_uc_om_source)

# uy_model <- multinom(
#   paste0(U, Y) ~ X + Ystar + C1 + C2 + C3,
#   data = df_uc_om_source
# )
# summary(uy_model, digits = 2)

df_observed <- data_observed(
  df_uc_om,
  bias = c("uc", "om"),
  exposure = "X_cont",
  outcome = "Ystar",
  confounders = c("C1", "C2", "C3")
)
list_for_uc_om <- list(
  u1y0 = c(-0.2, 0.62, 0.01, -0.08, 0.10, -0.15),
  u0y1 = c(-3.3, 0.63, 1.65, 0.42, -0.85, 0.26),
  u1y1 = c(-2.7, 1.22, 1.64, 0.32, -0.77, 0.09)

)
bp_uc_om <- bias_params(coef_list = list_for_uc_om)

single_run <- adjust_uc_om(
  df_observed,
  bias_params = bp_uc_om
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_om[sample(seq_len(n), n, replace = TRUE), ]
  df_observed <- data_observed(
    bdf,
    bias = c("uc", "om"),
    exposure = "X",
    outcome = "Ystar",
    confounders = c("C1", "C2", "C3")
  )
  results <- adjust_uc_om(
    df_observed,
    bias_params = bp_uc_om
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

or_val <- adjust_uc_om(
  data_observed = data_observed(
    df_uc_om,
    bias = c("uc", "om"),
    exposure = "X",
    outcome = "Ystar",
    confounders = c("C1", "C2", "C3")
  ),
  data_validation = data_validation(
    df_uc_om_source,
    true_exposure = "X",
    true_outcome = "Y",
    confounders = c("C1", "C2", "C3", "U"),
    misclassified_outcome = "Ystar"
  )
)

test_that("adjust_uc_om, validation data", {
  expect_gt(or_val$estimate, or_true - 0.1)
  expect_lt(or_val$estimate, or_true + 0.1)
})
