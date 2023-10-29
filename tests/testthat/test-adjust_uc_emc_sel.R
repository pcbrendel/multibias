set.seed(1234)
# 0 confounders

nobias_model <- glm(Y ~ X + U,
                    family = binomial(link = "logit"),
                    data = df_uc_emc_sel_source)

u_model <- glm(U ~ X + Y,
               family = binomial(link = "logit"),
               data = df_uc_emc_sel_source)
x_model <- glm(X ~ Xstar + Y,
               family = binomial(link = "logit"),
               data = df_uc_emc_sel_source)
s_model <- glm(S ~ Xstar + Y,
               family = binomial(link = "logit"),
               data = df_uc_emc_sel_source)

single_run <- adjust_uc_emc_sel(
  df_uc_emc_sel,
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
  ),
  s_model_coefs = c(
    s_model$coef[1],
    s_model$coef[2],
    s_model$coef[3]
  )
)

n <- 100000
nreps <- 10
est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_emc_sel[sample(seq_len(n), n, replace = TRUE), ]
  results <- adjust_uc_emc_sel(
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

nobias_model <- glm(Y ~ X + C1 + U,
                    family = binomial(link = "logit"),
                    data = df_uc_emc_sel_source)

u_model <- glm(U ~ X + Y,
               family = binomial(link = "logit"),
               data = df_uc_emc_sel_source)
x_model <- glm(X ~ Xstar + Y + C1,
               family = binomial(link = "logit"),
               data = df_uc_emc_sel_source)
s_model <- glm(S ~ Xstar + Y + C1,
               family = binomial(link = "logit"),
               data = df_uc_emc_sel_source)

single_run <- adjust_uc_emc_sel(
  df_uc_emc_sel,
  "Xstar",
  "Y",
  "C1",
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
  ),
  s_model_coefs = c(
    s_model$coef[1],
    s_model$coef[2],
    s_model$coef[3],
    s_model$coef[4]
  )
)

n <- 100000
nreps <- 10
est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_emc_sel[sample(seq_len(n), n, replace = TRUE), ]
  results <- adjust_uc_emc_sel(
    bdf,
    "Xstar",
    "Y",
    "C1",
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
    ),
    s_model_coefs = c(
      s_model$coef[1],
      s_model$coef[2],
      s_model$coef[3],
      s_model$coef[4]
    )
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

nobias_model <- glm(Y ~ X + C1 + C2 + U,
                    family = binomial(link = "logit"),
                    data = df_uc_emc_sel_source)

u_model <- glm(U ~ X + Y,
               family = binomial(link = "logit"),
               data = df_uc_emc_sel_source)
x_model <- glm(X ~ Xstar + Y + C1 + C2,
               family = binomial(link = "logit"),
               data = df_uc_emc_sel_source)
s_model <- glm(S ~ Xstar + Y + C1 + C2,
               family = binomial(link = "logit"),
               data = df_uc_emc_sel_source)

single_run <- adjust_uc_emc_sel(
  df_uc_emc_sel,
  "Xstar",
  "Y",
  c("C1", "C2"),
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
  ),
  s_model_coefs = c(
    s_model$coef[1],
    s_model$coef[2],
    s_model$coef[3],
    s_model$coef[4],
    s_model$coef[5]
  )
)

n <- 100000
nreps <- 10
est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_emc_sel[sample(seq_len(n), n, replace = TRUE), ]
  results <- adjust_uc_emc_sel(
    bdf,
    "Xstar",
    "Y",
    c("C1", "C2"),
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
    ),
    s_model_coefs = c(
      s_model$coef[1],
      s_model$coef[2],
      s_model$coef[3],
      s_model$coef[4],
      s_model$coef[5]
    )
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

nobias_model <- glm(Y ~ X + C1 + C2 + C3 + U,
                    family = binomial(link = "logit"),
                    data = df_uc_emc_sel_source)

u_model <- glm(U ~ X + Y,
               family = binomial(link = "logit"),
               data = df_uc_emc_sel_source)
x_model <- glm(X ~ Xstar + Y + C1 + C2 + C3,
               family = binomial(link = "logit"),
               data = df_uc_emc_sel_source)
s_model <- glm(S ~ Xstar + Y + C1 + C2 + C3,
               family = binomial(link = "logit"),
               data = df_uc_emc_sel_source)

single_run <- adjust_uc_emc_sel(
  df_uc_emc_sel,
  "Xstar",
  "Y",
  c("C1", "C2", "C3"),
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
  ),
  s_model_coefs = c(
    s_model$coef[1],
    s_model$coef[2],
    s_model$coef[3],
    s_model$coef[4],
    s_model$coef[5],
    s_model$coef[6]
  )
)

n <- 100000
nreps <- 10
est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_emc_sel[sample(seq_len(n), n, replace = TRUE), ]
  results <- adjust_uc_emc_sel(
    bdf,
    "Xstar",
    "Y",
    c("C1", "C2", "C3"),
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
    ),
    s_model_coefs = c(
      s_model$coef[1],
      s_model$coef[2],
      s_model$coef[3],
      s_model$coef[4],
      s_model$coef[5],
      s_model$coef[6]
    )
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