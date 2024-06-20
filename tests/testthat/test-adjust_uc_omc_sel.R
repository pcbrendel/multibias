set.seed(1234)
n <- 10000
nreps <- 10

# 0 confounders

nobias_model <- glm(Y ~ X + U,
                    family = binomial(link = "logit"),
                    data = df_uc_omc_sel_source)

u_model <- glm(U ~ X + Y,
               family = binomial(link = "logit"),
               data = df_uc_omc_sel_source)
y_model <- glm(Y ~ X + Ystar,
               family = binomial(link = "logit"),
               data = df_uc_omc_sel_source)
s_model <- glm(S ~ X + Ystar,
               family = binomial(link = "logit"),
               data = df_uc_omc_sel_source)

single_run <- adjust_uc_omc_sel(
  df_uc_omc_sel,
  exposure = "X",
  outcome = "Ystar",
  u_model_coefs = c(
    u_model$coef[1],
    u_model$coef[2],
    u_model$coef[3]
  ),
  y_model_coefs = c(
    y_model$coef[1],
    y_model$coef[2],
    y_model$coef[3]
  ),
  s_model_coefs = c(
    s_model$coef[1],
    s_model$coef[2],
    s_model$coef[3]
  )
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_omc_sel[sample(seq_len(n), n, replace = TRUE), ]
  results <- adjust_uc_omc_sel(
    bdf,
    exposure = "X",
    outcome = "Ystar",
    u_model_coefs = c(
      u_model$coef[1],
      u_model$coef[2],
      u_model$coef[3]
    ),
    y_model_coefs = c(
      y_model$coef[1],
      y_model$coef[2],
      y_model$coef[3]
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
                    data = df_uc_omc_sel_source)

u_model <- glm(U ~ X + Y,
               family = binomial(link = "logit"),
               data = df_uc_omc_sel_source)
y_model <- glm(Y ~ X + Ystar + C1,
               family = binomial(link = "logit"),
               data = df_uc_omc_sel_source)
s_model <- glm(S ~ X + Ystar + C1,
               family = binomial(link = "logit"),
               data = df_uc_omc_sel_source)

single_run <- adjust_uc_omc_sel(
  df_uc_omc_sel,
  "X",
  "Ystar",
  "C1",
  u_model_coefs = c(
    u_model$coef[1],
    u_model$coef[2],
    u_model$coef[3]
  ),
  y_model_coefs = c(
    y_model$coef[1],
    y_model$coef[2],
    y_model$coef[3],
    y_model$coef[4]
  ),
  s_model_coefs = c(
    s_model$coef[1],
    s_model$coef[2],
    s_model$coef[3],
    s_model$coef[4]
  )
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_omc_sel[sample(seq_len(n), n, replace = TRUE), ]
  results <- adjust_uc_omc_sel(
    bdf,
    "X",
    "Ystar",
    "C1",
    u_model_coefs = c(
      u_model$coef[1],
      u_model$coef[2],
      u_model$coef[3]
    ),
    y_model_coefs = c(
      y_model$coef[1],
      y_model$coef[2],
      y_model$coef[3],
      y_model$coef[4]
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
                    data = df_uc_omc_sel_source)

u_model <- glm(U ~ X + Y,
               family = binomial(link = "logit"),
               data = df_uc_omc_sel_source)
y_model <- glm(Y ~ X + Ystar + C1 + C2,
               family = binomial(link = "logit"),
               data = df_uc_omc_sel_source)
s_model <- glm(S ~ X + Ystar + C1 + C2,
               family = binomial(link = "logit"),
               data = df_uc_omc_sel_source)

single_run <- adjust_uc_omc_sel(
  df_uc_omc_sel,
  "X",
  "Ystar",
  c("C1", "C2"),
  u_model_coefs = c(
    u_model$coef[1],
    u_model$coef[2],
    u_model$coef[3]
  ),
  y_model_coefs = c(
    y_model$coef[1],
    y_model$coef[2],
    y_model$coef[3],
    y_model$coef[4],
    y_model$coef[5]
  ),
  s_model_coefs = c(
    s_model$coef[1],
    s_model$coef[2],
    s_model$coef[3],
    s_model$coef[4],
    s_model$coef[5]
  )
)

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_omc_sel[sample(seq_len(n), n, replace = TRUE), ]
  results <- adjust_uc_omc_sel(
    bdf,
    "X",
    "Ystar",
    c("C1", "C2"),
    u_model_coefs = c(
      u_model$coef[1],
      u_model$coef[2],
      u_model$coef[3]
    ),
    y_model_coefs = c(
      y_model$coef[1],
      y_model$coef[2],
      y_model$coef[3],
      y_model$coef[4],
      y_model$coef[5]
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
                    data = df_uc_omc_sel_source)

u_model <- glm(U ~ X + Y,
               family = binomial(link = "logit"),
               data = df_uc_omc_sel_source)
y_model <- glm(Y ~ X + Ystar + C1 + C2 + C3,
               family = binomial(link = "logit"),
               data = df_uc_omc_sel_source)
s_model <- glm(S ~ X + Ystar + C1 + C2 + C3,
               family = binomial(link = "logit"),
               data = df_uc_omc_sel_source)

single_run <- adjust_uc_omc_sel(
  df_uc_omc_sel,
  "X",
  "Ystar",
  c("C1", "C2", "C3"),
  u_model_coefs = c(
    u_model$coef[1],
    u_model$coef[2],
    u_model$coef[3]
  ),
  y_model_coefs = c(
    y_model$coef[1],
    y_model$coef[2],
    y_model$coef[3],
    y_model$coef[4],
    y_model$coef[5],
    y_model$coef[6]
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

est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_omc_sel[sample(seq_len(n), n, replace = TRUE), ]
  results <- adjust_uc_omc_sel(
    bdf,
    "X",
    "Ystar",
    c("C1", "C2", "C3"),
    u_model_coefs = c(
      u_model$coef[1],
      u_model$coef[2],
      u_model$coef[3]
    ),
    y_model_coefs = c(
      y_model$coef[1],
      y_model$coef[2],
      y_model$coef[3],
      y_model$coef[4],
      y_model$coef[5],
      y_model$coef[6]
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