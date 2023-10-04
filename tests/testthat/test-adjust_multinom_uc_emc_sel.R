library(nnet)

nobias_model <- glm(Y ~ X + C1 + C2 + C3 + U,
                    family = binomial(link = "logit"),
                    data = df_uc_emc_sel_source)

xu_model <- multinom(
  paste(X, U) ~ Xstar + Y + C1 + C2 + C3,
  data = df_uc_emc_sel_source
)

s_model <- glm(S ~ Xstar + Y + C1 + C2 + C3,
               family = binomial(link = "logit"),
               data = df_uc_emc_sel_source)

single_run <- adjust_multinom_uc_emc_sel(
  df_uc_emc_sel,
  exposure = "Xstar",
  outcome = "Y",
  confounders = c("C1", "C2", "C3"),
  x1u0_model_coefs = c(
    summary(xu_model)$coefficients[2, 1],
    summary(xu_model)$coefficients[2, 2],
    summary(xu_model)$coefficients[2, 3],
    summary(xu_model)$coefficients[2, 4],
    summary(xu_model)$coefficients[2, 5],
    summary(xu_model)$coefficients[2, 6]
  ),
  x0u1_model_coefs = c(
    summary(xu_model)$coefficients[1, 1],
    summary(xu_model)$coefficients[1, 2],
    summary(xu_model)$coefficients[1, 3],
    summary(xu_model)$coefficients[1, 4],
    summary(xu_model)$coefficients[1, 5],
    summary(xu_model)$coefficients[1, 6]
  ),
  x1u1_model_coefs = c(
    summary(xu_model)$coefficients[3, 1],
    summary(xu_model)$coefficients[3, 2],
    summary(xu_model)$coefficients[3, 3],
    summary(xu_model)$coefficients[3, 4],
    summary(xu_model)$coefficients[3, 5],
    summary(xu_model)$coefficients[3, 6]
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
  results <- adjust_multinom_uc_emc_sel(
    bdf,
    exposure = "Xstar",
    outcome = "Y",
    confounders = c("C1", "C2", "C3"),
    x1u0_model_coefs = c(
      summary(xu_model)$coefficients[2, 1],
      summary(xu_model)$coefficients[2, 2],
      summary(xu_model)$coefficients[2, 3],
      summary(xu_model)$coefficients[2, 4],
      summary(xu_model)$coefficients[2, 5],
      summary(xu_model)$coefficients[2, 6]
    ),
    x0u1_model_coefs = c(
      summary(xu_model)$coefficients[1, 1],
      summary(xu_model)$coefficients[1, 2],
      summary(xu_model)$coefficients[1, 3],
      summary(xu_model)$coefficients[1, 4],
      summary(xu_model)$coefficients[1, 5],
      summary(xu_model)$coefficients[1, 6]
    ),
    x1u1_model_coefs = c(
      summary(xu_model)$coefficients[3, 1],
      summary(xu_model)$coefficients[3, 2],
      summary(xu_model)$coefficients[3, 3],
      summary(xu_model)$coefficients[3, 4],
      summary(xu_model)$coefficients[3, 5],
      summary(xu_model)$coefficients[3, 6]
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

test_that("odds ratio and confidence interval output", {
  expect_gt(or_adjusted, or_true - 0.15)
  expect_lt(or_adjusted, or_true + 0.15)
  expect_vector(
    single_run$ci,
    ptype = double(),
    size = 2
  )
})
