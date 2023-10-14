nobias_model <- glm(Y ~ X + C1 + U,
                    family = binomial(link = "logit"),
                    data = df_uc_emc_source)

# library(nnet)
# xu_model <- multinom(
#   paste0(X, U) ~ Xstar + Y + C1,
#   data = df_uc_emc_source
# )
# summary(xu_model)

single_run <- adjust_multinom_uc_emc(
  df_uc_emc,
  exposure = "Xstar",
  outcome = "Y",
  confounders = "C1",
  x1u0_model_coefs = c(-2.82, 1.62, 0.68, -0.06),
  x0u1_model_coefs = c(-0.20, 0.00, 0.68, -0.05),
  x1u1_model_coefs = c(-2.36, 1.62, 1.29, 0.27)
)

n <- 100000
nreps <- 10
est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_emc[sample(seq_len(n), n, replace = TRUE), ]
  results <- adjust_multinom_uc_emc(
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
  expect_gt(or_adjusted, or_true - 0.15)
  expect_lt(or_adjusted, or_true + 0.15)
  expect_vector(
    single_run$ci,
    ptype = double(),
    size = 2
  )
})