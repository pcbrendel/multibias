nobias_model <- glm(Y ~ X + C1 + C2 + C3 + U,
                    family = binomial(link = "logit"),
                    data = df_uc_emc_sel_source)

# library(nnet)
# xu_model <- multinom(
#   paste(X, U) ~ Xstar + Y + C1 + C2 + C3,
#   data = df_uc_emc_sel_source
# )

s_model <- glm(S ~ Xstar + Y + C1 + C2 + C3,
               family = binomial(link = "logit"),
               data = df_uc_emc_sel_source)

single_run <- adjust_multinom_uc_emc_sel(
  df_uc_emc_sel,
  exposure = "Xstar",
  outcome = "Y",
  confounders = c("C1", "C2", "C3"),
  x1u0_model_coefs = c(-2.78, 1.62, 0.61, 0.36, -0.27, 0.88),
  x0u1_model_coefs = c(-0.17, -0.01, 0.71, -0.08, 0.07, -0.15),
  x1u1_model_coefs = c(-2.36, 1.62, 1.29, 0.25, -0.06, 0.74),
  s_model_coefs = c(0.00, 0.26, 0.78, 0.03, -0.02, 0.10)
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
    x1u0_model_coefs = c(-2.78, 1.62, 0.61, 0.36, -0.27, 0.88),
    x0u1_model_coefs = c(-0.17, -0.01, 0.71, -0.08, 0.07, -0.15),
    x1u1_model_coefs = c(-2.36, 1.62, 1.29, 0.25, -0.06, 0.74),
    s_model_coefs = c(0.00, 0.26, 0.78, 0.03, -0.02, 0.10)
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
