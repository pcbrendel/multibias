nobias_model <- glm(Y ~ X + C1 + U,
                    family = binomial(link = "logit"),
                    data = df_uc_omc_source)

# library(nnet)
# uy_model <- multinom(
#   paste0(U, Y) ~ X + Ystar + C1,
#   data = df_uc_omc_source
# )
# summary(uy_model)

single_run <- adjust_multinom_uc_omc(
  df_uc_omc,
  exposure = "X",
  outcome = "Ystar",
  confounders = "C1",
  u1y0_model_coefs = c(-0.19, 0.61, 0.00, -0.07),
  u0y1_model_coefs = c(-3.21, 0.60, 1.60, 0.36),
  u1y1_model_coefs = c(-2.72, 1.24, 1.59, 0.34)
)

n <- 100000
nreps <- 10
est <- vector()
for (i in 1:nreps) {
  bdf <- df_uc_omc[sample(seq_len(n), n, replace = TRUE), ]
  results <- adjust_multinom_uc_omc(
    bdf,
    exposure = "X",
    outcome = "Ystar",
    confounders = "C1",
    u1y0_model_coefs = c(-0.19, 0.61, 0.00, -0.07),
    u0y1_model_coefs = c(-3.21, 0.60, 1.60, 0.36),
    u1y1_model_coefs = c(-2.72, 1.24, 1.59, 0.34)
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