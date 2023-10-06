nobias_model <- glm(Y ~ X + C1,
                    family = binomial(link = "logit"),
                    data = df_emc_source)

x_model <- glm(X ~ Xstar + Y + C1,
               family = binomial(link = "logit"),
               data = df_emc_source)

single_run <- adjust_emc(
  df_emc,
  "Xstar",
  "Y",
  "C1",
  x_model_coefs = c(
    x_model$coef[1],
    x_model$coef[2],
    x_model$coef[3],
    x_model$coef[4]
  )
)

n <- 100000
nreps <- 10
est <- vector()
for (i in 1:nreps) {
  bdf <- df_emc[sample(seq_len(n), n, replace = TRUE), ]
  results <- adjust_emc(
    bdf,
    exposure = "Xstar",
    outcome = "Y",
    confounders = "C1",
    x_model_coefs = c(
      x_model$coef[1],
      x_model$coef[2],
      x_model$coef[3],
      x_model$coef[4]
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