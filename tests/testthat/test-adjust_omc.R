nobias_model <- glm(Y ~ X + C1,
                    family = binomial(link = "logit"),
                    data = df_omc_source)

y_model <- glm(Y ~ X + Ystar + C1,
               family = binomial(link = "logit"),
               data = df_omc_source)

single_run <- adjust_omc(
  df_omc,
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

n <- 100000
nreps <- 10
est <- vector()
for (i in 1:nreps) {
  bdf <- df_omc[sample(seq_len(n), n, replace = TRUE), ]
  results <- adjust_omc(
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
  expect_gt(or_adjusted, or_true - 0.15)
  expect_lt(or_adjusted, or_true + 0.15)
  expect_vector(
    single_run$ci,
    ptype = double(),
    size = 2
  )
})