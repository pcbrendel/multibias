nobias_model <- glm(Y ~ X + C1,
                    family = binomial(link = "logit"),
                    data = df_sel_source)

s_model <- glm(S ~ X + Y,
               family = binomial(link = "logit"),
               data = df_sel_source)

single_run <- adjust_sel(
  df_sel,
  "X",
  "Y",
  "C1",
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
  bdf <- df_sel[sample(seq_len(n), n, replace = TRUE), ]
  results <- adjust_sel(
    bdf,
    "X",
    "Y",
    "C1",
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

test_that("odds ratio and confidence interval output", {
  expect_gt(or_adjusted, or_true - 0.15)
  expect_lt(or_adjusted, or_true + 0.15)
  expect_vector(
    single_run$ci,
    ptype = double(),
    size = 2
  )
})