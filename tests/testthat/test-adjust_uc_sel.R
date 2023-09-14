results <- adjust_uc_sel(
  df_uc_sel,
  exposure = "X",
  outcome = "Y",
  confounders = "C1",
  u_model_coefs = c(-0.36, 0.96, 1.11, -0.12),
  s_model_coefs = c(0.02, 0.89, 0.90)
)

test_that("odds ratio and confidence interval output", {
  expect_gt(results$estimate, 1.5)
  expect_lt(results$estimate, 2.5)
  expect_vector(
    results$ci,
    ptype = double(),
    size = 2
  )
})
