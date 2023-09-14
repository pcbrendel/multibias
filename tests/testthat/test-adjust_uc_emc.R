results <- adjust_uc_emc(
  df_uc_emc,
  exposure = "Xstar",
  outcome = "Y",
  confounders = "C1",
  u_model_coefs = c(-0.23, 0.63, 0.66),
  x_model_coefs = c(-2.47, 1.62, 0.73, 0.32)
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
