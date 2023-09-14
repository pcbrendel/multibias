results <- adjust_multinom_uc_emc(
  df_uc_emc,
  exposure = "Xstar",
  outcome = "Y",
  confounders = "C1",
  x1u0_model_coefs = c(-2.82, 1.62, 0.68, -0.06),
  x0u1_model_coefs = c(-0.20, 0.00, 0.68, -0.05),
  x1u1_model_coefs = c(-2.36, 1.62, 1.29, 0.27)
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