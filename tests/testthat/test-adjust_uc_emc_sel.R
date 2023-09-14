results <- adjust_uc_emc_sel(
  df_uc_emc_sel,
  exposure = "Xstar",
  outcome = "Y",
  confounders = c("C1", "C2", "C3"),
  u_model_coefs = c(-0.32, 0.59, 0.69),
  x_model_coefs = c(-2.44, 1.62, 0.72, 0.32, -0.15, 0.85),
  s_model_coefs = c(0.00, 0.26, 0.78, 0.03, -0.02, 0.10)
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
