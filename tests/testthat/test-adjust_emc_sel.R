results <- adjust_emc_sel(
  df_emc_sel,
  exposure = "Xstar",
  outcome = "Y",
  confounders = "C1",
  x_model_coefs = c(-2.78, 1.62, 0.58, 0.34),
  s_model_coefs = c(0.04, 0.18, 0.92, 0.05)
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
