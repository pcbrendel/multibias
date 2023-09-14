results <- adjust_sel(
  evans,
  exposure = "SMK",
  outcome = "CHD",
  confounders = "HPT",
  s_model_coefs = c(qlogis(0.25), log(0.75), log(0.75))
)

test_that("odds ratio and confidence interval output", {
  expect_vector(
    results$ci,
    ptype = double(),
    size = 2
  )
})
