results <- adjust_emc(
  evans,
  exposure = "SMK",
  outcome = "CHD",
  confounders = "HPT",
  x_model_coefs = c(qlogis(0.01), log(6), log(2), log(2))
)

test_that("odds ratio and confidence interval output", {
  expect_vector(
    results$ci,
    ptype = double(),
    size = 2
  )
})
