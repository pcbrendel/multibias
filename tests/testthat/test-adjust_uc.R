results <- adjust_uc(
  evans,
  exposure = "SMK",
  outcome = "CHD",
  confounders = "AGE",
  u_model_coefs = c(qlogis(0.25), log(0.5), log(2.5), log(2)),
)

test_that("odds ratio and confidence interval output", {
  expect_vector(
    results$ci,
    ptype = double(),
    size = 2
  )
})
