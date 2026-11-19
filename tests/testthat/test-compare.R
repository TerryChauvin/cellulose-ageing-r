test_that("before_after detects a clear decrease, paired", {
  df <- data.frame(materiau = "A", temperature_C = 80,
                   jours = rep(c(0, 28), each = 4),
                   dpv = c(1000, 1010, 990, 1005, 600, 615, 585, 604))
  r <- before_after(df, "dpv", 80)
  expect_equal(r$test, "paired t")
  expect_lt(r$diff, -300)
  expect_lt(r$p_value, 0.001)
})
