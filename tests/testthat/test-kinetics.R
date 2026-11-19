test_that("Mark-Houwink (ISO 5351 constants) on a hand-computed value", {
  # (0.75 * 1000)^(1/0.905) = 1502.7
  expect_equal(dpv_from_viscosity(1000), 1502.7, tolerance = 1e-3)
  expect_equal(dpv_from_viscosity(0), 0)
  expect_error(dpv_from_viscosity(-1))
})

test_that("chain scissions: halving DP gives one scission per chain", {
  expect_equal(chain_scissions(500, 1000), 1)
})

test_that("Ekenstam recovers k exactly on noise-free data", {
  days <- rep(c(3, 7, 14, 28), each = 3)
  dp <- 1 / (1 / 1000 + 5e-5 * days)
  # exact fit: lm warns about a perfect fit, expected here
  fit <- suppressWarnings(ekenstam_fit(days, dp, 1000))
  expect_equal(fit$k, 5e-5, tolerance = 1e-8)
})

test_that("Ekenstam recovers k on noisy data and the CI covers it", {
  set.seed(1)
  days <- rep(c(3, 7, 14, 28), each = 5)
  dp <- 1 / (1 / 1000 + 5e-5 * days) * rnorm(length(days), 1, 0.02)
  fit <- ekenstam_fit(days, dp, 1000)
  expect_equal(fit$k, 5e-5, tolerance = 0.1)
  expect_true(fit$k_low < 5e-5 && 5e-5 < fit$k_high)
})

test_that("Arrhenius recovers Ea", {
  Ea <- 100e3
  temp <- c(70, 80, 90)
  k <- 1e-5 * exp(-Ea / 8.314462618 * (1 / (temp + 273.15) - 1 / 353.15))
  # exact fit on 3 points: lm warns about a perfect fit, expected here
  expect_equal(suppressWarnings(arrhenius_fit(k, temp))$Ea_kJ_mol, 100, tolerance = 1e-6)
  expect_error(arrhenius_fit(k[1:2], temp[1:2]), "3 distinct temperatures")
})
