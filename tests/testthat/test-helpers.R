test_that("buderer_n reproduces known values", {
  # Se=0.85, d=0.07 -> n=100
  expect_equal(buderer_n(0.85, 0.07), 100)
  # Se=0.90, d=0.05 -> n=139
  expect_equal(buderer_n(0.90, 0.05), 139)
})

test_that("wilson_ci returns valid interval", {
  ci <- wilson_ci(85, 100)
  expect_true(ci["lower"] >= 0)
  expect_true(ci["upper"] <= 1)
  expect_true(ci["lower"] < ci["upper"])
  expect_true(ci["width"] > 0)
})

test_that("wald_ci returns valid interval", {
  ci <- wald_ci(85, 100)
  expect_true(ci["lower"] >= 0)
  expect_true(ci["upper"] <= 1)
  expect_true(ci["lower"] < ci["upper"])
  expect_true(ci["width"] > 0)
})

test_that("hanley_mcneil_var is positive", {
  v <- hanley_mcneil_var(0.75, 100, 100)
  expect_true(v > 0)
})
