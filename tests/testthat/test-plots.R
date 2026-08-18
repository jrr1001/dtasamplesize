
test_that("plot_assurance_curve returns a ggplot object", {
  skip_if_not_installed("ggplot2")

  result <- ss_unified(B = 300, seed = 2026,
                       delta_se = 0.07, delta_sp = 0.05, delta_auc = 0,
                       N_range = seq(300, 1200, by = 100), full_grid = TRUE)
  p <- plot_assurance_curve(result)
  expect_s3_class(p, "ggplot")
})

test_that("plot_assurance_curve errors informatively without full_grid", {
  skip_if_not_installed("ggplot2")

  # A single-N search never reaches the target, hence the expected
  # convergence warning suppressed below; what is under test here is that
  # plot_assurance_curve() then refuses a 1-row grid with a clear message.
  result <- suppressWarnings(ss_unified(
    B = 300, seed = 2026,
    delta_se = 0.07, delta_sp = 0.05, delta_auc = 0,
    N_range = 300, full_grid = FALSE
  ))
  expect_error(plot_assurance_curve(result), "full_grid")
})

test_that("plot_inflation_heatmap returns a ggplot object", {
  skip_if_not_installed("ggplot2")

  result <- ss_imperfect_ref(B = 0, sensitivity_table = TRUE)
  p <- plot_inflation_heatmap(result)
  expect_s3_class(p, "ggplot")
})

test_that("plot_inflation_heatmap errors informatively without sensitivity_table", {
  skip_if_not_installed("ggplot2")

  result <- ss_imperfect_ref(B = 0, sensitivity_table = FALSE)
  expect_error(plot_inflation_heatmap(result), "sensitivity_table")
})

test_that("plot_method_comparison returns a ggplot object", {
  skip_if_not_installed("ggplot2")

  result <- ss_unified(B = 300, seed = 2026,
                       delta_se = 0.07, delta_sp = 0.05, delta_auc = 0,
                       N_range = seq(300, 1200, by = 100))
  p <- plot_method_comparison(result)
  expect_s3_class(p, "ggplot")
})

test_that("the plot_* functions require a dtasamplesize object", {
  skip_if_not_installed("ggplot2")

  expect_error(plot_assurance_curve(list()), "dtasamplesize")
  expect_error(plot_inflation_heatmap(list()), "dtasamplesize")
  expect_error(plot_method_comparison(list()), "dtasamplesize")
})
