
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

# --- (v0.7.0) plot_assurance_curve() for dta_bam_curve objects ---------------

test_that("plot_assurance_curve returns a ggplot for a dta_bam_curve", {
  skip_if_not_installed("ggplot2")
  cv <- bam_assurance_curve(prior_se = c(17, 3), prior_sp = c(2, 2),
                            prior_prev = c(4, 16), N_range = seq(20, 120, by = 20))
  p <- plot_assurance_curve(cv)
  expect_s3_class(p, "ggplot")
  expect_null(p$labels$title)
  expect_match(p$labels$y, "exact", ignore.case = TRUE)
  expect_match(p$labels$y, "joint", ignore.case = TRUE)
  # reaching the target adds a vertical line at N_first
  # (delta_se = 0.4 is met by the Se prior alone with no data: the warning
  # about the effective analysis prior is expected)
  expect_warning(
    cv2 <- bam_assurance_curve(prior_se = c(17, 3), prior_sp = c(2, 2),
                               prior_prev = c(4, 16), delta_se = 0.4,
                               delta_sp = 0.4, N_range = 5:60),
    "analysis prior for Se"
  )
  expect_false(is.na(cv2$N_first))
  p2 <- plot_assurance_curve(cv2, title = "T")
  expect_s3_class(p2, "ggplot")
  expect_identical(p2$labels$title, "T")
  geoms <- vapply(p2$layers, function(l) class(l$geom)[1], "")
  expect_true("GeomVline" %in% geoms)
  expect_true("GeomHline" %in% geoms)
})

test_that("plot_assurance_curve names the convention in the y axis", {
  skip_if_not_installed("ggplot2")
  cv <- bam_assurance_curve(prior_se = c(25, 7), prior_sp = c(30, 2),
                            prior_prev = c(13.56, 122.06),
                            width_type = "lower_half_width", delta_se = 0.10,
                            delta_sp = 0.05, assurance_for = "se",
                            N_range = c(100L, 200L, 300L))
  p <- plot_assurance_curve(cv)
  expect_match(p$labels$y, "marginal", ignore.case = TRUE)
  expect_match(p$labels$y, "lower half-width", ignore.case = TRUE)
})

test_that("plot_assurance_curve refuses a bam_sample_size() object and says what to use", {
  skip_if_not_installed("ggplot2")
  r <- suppressWarnings(bam_sample_size(B = 100, n_range = 20:60,
                                        N_range = c(100L, 200L)))
  expect_error(plot_assurance_curve(r), "bam_assurance_curve")
})
