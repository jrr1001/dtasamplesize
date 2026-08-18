
scenarios_list <- list(
  optimistic = list(prior_se = c(18, 2), prior_sp = c(19, 1),
                     prior_prev = c(5, 15)),
  neutral = list(prior_se = c(17, 3), prior_sp = c(2, 2),
                  prior_prev = c(5, 15)),
  pessimistic = list(prior_se = c(14, 6), prior_sp = c(16, 4),
                      prior_prev = c(5, 15))
)

common_args <- list(delta_se = 0.07, delta_sp = 0.05, delta_auc = 0,
                    N_range = seq(300, 1200, by = 100), B = 300, seed = 2026)

test_that("sensitivity_analysis returns the promised columns, one row per scenario, in order", {
  result <- do.call(sensitivity_analysis, c(list(scenarios_list), common_args))

  expect_true(is.data.frame(result))
  expect_identical(
    names(result),
    c("scenario", "N_total", "joint_assurance", "N_buderer", "n_diseased")
  )
  expect_equal(nrow(result), length(scenarios_list))
  expect_identical(result$scenario, names(scenarios_list))
})

test_that("sensitivity_analysis accepts a named list of scenarios", {
  result <- do.call(sensitivity_analysis, c(list(scenarios_list), common_args))
  expect_true(all(is.finite(result$N_total)))
  expect_true(all(result$joint_assurance >= 0 & result$joint_assurance <= 1))
  expect_true(all(is.finite(result$N_buderer)))
  expect_true(all(is.finite(result$n_diseased)))
})

test_that("sensitivity_analysis accepts a data frame of scenarios and agrees with the list form", {
  df_scenarios <- data.frame(
    scenario = c("optimistic", "pessimistic"),
    prior_se_shape1 = c(18, 14), prior_se_shape2 = c(2, 6),
    prior_sp_shape1 = c(19, 16), prior_sp_shape2 = c(1, 4),
    prior_prev_shape1 = c(5, 5), prior_prev_shape2 = c(15, 15),
    stringsAsFactors = FALSE
  )
  from_df <- do.call(sensitivity_analysis, c(list(df_scenarios), common_args))
  from_list <- do.call(
    sensitivity_analysis,
    c(list(scenarios_list[c("optimistic", "pessimistic")]), common_args)
  )

  expect_identical(names(from_df), names(from_list))
  expect_equal(nrow(from_df), 2)
  expect_identical(from_df$scenario, c("optimistic", "pessimistic"))
  # Same priors, same N_range/B/seed -> identical numeric results either way
  # the scenarios were supplied.
  expect_equal(from_df$N_total, from_list$N_total)
  expect_equal(from_df$joint_assurance, from_list$joint_assurance)
})

test_that("an informative error is raised when a scenario element is missing", {
  bad_list <- list(
    optimistic = list(prior_se = c(18, 2), prior_prev = c(5, 15))  # prior_sp missing
  )
  expect_error(
    do.call(sensitivity_analysis, c(list(bad_list), common_args)),
    "prior_sp"
  )
})

test_that("an informative error is raised when a data frame column is missing", {
  bad_df <- data.frame(
    scenario = "optimistic",
    prior_se_shape1 = 18, prior_se_shape2 = 2,
    prior_sp_shape1 = 19, prior_sp_shape2 = 1
    # prior_prev_shape1/2 missing
  )
  expect_error(
    do.call(sensitivity_analysis, c(list(bad_df), common_args)),
    "prior_prev"
  )
})

test_that("sensitivity_analysis rejects prior_se/prior_sp/prior_prev in ...", {
  expect_error(
    sensitivity_analysis(scenarios_list, prior_se = c(1, 1), B = 300),
    "scenarios"
  )
})

test_that("sensitivity_analysis requires named scenarios and unique names", {
  expect_error(
    sensitivity_analysis(list(scenarios_list[[1]]), B = 300),
    "named"
  )
  dup <- scenarios_list
  names(dup)[2] <- names(dup)[1]
  expect_error(
    do.call(sensitivity_analysis, c(list(dup), common_args)),
    "unique"
  )
})
