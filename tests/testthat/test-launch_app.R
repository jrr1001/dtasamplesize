# launch_app() had zero coverage. We cannot start a Shiny server in a
# test, but we can guard the two things that actually break: the app files
# must ship, and app.R must be syntactically valid R (it was hand-edited in
# the v0.3.0 correction to add error handling and fix the impossible AUC
# default, so a syntax slip there would ship a broken app silently).

test_that("the Shiny app directory and app.R ship with the package", {
  app_dir <- system.file("shiny", package = "dtasamplesize")
  expect_true(nzchar(app_dir))
  expect_true(file.exists(file.path(app_dir, "app.R")))
})

test_that("inst/shiny/app.R parses without error", {
  app_file <- system.file("shiny", "app.R", package = "dtasamplesize")
  skip_if(app_file == "" || !file.exists(app_file))
  expect_no_error(parse(app_file))
})

test_that("launch_app errors clearly when shiny is unavailable", {
  # Only meaningful if shiny is genuinely not installed; otherwise skip.
  skip_if(requireNamespace("shiny", quietly = TRUE),
          "shiny is installed; cannot test the missing-dependency path")
  expect_error(launch_app(), "shiny")
})
