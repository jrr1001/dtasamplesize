# Tests for bam_prior_from_pilot()

vap_args <- function(...) {
  base <- list(tp = 16, fn = 1, tn = 20, fp = 35,
               base_prior_se = c(9.9, 1.1),
               base_prior_sp = c(1, 1),
               base_prior_prev = c(12, 43))
  utils::modifyList(base, list(...))
}

test_that("VAP example of Wilson et al. 2022 is reproduced (a0 = 1, cohort)", {
  p <- do.call(bam_prior_from_pilot, vap_args(a0 = 1, design = "cohort"))
  expect_s3_class(p, "dta_bam_priors")
  expect_equal(p$prior_se,   c(25.9, 2.1), tolerance = 1e-12)
  expect_equal(p$prior_prev, c(29, 98),    tolerance = 1e-12)
  expect_equal(p$prior_sp,   c(21, 36),    tolerance = 1e-12)
  # use_in = "both" by default: analysis priors equal the updated ones
  expect_equal(p$analysis_prior_se, p$prior_se)
  expect_equal(p$analysis_prior_sp, p$prior_sp)
})

test_that("contract: exact component names and structure", {
  p <- do.call(bam_prior_from_pilot, vap_args())
  expect_identical(
    names(p),
    c("prior_se", "prior_sp", "prior_prev",
      "analysis_prior_se", "analysis_prior_sp",
      "base", "pilot", "a0", "design", "use_in"))
  expect_identical(names(p$base), c("prior_se", "prior_sp", "prior_prev"))
  expect_identical(names(p$pilot), c("tp", "fn", "tn", "fp"))
  expect_equal(p$base$prior_se, c(9.9, 1.1))
  expect_equal(p$pilot$tp, 16)
  expect_identical(p$design, "cohort")
  expect_identical(p$use_in, "both")
  expect_equal(p$a0, 1)
})

test_that("case_control leaves prevalence equal to the base prior", {
  p <- do.call(bam_prior_from_pilot, vap_args(design = "case_control"))
  expect_equal(p$prior_prev, c(12, 43))
  expect_equal(p$prior_se, c(25.9, 2.1), tolerance = 1e-12)
  expect_equal(p$prior_sp, c(21, 36), tolerance = 1e-12)
})

test_that("a0 = 0.5 halves the increments", {
  p <- do.call(bam_prior_from_pilot, vap_args(a0 = 0.5))
  expect_equal(p$prior_se,   c(9.9 + 8, 1.1 + 0.5), tolerance = 1e-12)
  expect_equal(p$prior_sp,   c(1 + 10, 1 + 17.5),  tolerance = 1e-12)
  expect_equal(p$prior_prev, c(12 + 8.5, 43 + 27.5), tolerance = 1e-12)
  expect_equal(p$a0, 0.5)
})

test_that("use_in = 'design' keeps the base prior as analysis prior", {
  p <- do.call(bam_prior_from_pilot, vap_args(use_in = "design"))
  expect_equal(p$prior_se, c(25.9, 2.1), tolerance = 1e-12)
  expect_equal(p$analysis_prior_se, c(9.9, 1.1))
  expect_equal(p$analysis_prior_sp, c(1, 1))
  expect_false(is.null(p$analysis_prior_se))
  expect_identical(p$use_in, "design")
})

test_that("input errors", {
  f <- function(...) do.call(bam_prior_from_pilot, vap_args(...))
  expect_error(f(tp = -1), "tp")
  expect_error(f(fn = 1.5), "fn")
  expect_error(f(tn = NA), "tn")
  expect_error(f(fp = c(1, 2)), "fp")
  expect_error(f(tp = "a"), "tp")
  expect_error(f(a0 = 0), "a0")
  expect_error(f(a0 = 1.01), "a0")
  expect_error(f(a0 = -0.2), "a0")
  expect_error(f(a0 = c(0.5, 1)), "a0")
  expect_error(f(base_prior_se = c(1, 0)), "base_prior_se")
  expect_error(f(base_prior_sp = c(1, 2, 3)), "base_prior_sp")
  expect_error(f(base_prior_prev = c(-1, 2)), "base_prior_prev")
  expect_error(f(base_prior_prev = NULL), "base_prior_prev")
  expect_error(f(design = "xx"))
  expect_error(f(use_in = "analysis"))
})

test_that("print method runs and shows key information", {
  p <- do.call(bam_prior_from_pilot, vap_args(use_in = "design"))
  txt <- paste(capture.output(print(p)), collapse = "\n")
  expect_match(txt, "a0")
  expect_match(txt, "cohort")
  expect_match(txt, "design prior", ignore.case = TRUE)
  expect_match(txt, "analysis prior", ignore.case = TRUE)
  expect_match(txt, "a \\+ b")
  expect_identical(class(suppressWarnings(print(p)) |> invisible()),
                   "dta_bam_priors")
})

test_that("integration with bam_sample_size(priors = obj)", {
  skip_if_not("priors" %in% names(formals(bam_sample_size)),
              "bam_sample_size() has no `priors` argument yet")
  p <- do.call(bam_prior_from_pilot, vap_args(use_in = "design"))
  common <- list(delta_se = 0.15, delta_sp = 0.15, target_assurance = 0.8,
                 method = "exact", N_range = c(40, 400))
  r1 <- do.call(bam_sample_size, c(list(priors = p), common))
  r2 <- do.call(bam_sample_size,
                c(list(prior_se = p$prior_se, prior_sp = p$prior_sp,
                       prior_prev = p$prior_prev,
                       analysis_prior_se = p$analysis_prior_se,
                       analysis_prior_sp = p$analysis_prior_sp), common))
  keys <- c("N_total", "target_reached", "joint_assurance",
            "max_assurance_evaluated", "N_at_max_assurance")
  expect_identical(r1[keys], r2[keys])
})
