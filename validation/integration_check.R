# Integration and reproducibility check for dtasamplesize
# -------------------------------------------------------
# Runs all exported estimators end-to-end, prints the sample size each one
# returns so the values can be inspected, verifies that results are
# reproducible under a fixed seed, and asserts the "more uncertainty ->
# larger N" ordering for bam_sample_size(method = "exact") across the Se
# prior triple Beta(34,6), Beta(17,3), Beta(8.5,1.5) (all other arguments
# default). That ordering check is executable code, in the dedicated section
# below headed "Ordering check", not just a claim in this header: the three
# N_total values are computed and printed by the script, and the script
# stop()s if they are not strictly increasing.
#
# Run after installing the package:  R -f validation/integration_check.R

library(dtasamplesize)
options(width = 100)
sep <- function(t) cat("\n----------", t, "----------\n")

sep("Integration: all estimators run")
r1 <- mc_validate_buderer(B = 800, seed = 1)
r2 <- suppressWarnings(bam_sample_size(B = 800, seed = 1, n_range = 30:300))
r3 <- suppressWarnings(joint_sample_size(B = 800, seed = 1, N_range = seq(100, 700, 50)))
r4 <- ss_imperfect_ref(B = 0)
r5 <- ss_adaptive_prevalence(B = 400, seed = 1, prev_true_range = c(0.20, 0.30))
r6 <- suppressWarnings(ss_net_benefit(B = 800, seed = 1))
r8 <- suppressWarnings(ss_unified(B = 400, seed = 1, N_range = seq(200, 800, 50),
                                  delta_auc = 0, check_nb = FALSE))
cat("buderer Se-only :", buderer_n(0.85, 0.07), "\n")
cat("BAM n_diseased  :", r2$n_diseased, " N_total_median:", r2$N_total_median, "\n")
cat("Joint N         :", r3$n_total, "\n")
cat("Imperfect N     :", r4$n_total, "\n")
cat("Adaptive N_init :", r5$N_initial_adj, "\n")
cat("NetBenefit Ncons:", r6$N_conservative, "\n")
cat("Unified N       :", r8$n_total, "\n")
if (requireNamespace("timeROC", quietly = TRUE)) {
  r7 <- suppressWarnings(ss_time_dependent_roc(B = 30, seed = 1,
          N_range = seq(150, 350, 50), censoring_rates = 0.2))
  cat("timeROC N       :", r7$n_total, "(B=30 demo)\n")
}

sep("Reproducibility: identical results on repeat (same seed)")
a <- suppressWarnings(ss_net_benefit(B = 1000, seed = 2026, pt_range = c(0.1, 0.3, 0.5)))
b <- suppressWarnings(ss_net_benefit(B = 1000, seed = 2026, pt_range = c(0.1, 0.3, 0.5)))
cat("net_benefit identical:", identical(a$N_by_pt, b$N_by_pt), "\n")
c1 <- suppressWarnings(bam_sample_size(B = 1000, seed = 2026, n_range = 30:200))
c2 <- suppressWarnings(bam_sample_size(B = 1000, seed = 2026, n_range = 30:200))
cat("bam identical N_total:", identical(c1$N_total_median, c2$N_total_median), "\n")

# ============================================================================
# Ordering check: more uncertainty -> larger N (exact mode)
# ----------------------------------------------------------------------------
# Computes bam_sample_size(prior_se = p, method = "exact") (all other args
# default) for p in Beta(34,6), Beta(17,3), Beta(8.5,1.5) -- a fixed prior
# mean (E[Se] = 0.85 in all three) with decreasing prior sample size, i.e.
# increasing prior uncertainty about Se -- and asserts that the resulting
# N_total values are strictly increasing. method = "exact" is deterministic
# (see bam_sample_size()'s @details), so this check has no sampling error and
# needs no seed. This is the only place in this script that verifies the
# ordering with code rather than asserting it in prose.
# ============================================================================
sep("Ordering check: more uncertainty -> larger N (exact mode)")
priors_se <- list(c(34, 6), c(17, 3), c(8.5, 1.5))
n_ordering <- vapply(priors_se, function(p) {
  bam_sample_size(prior_se = p, method = "exact")$N_total
}, numeric(1))
cat("N_total for prior_se = Beta(34,6), Beta(17,3), Beta(8.5,1.5):",
    paste(n_ordering, collapse = ", "), "\n")
if (!all(diff(n_ordering) > 0)) {
  stop("Ordering check FAILED: N_total is not strictly increasing as the Se ",
       "prior widens (Beta(34,6) -> Beta(17,3) -> Beta(8.5,1.5)). Got N_total = ",
       paste(n_ordering, collapse = ", "), ".")
}
cat("PASS: more uncertainty -> larger N (exact mode)\n")
# ============================================================================

cat("\n=== integration_check done ===\n")
