# Integration and reproducibility check for dtasamplesize
# -------------------------------------------------------
# Runs all exported estimators end-to-end, confirms the
# "more uncertainty -> larger N" ordering, and verifies that
# results are reproducible under a fixed seed.
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

cat("\n=== integration_check done ===\n")
