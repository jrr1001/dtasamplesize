# HISTORY (kept as a labelled historical note, not current validation
# evidence): the installed package was originally rebuilt (0.6.2 -> 0.6.3)
# by a parallel agent while the confirmatory grid was running, which is why
# conf_package.json's baseline predates later releases.  Re-verify the
# headline numbers against whatever is installed NOW (as of Lote 02, a
# privately installed dtasamplesize 0.6.6), so the report states a single
# reproducible provenance.  Compares against the baseline values recorded in
# conf_package.json.
suppressPackageStartupMessages({library(dtasamplesize); library(jsonlite)})
options(dtasamplesize.warn_small_B = FALSE)
v <- as.character(packageVersion("dtasamplesize"))
cat("dtasamplesize now installed:", v, "\n\n")

grid <- fromJSON("LOCKED_confirmatory_grid.json", simplifyDataFrame = FALSE)
prev <- fromJSON("conf_package.json", simplifyDataFrame = FALSE)

run1 <- function(s) suppressWarnings(bam_sample_size(
  prior_se = unlist(s$prior_se), prior_sp = unlist(s$prior_sp),
  prior_prev = unlist(s$prior_prev), delta_se = s$delta_se,
  delta_sp = s$delta_sp, target_assurance = s$target, alpha_ci = s$level,
  n_range = 20:200, B = 1000, N_range = s$N_lo:s$N_hi,
  method = "exact", seed = 2026))

cat("-- published scenarios --\n")
pub <- list(list(tag="P1", psp=c(2,2), N=678, A=0.800349, Am1=0.799685),
            list(tag="P2", psp=c(18,2), N=672, A=0.800269, Am1=NA))
for (p in pub) {
  r <- suppressWarnings(bam_sample_size(
    prior_se = c(17,3), prior_sp = p$psp, prior_prev = c(4,16),
    delta_se = 0.14, delta_sp = 0.10, target_assurance = 0.80,
    n_range = 20:200, N_range = 300:900, method = "exact"))
  r1 <- suppressWarnings(bam_sample_size(
    prior_se = c(17,3), prior_sp = p$psp, prior_prev = c(4,16),
    delta_se = 0.14, delta_sp = 0.10, target_assurance = 0.80,
    n_range = 20:200, N_range = p$N - 1, method = "exact"))
  cat(sprintf("%s: N_total=%d (claim %d)  A=%.10f (claim %.6f)  A(N-1)=%.10f\n",
              p$tag, r$N_total, p$N, r$joint_assurance, p$A, r1$joint_assurance))
}

cat("\n-- locked confirmatory grid, re-run on", v, "--\n")
cat(sprintf("%-5s %8s %8s %8s %15s %15s\n",
            "id", "N(baseline)", "N(now)", "same", "A(baseline)", "A(now)"))
allsame <- TRUE
for (s in grid$scenarios) {
  r <- run1(s)
  old <- prev[[s$id]]$main
  same <- identical(as.numeric(r$N_total), as.numeric(old$N_total)) &&
          isTRUE(all.equal(r$joint_assurance, old$joint_assurance, tolerance = 1e-12))
  allsame <- allsame && same
  cat(sprintf("%-5s %8s %8s %8s %15.10f %15.10f\n", s$id,
              format(old$N_total), format(r$N_total), same,
              old$joint_assurance, r$joint_assurance))
}
cat("\nALL CONFIRMATORY RESULTS UNCHANGED ACROSS THE REBUILD:", allsame, "\n")
