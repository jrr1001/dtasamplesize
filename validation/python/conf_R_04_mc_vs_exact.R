# CONFIRMATORY extra: does the PACKAGE's own Monte Carlo engine agree with the
# PACKAGE's closed form at the published decision points?  This tests the
# closed-form derivation against a simulation of the generative model that the
# package itself implements -- an internal-consistency check that is
# independent of my Python reference.
# MC budget fixed by the locked rule: B = 160000 -> MCSE ~ 0.0010 near p=0.8.
suppressPackageStartupMessages(library(dtasamplesize))
options(dtasamplesize.warn_small_B = FALSE)

B <- 160000
cases <- list(
  list(tag = "P1", psp = c(2, 2),  N = 678, exact = 0.8003489948),
  list(tag = "P2", psp = c(18, 2), N = 672, exact = 0.8002692084))

cat(sprintf("%-4s %-5s %-6s %13s %13s %9s %7s\n",
            "case", "seed", "N", "pkg exact", "pkg MC", "mcse", "z"))
for (cs in cases) for (sd in c(70117, 70118, 70119)) {
  common <- list(prior_se = c(17, 3), prior_sp = cs$psp,
                 prior_prev = c(4, 16), delta_se = 0.14, delta_sp = 0.10,
                 target_assurance = 0.80, alpha_ci = 0.95,
                 n_range = 20:200, N_range = cs$N, seed = sd)
  e <- suppressWarnings(do.call(bam_sample_size,
                                c(common, list(method = "exact", B = 1000))))
  m <- suppressWarnings(do.call(bam_sample_size,
                                c(common, list(method = "monte_carlo", B = B))))
  z <- (m$joint_assurance - e$joint_assurance) / m$assurance_mcse
  cat(sprintf("%-4s %-5d %-6d %13.10f %13.10f %9.5f %+7.2f\n",
              cs$tag, sd, cs$N, e$joint_assurance, m$joint_assurance,
              m$assurance_mcse, z))
  stopifnot(abs(e$joint_assurance - cs$exact) < 1e-9)
}
cat("\n(pkg exact values matched the published figures to <1e-9 in all rows)\n")
