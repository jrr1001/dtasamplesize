# CONFIRMATORY run, package side.
# Reads the LOCKED grid; queries dtasamplesize::bam_sample_size(method="exact").
# Writes conf_package.json.  The locked grid is NOT modified.
suppressPackageStartupMessages({library(dtasamplesize); library(jsonlite)})
options(dtasamplesize.warn_small_B = FALSE)

grid <- fromJSON("LOCKED_confirmatory_grid.json", simplifyDataFrame = FALSE)
sc <- grid$scenarios

call_pkg <- function(s, N_range, target = NULL, dmult = 1) {
  warns <- character(0)
  r <- withCallingHandlers(
    bam_sample_size(
      prior_se = unlist(s$prior_se), prior_sp = unlist(s$prior_sp),
      prior_prev = unlist(s$prior_prev),
      delta_se = s$delta_se * dmult, delta_sp = s$delta_sp * dmult,
      target_assurance = if (is.null(target)) s$target else target,
      alpha_ci = s$level, n_range = 20:200, B = 1000,
      N_range = N_range, method = "exact", seed = 2026),
    warning = function(w) { warns <<- c(warns, conditionMessage(w))
                            invokeRestart("muffleWarning") })
  list(N_total = r$N_total, joint_assurance = r$joint_assurance,
       mcse = r$assurance_mcse, assurance_method = r$assurance_method,
       warnings = warns)
}

out <- list()
for (s in sc) {
  t0 <- Sys.time()
  cat(sprintf("[%s] searching N_range=%d:%d ... ", s$id, s$N_lo, s$N_hi))
  flush.console()
  main <- call_pkg(s, s$N_lo:s$N_hi)
  # the package's own assurance at N-1, N, N+1 (single-N N_range probes)
  trio <- list()
  for (d in c(-1, 0, 1)) {
    Nq <- main$N_total + d
    if (Nq >= 1) {
      q <- call_pkg(s, Nq)
      trio[[as.character(d)]] <- list(N = Nq, assurance = q$joint_assurance)
    }
  }
  el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("N_total=%s  A=%.10f  (%.0fs)  warns=%d\n",
              format(main$N_total), main$joint_assurance, el,
              length(main$warnings)))
  out[[s$id]] <- list(id = s$id, main = main, trio = trio, secs = el)
}

write_json(out, "conf_package.json", auto_unbox = TRUE, digits = 16,
           pretty = TRUE)
cat("\nwrote conf_package.json\n")
cat("sessionInfo R:", R.version.string, "\n")
cat("dtasamplesize:", as.character(packageVersion("dtasamplesize")), "\n")
