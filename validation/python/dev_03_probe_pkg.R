# DEVELOPMENT step 3: probe how to query the package for the assurance at an
# arbitrary single N, so that N-1 / N / N+1 can be read off the package too.
suppressPackageStartupMessages(library(dtasamplesize))
options(dtasamplesize.warn_small_B = FALSE)
cat("package version:", as.character(packageVersion("dtasamplesize")), "\n\n")

args_common <- list(prior_se = c(17, 3), prior_sp = c(2, 2),
                    prior_prev = c(4, 16), delta_se = 0.14, delta_sp = 0.10,
                    target_assurance = 0.80, method = "exact",
                    n_range = 20:400, B = 1000)

probe <- function(Nr, label) {
  r <- tryCatch(
    withCallingHandlers(
      do.call(bam_sample_size, c(args_common, list(N_range = Nr))),
      warning = function(w) { cat("   [warn]", conditionMessage(w), "\n")
                              invokeRestart("muffleWarning") }),
    error = function(e) { cat("   [error]", conditionMessage(e), "\n"); NULL })
  if (is.null(r)) return(invisible(NULL))
  cat(sprintf("%-28s N_total=%s  joint_assurance=%.10f  mcse=%s\n",
              label, format(r$N_total), r$joint_assurance,
              format(r$assurance_mcse)))
  invisible(r)
}

cat("-- single-N probes around the published crossing (P1, N*=678) --\n")
for (N in c(676, 677, 678, 679)) probe(N, paste0("N_range=", N))

cat("\n-- full-range search --\n")
probe(600:750, "N_range=600:750")

cat("\n-- names of the returned object --\n")
r <- probe(600:750, "again")
cat(paste(names(r), collapse = ", "), "\n")
