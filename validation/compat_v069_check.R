# Identity check against v0.6.9.
#
# With the new arguments of bam_sample_size() at their default values, every
# numerical output must be bit-for-bit identical to that of v0.6.9.
#
# Usage (each step in its own Rscript process):
#
#   Rscript compat_v069_check.R <library_path> <output_directory>
#       Loads dtasamplesize from that library, runs the fixed battery with the
#       default arguments and saves <output_directory>/compat_battery.rds
#       (results, warnings, sessionInfo() and package version).
#
#   Rscript compat_v069_check.R --compare <dir_v069> <dir_new> <log>
#       Compares the two RDS files with identical(), field by field and
#       excluding only `call`, and writes the result to <log>. Fields that
#       exist only in the new version are listed separately (not compared).
#
# Battery: bam_sample_size(), the published P1 and P2; the 12 scenarios of
# validation/python/LOCKED_confirmatory_grid.json (read with jsonlite, not
# modified); and 3 short calls in Monte Carlo mode. Scenario C12 is infeasible
# by design: in both versions N_total must be NA. mc_validate_buderer(): 6
# calls in Monte Carlo mode (MV1 to MV6, with the arguments of version 0.6.9;
# the argument `method` and the field `method_type` exist only since 0.7.0 and
# are not compared). Each field is compared with identical(), excluding only
# `call`.

args <- commandArgs(trailingOnly = TRUE)

script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[1], winslash = "/")) else getwd()
}

# ---------------------------------------------------------------- compare --
if (length(args) >= 1 && identical(args[1], "--compare")) {
  stopifnot(length(args) == 4)
  dir_old <- args[2]; dir_new <- args[3]; log_file <- args[4]
  old <- readRDS(file.path(dir_old, "compat_battery.rds"))
  new <- readRDS(file.path(dir_new, "compat_battery.rds"))
  lines <- character()
  say <- function(...) lines <<- c(lines, paste0(...))
  say("compat_v069_check: comparison written ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  say("old: dtasamplesize ", old$package_version, " (", old$lib, ")")
  say("new: dtasamplesize ", new$package_version, " (", new$lib, ")")
  say("")
  stopifnot(identical(names(old$results), names(new$results)))
  n_ok <- 0L; n_bad <- 0L; n_fields <- 0L; extra_shown <- character()
  for (id in names(old$results)) {
    ro <- old$results[[id]]; rn <- new$results[[id]]
    fo <- setdiff(names(ro), "call")
    fn <- setdiff(names(rn), "call")
    extra <- setdiff(fn, fo)
    missing_new <- setdiff(fo, fn)
    bad <- character()
    for (f in intersect(fo, fn)) {
      n_fields <- n_fields + 1L
      if (!identical(ro[[f]], rn[[f]])) bad <- c(bad, f)
    }
    wo <- old$warnings[[id]]; wn <- new$warnings[[id]]
    warn_same <- identical(wo, wn)
    ok <- length(bad) == 0L && length(missing_new) == 0L
    if (ok) n_ok <- n_ok + 1L else n_bad <- n_bad + 1L
    ntot <- function(r) if (!is.null(r$N_total)) r$N_total else r$n_total
    say(sprintf("%-14s %s  N_total old=%s new=%s | fields compared=%d | warnings identical=%s",
                id, if (ok) "IDENTICAL" else "DIFFERENT",
                format(ntot(ro)), format(ntot(rn)),
                length(intersect(fo, fn)), warn_same))
    if (length(bad)) say("    differing fields: ", paste(bad, collapse = ", "))
    if (length(missing_new)) say("    fields missing in new: ", paste(missing_new, collapse = ", "))
    kind <- if (startsWith(id, "MV")) "mc_validate_buderer" else "bam_sample_size"
    if (length(extra) && !(kind %in% extra_shown)) {
      extra_shown <- c(extra_shown, kind)
      say("    (", kind, ": fields only in new version, not compared: ",
          paste(extra, collapse = ", "), ")")
    }
    if (!warn_same) {
      say("    warnings old: ", paste(wo, collapse = " || "))
      say("    warnings new: ", paste(wn, collapse = " || "))
    }
  }
  c12o <- old$results[["C12"]]$N_total; c12n <- new$results[["C12"]]$N_total
  say("")
  say("C12 (infeasible by design): N_total old = ", format(c12o), ", new = ", format(c12n),
      "; both NA: ", is.na(c12o) && is.na(c12n))
  say(sprintf("SUMMARY: %d/%d calls identical field by field (%d fields compared, `call` excluded); %d different.",
              n_ok, n_ok + n_bad, n_fields, n_bad))
  say(if (n_bad == 0L && is.na(c12o) && is.na(c12n)) "RESULT: PASS" else "RESULT: FAIL")
  writeLines(lines, log_file)
  cat(paste(lines, collapse = "\n"), "\n")
  quit(status = if (n_bad == 0L) 0L else 1L)
}

# -------------------------------------------------------------------- run --
stopifnot(length(args) == 2)
lib <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
out_dir <- args[2]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressMessages(library(dtasamplesize, lib.loc = lib))
pv <- as.character(utils::packageVersion("dtasamplesize", lib.loc = lib))
cat("dtasamplesize", pv, "loaded from", lib, "\n")

grid_path <- file.path(script_dir(), "python", "LOCKED_confirmatory_grid.json")
grid <- jsonlite::fromJSON(grid_path, simplifyVector = TRUE)
sc <- grid$scenarios

run_one <- function(args) {
  w <- character()
  res <- withCallingHandlers(
    do.call(bam_sample_size, args),
    warning = function(cond) {
      w <<- c(w, conditionMessage(cond))
      invokeRestart("muffleWarning")
    }
  )
  list(res = res, warnings = w)
}

# n_range/B: arguments of the per-arm legacy fields; they are fixed so that
# the battery is the same in both versions (the NEW arguments, which are not
# passed, keep their default values).
legacy <- list(n_range = 20:500, B = 2000)

calls <- list(
  P1 = c(list(prior_se = c(17, 3), prior_sp = c(2, 2), prior_prev = c(4, 16),
              delta_se = 0.14, delta_sp = 0.10, target_assurance = 0.80), legacy),
  P2 = c(list(prior_se = c(17, 3), prior_sp = c(18, 2), prior_prev = c(4, 16),
              delta_se = 0.14, delta_sp = 0.10, target_assurance = 0.80), legacy)
)
for (i in seq_len(nrow(sc))) {
  s <- sc[i, ]
  calls[[s$id]] <- c(list(
    prior_se = as.numeric(s$prior_se[[1]]), prior_sp = as.numeric(s$prior_sp[[1]]),
    prior_prev = as.numeric(s$prior_prev[[1]]),
    delta_se = s$delta_se, delta_sp = s$delta_sp,
    target_assurance = s$target, alpha_ci = s$level,
    N_max = as.integer(s$N_hi)
  ), legacy)
}
# Tres llamadas cortas en modo Monte Carlo (B y N_range pequenos)
mc_base <- list(prior_se = c(17, 3), prior_sp = c(2, 2), prior_prev = c(4, 16),
                method = "monte_carlo")
calls$MC1 <- c(mc_base, list(delta_se = 0.14, delta_sp = 0.10, B = 400,
                             n_range = 20:120, N_range = seq(300, 700, by = 50), seed = 2026))
calls$MC2 <- c(mc_base, list(delta_se = 0.30, delta_sp = 0.30, target_assurance = 0.85,
                             B = 300, n_range = 20:60, N_range = NULL, seed = 7))
calls$MC3 <- c(mc_base, list(delta_se = 0.20, delta_sp = 0.20, B = 500, n_range = 20:100,
                             N_range = 40:90, seed = 99, alpha_ci = 0.90))

# Six calls of mc_validate_buderer() in Monte Carlo mode (the arguments that
# exist in 0.6.9; `method` is not passed, its default is the Monte Carlo mode
# of 0.6.9).
bud_calls <- list(
  MV1 = list(),
  MV2 = list(Se = 0.85, d = 0.07, B = 20000),
  MV3 = list(Se = 0.98, d = 0.05, B = 4000, ci_method = "wald"),
  MV4 = list(Se = 0.70, d = 0.10, n_diseased = 60, B = 3000,
             ci_method = "wilson", seed = 11),
  MV5 = list(Se = 0.90, d = 0.03, alpha = 0.10, B = 2500),
  MV6 = list(Se = 0.60, d = 0.08, n_diseased = 150, B = 1000, seed = 7)
)

run_bud <- function(args) {
  w <- character()
  res <- withCallingHandlers(
    do.call(mc_validate_buderer, args),
    warning = function(cond) {
      w <<- c(w, conditionMessage(cond))
      invokeRestart("muffleWarning")
    }
  )
  list(res = res, warnings = w)
}

results <- list(); warns <- list()
for (id in names(calls)) {
  t0 <- Sys.time()
  r <- run_one(calls[[id]])
  results[[id]] <- unclass(r$res)
  warns[[id]] <- r$warnings
  cat(sprintf("%-4s N_total = %s  (%.1f s)\n", id, format(r$res$N_total),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
for (id in names(bud_calls)) {
  t0 <- Sys.time()
  r <- run_bud(bud_calls[[id]])
  results[[id]] <- unclass(r$res)
  warns[[id]] <- r$warnings
  calls[[id]] <- bud_calls[[id]]
  cat(sprintf("%-4s n_total = %s  (%.1f s)\n", id, format(r$res$n_total),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

saveRDS(list(
  results = results, warnings = warns, calls = calls,
  package_version = pv, lib = lib, R_version = R.version.string,
  session_info = utils::capture.output(print(sessionInfo())),
  written = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
), file.path(out_dir, "compat_battery.rds"))
cat("saved", file.path(out_dir, "compat_battery.rds"), "\n")
