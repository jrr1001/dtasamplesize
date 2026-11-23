#!/usr/bin/env Rscript
# comparator_calls.R -- calls to external comparators and the two integer-grid
# checks quoted in the manuscript.
#
# Usage:   Rscript validation/comparator_calls.R <output_directory>
# Writes ONLY to <output_directory> (comparator_calls_results.csv and
# comparator_calls_sessioninfo.txt); nothing else is written anywhere.
#
# Conventions (named explicitly, rule 11):
#   * power.diagnostic.test(): one arm at a time; delta is the half-width
#     margin on the lower confidence bound (marginal assurance); prevalence is
#     a fixed, known constant.
#   * SampleSizeDiagnostics(): classical (Buderer) total sample size; w is the
#     half-width.
#   * joint_sample_size(): joint Se/Sp event, half-width targets delta_se and
#     delta_sp, Monte Carlo (B draws, seed), grid N_range.
#   * Wilson enumeration: two-sided 95% Wilson interval, FULL widths
#     (upper - lower) compared with 2*delta_se and 2*delta_sp, with the
#     declared tolerance TOL_WIDTH = 1e-12 added to the target
#     (width <= target + TOL_WIDTH). Joint event = both widths attained, at
#     point values Se, Sp, prevalence: D ~ Binomial(N, prev);
#     TP | D ~ Binomial(D, Se); TN | N - D ~ Binomial(N - D, Sp); an arm with
#     0 subjects never attains the target.
#   * Selection: "first crossing" = smallest N in the grid with A(N) >= 0.80;
#     "stable crossing" = smallest N such that A(N') >= 0.80 for every N' in
#     the grid from N upward.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: Rscript validation/comparator_calls.R <output_directory>", call. = FALSE)
out_dir <- args[1]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("comparator_calls.R\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("R:", R.version.string, "\n")
pv <- function(p) if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else NA_character_
for (p in c("dtasamplesize", "MKpower", "MKmisc", "SampleSizeDiagnostics")) {
  cat(sprintf("  %-24s %s\n", p, pv(p)))
}
cat("\n")

res <- list()
add <- function(id, quantity, manuscript, value, note = "") {
  res[[length(res) + 1L]] <<- data.frame(id = id, quantity = quantity, manuscript = manuscript,
                                         value = as.character(value), note = note,
                                         stringsAsFactors = FALSE)
}
TOL_WIDTH <- 1e-12

## (a) power.diagnostic.test: MKpower and MKmisc -------------------------------
cat("== (a) power.diagnostic.test(sens = 0.85, delta = 0.07, power = 0.80, prev = 0.20) ==\n")
for (pkg in c("MKpower", "MKmisc")) {
  if (!requireNamespace(pkg, quietly = TRUE)) { cat(pkg, "not installed\n"); next }
  f <- getExportedValue(pkg, "power.diagnostic.test")
  r <- suppressMessages(f(sens = 0.85, delta = 0.07, power = 0.80, prev = 0.20))
  cat("-- ", pkg, " ", pv(pkg), " --\n", sep = "")
  print(r)
  cat("names(result):", paste(names(r), collapse = ", "), "\n")
  nums <- unlist(r[vapply(r, function(z) is.numeric(z) && length(z) == 1L, logical(1))])
  cat("numeric scalars:", paste(names(nums), format(nums), sep = " = ", collapse = "; "), "\n\n")
  add(paste0("a_", pkg), paste0(pkg, "::power.diagnostic.test"), "216 diseased / 864 non-diseased",
      paste(names(nums), format(nums), sep = "=", collapse = " | "),
      paste0("version ", pv(pkg), "; call sens=0.85, delta=0.07, power=0.80, prev=0.20; other arguments at defaults"))
}

## (b) SampleSizeDiagnostics -----------------------------------------------------
cat("== (b) SampleSizeDiagnostics::SampleSizeDiagnostics(sn = 0.85, sp = 0.90, p = 0.20, w = 0.07, CI = 0.95) ==\n")
if (requireNamespace("SampleSizeDiagnostics", quietly = TRUE)) {
  r <- SampleSizeDiagnostics::SampleSizeDiagnostics(sn = 0.85, sp = 0.90, p = 0.20, w = 0.07, CI = 0.95)
  print(r)
  cat("\n")
  add("b_SSD", "SampleSizeDiagnostics::SampleSizeDiagnostics", "500",
      paste(capture.output(print(r)), collapse = " ; "),
      paste0("version ", pv("SampleSizeDiagnostics"),
             "; exact call SampleSizeDiagnostics(sn=0.85, sp=0.90, p=0.20, w=0.07, CI=0.95); w is the half-width"))
} else cat("SampleSizeDiagnostics not installed\n")

## (c) joint_sample_size, integer grid, B = 20000, seed = 2026 ----------------------
cat("== (c) dtasamplesize::joint_sample_size(Se = 0.85, Sp = 0.90, delta_se = 0.07, delta_sp = 0.05, prev = 0.20,\n",
    "      design = 'cohort', target_prob = 0.80, N_range = 100:800, B = 20000, seed = 2026) ==\n", sep = "")
if (requireNamespace("dtasamplesize", quietly = TRUE)) {
  jr <- suppressWarnings(dtasamplesize::joint_sample_size(
    Se = 0.85, Sp = 0.90, delta_se = 0.07, delta_sp = 0.05, prev = 0.20, design = "cohort",
    target_prob = 0.80, N_range = 100:800, B = 20000, seed = 2026))
  cat("n_total:", jr$n_total, "\n")
  cat("joint_prob_se_sp:", format(jr$joint_prob_se_sp, digits = 10), "\n")
  cat("joint_prob_mcse:", format(jr$joint_prob_mcse, digits = 6), "\n")
  cat("search_type:", jr$search_type, "\n\n")
  add("c_joint", "joint_sample_size n_total (N_range = 100:800, B = 20000, seed = 2026)", "579",
      jr$n_total,
      sprintf("dtasamplesize %s; joint_prob_se_sp = %.6f; mcse = %.6f; search_type = %s",
              pv("dtasamplesize"), jr$joint_prob_se_sp, jr$joint_prob_mcse, jr$search_type))
}

## (d) deterministic enumeration of the Wilson joint event --------------------------
cat("== (d) Deterministic enumeration of the Wilson joint event (point Se, Sp, prevalence) ==\n")
wilson_width <- function(x, n, z = qnorm(0.975)) {
  # full width (upper - lower) of the two-sided Wilson score interval
  ph <- x / n; z2 <- z^2
  2 * z * sqrt(ph * (1 - ph) / n + z2 / (4 * n^2)) / (1 + z2 / n)
}
p_arm_ok <- function(n, p, target) {
  # P(Wilson full width <= target + TOL_WIDTH), X ~ Binomial(n, p); n = 0 never attains
  if (n < 1L) return(0)
  x <- 0:n
  sum(dbinom(x, n, p)[wilson_width(x, n) <= target + TOL_WIDTH])
}
arm_table <- function(nmax, p, target) vapply(0:nmax, p_arm_ok, numeric(1), p = p, target = target)
joint_A <- function(N, pse_tab, psp_tab, prev) {
  D <- 0:N
  sum(dbinom(D, N, prev) * pse_tab[D + 1L] * psp_tab[N - D + 1L])
}
SE <- 0.85; SP <- 0.90; PREV <- 0.20; DSE <- 0.07; DSP <- 0.05; TARGET <- 0.80
grid <- 100:800
pse_tab <- arm_table(max(grid), SE, 2 * DSE)
psp_tab <- arm_table(max(grid), SP, 2 * DSP)
A <- vapply(grid, joint_A, numeric(1), pse_tab = pse_tab, psp_tab = psp_tab, prev = PREV)
ok <- A >= TARGET
first <- grid[which(ok)[1]]
stable <- grid[min(which(vapply(seq_along(grid), function(i) all(ok[i:length(grid)]), logical(1))))]
cat(sprintf("Se = %.2f, Sp = %.2f, prevalence = %.2f; full-width targets %.2f (Se) and %.2f (Sp); TOL_WIDTH = %g\n",
            SE, SP, PREV, 2 * DSE, 2 * DSP, TOL_WIDTH))
cat("First integer N with A(N) >= 0.80:", first, " A =", format(A[grid == first], digits = 10), "\n")
cat("A(N-1) at first crossing:", format(A[grid == first - 1L], digits = 10), "\n")
cat("Stable crossing (A >= 0.80 for all N' >= N up to 800):", stable, " A =", format(A[grid == stable], digits = 10), "\n")
cat("Grid points with A >= 0.80 followed by a drop below 0.80:", sum(ok[-length(ok)] & !ok[-1]), "\n")
cat("First of seq(100, 800, 10) with A >= 0.80:", grid[which(grid %% 10 == 0 & ok)[1]], "\n\n")
add("d_enum_first", "Wilson joint-event enumeration: first N with A >= 0.80", "578", first,
    sprintf("A(first) = %.7f; A(N-1) = %.7f; tolerance %g on widths; n = 0 arm fails",
            A[grid == first], A[grid == first - 1L], TOL_WIDTH))
add("d_enum_stable", "Wilson joint-event enumeration: stable crossing", "578", stable,
    sprintf("A(stable) = %.7f", A[grid == stable]))

res <- do.call(rbind, res)
utils::write.csv(res, file.path(out_dir, "comparator_calls_results.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "comparator_calls_sessioninfo.txt"))
cat("Wrote results to", normalizePath(out_dir), "\n")
