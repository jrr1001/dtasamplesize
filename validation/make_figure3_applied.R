#!/usr/bin/env Rscript
# =============================================================================
# make_figure3_applied.R
#
# Data and figure of Figure 3 of the R1 manuscript: exact assurance curves of
# the applied example (Panbio validation study planned from the published
# pilot of Albert et al. 2021). The curves are
# computed with the package, with the same calls as section 7 of
# validation/applied_example.R:
#   bam_prior_from_pilot(tp = 43, fn = 11, tn = 358, fp = 0, base priors
#   Beta(1, 1), a0, use_in) and bam_assurance_curve(priors, N_range = 2:1450,
#   delta_se = 0.10, delta_sp = 0.05, width_type = "lower_half_width",
#   alpha_ci = 0.95, target_assurance = 0.80).
# Curves: a0 = 0.5, 0.25 and 1 with use_in = "design" (flat analysis prior), and
# a0 = 0.5 with use_in = "both" (the updated prior also used in the analysis).
#
# Usage
#   Rscript validation/make_figure3_applied.R <out_dir> --applied=<dir> [--name=Figure_3]
#   <out_dir> is the only place where the script writes (default:
#   validation/output/ inside the package tree). --applied is the output
#   directory of validation/applied_example.R: its main table and
#   sensitivity table are read for the anchors; the log gives path and SHA-256.
#   --name is the base name of the figure files (Figure_3 in the main text or
#   Figure_S3.2 in Additional file 3; default Figure_3).
#   Libraries: DTA_PKG_LIB and DTA_EXTRA_LIB, if set, are prepended to
#   .libPaths(). Single core. The calculation (4 curves of 1449 points, exact
#   enumeration, O(N^2) per curve) takes several minutes: run it in the background.
#
# Conventions: total sample size N; one-sided 95% lower semi-widths 0.10 (Se)
# and 0.05 (Sp), joint objective, assurance target 0.80; exact (no Monte
# Carlo). Reference lines: 927 (sample size planned by Bulilete et al. with the
# frequentist formula) and 1369 (participants enrolled), both from Bulilete et
# al. 2021.
#
# Outputs: R1_Figure3_applied_curves.csv (curve, a0, use_in, N, assurance),
#   <name>.pdf (vector), <name>.tiff (600 dpi, 170 mm wide, no internal title),
#   make_figure3_applied.log.
# Anchors: first N with A(N) >= 0.80 for the four curves (497, 520, 480, 258),
# and A(497) of the main curve equal to the main table of applied_example.R (<= 1e-10).
# =============================================================================

options(width = 140, stringsAsFactors = FALSE, dtasamplesize.warn_small_B = FALSE)
libs <- c(Sys.getenv("DTA_PKG_LIB"), Sys.getenv("DTA_EXTRA_LIB")); libs <- libs[nzchar(libs)]
if (length(libs)) .libPaths(c(libs, .libPaths()))

locate_script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) == 1) return(dirname(normalizePath(sub("^--file=", "", f[1]), winslash = "/")))
  if (basename(getwd()) == "validation") return(normalizePath(getwd(), winslash = "/"))
  if (dir.exists("validation")) return(normalizePath(file.path(getwd(), "validation"), winslash = "/"))
  stop("Cannot locate the validation/ folder; run with Rscript validation/make_figure3_applied.R")
}
script_dir <- locate_script_dir(); pkg_root <- dirname(script_dir)
args_all <- commandArgs(trailingOnly = TRUE)
named <- list(); pos <- character(0)
for (a in args_all) {
  if (grepl("^--[A-Za-z0-9_]+=", a)) named[[sub("^--([A-Za-z0-9_]+)=.*$", "\\1", a)]] <- sub("^--[A-Za-z0-9_]+=", "", a)
  else pos <- c(pos, a)
}
out_dir <- if (length(pos) >= 1 && nzchar(pos[1])) pos[1] else file.path(script_dir, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_dir <- normalizePath(out_dir, winslash = "/")
if (is.null(named$applied) || !dir.exists(named$applied)) stop("Missing or invalid --applied=<dir>")
FIG_NAME <- if (!is.null(named$name) && nzchar(named$name)) named$name else "Figure_3"

logcon <- file(file.path(out_dir, "make_figure3_applied.log"), open = "wt", encoding = "UTF-8")
sink(logcon, split = TRUE)
suppressMessages(library(dtasamplesize))
git_commit <- tryCatch({
  h <- system2("git", c("-C", shQuote(pkg_root), "rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE)
  n <- length(system2("git", c("-C", shQuote(pkg_root), "status", "--porcelain"), stdout = TRUE, stderr = FALSE))
  sprintf("%s%s", h[1], if (n > 0) sprintf(" (+%d uncommitted paths)", n) else "")
}, error = function(e) "not available", warning = function(w) "not available")
cat(sprintf("%s | %s | dtasamplesize %s at %s | tree commit %s | script make_figure3_applied.R | out %s\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), R.version.string,
            as.character(utils::packageVersion("dtasamplesize")), find.package("dtasamplesize"), git_commit, out_dir))
T_START <- Sys.time()

rd <- function(file) {
  p <- normalizePath(file.path(named$applied, file), winslash = "/", mustWork = FALSE)
  if (!file.exists(p)) stop("Input not found: ", p)
  cat(sprintf("INPUT --applied %s sha256=%s\n", p, unname(tools::sha256sum(p))))
  utils::read.csv(p, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
}
ANCH <- new.env(); ANCH$n <- 0L; ANCH$pass <- 0L; ANCH$fail <- character(0)
anchor <- function(id, desc, obs, exp, tol = 0) {
  ok <- isTRUE(length(obs) == length(exp) && all(abs(as.numeric(obs) - exp) <= tol))
  ANCH$n <- ANCH$n + 1L
  if (ok) ANCH$pass <- ANCH$pass + 1L else ANCH$fail <- c(ANCH$fail, id)
  cat(sprintf("ANCHOR %s %s %s (%s)\n", id, desc, if (ok) "PASS" else "**CHECK**", paste(format(obs, digits = 12, trim = TRUE), collapse = ", ")))
  invisible(ok)
}

main_tab <- rd("applied_example_main_table.csv")
sens_tab <- rd("applied_example_sensitivity.csv")
pilot <- list(tp = 43, fn = 11, tn = 358, fp = 0)   # Albert et al. 2021; G3b_conteos_verificados.md L7-L9
flat <- c(1, 1)
pri <- function(a0, use_in) do.call(bam_prior_from_pilot, c(pilot, list(
  base_prior_se = flat, base_prior_sp = flat, base_prior_prev = flat, a0 = a0, use_in = use_in)))
COMMON <- list(delta_se = 0.10, delta_sp = 0.05, width_type = "lower_half_width", alpha_ci = 0.95, target_assurance = 0.80)
NR <- 2:1450
CURVES <- data.frame(
  curve = c("a0 = 0.5, use_in design (main)", "a0 = 0.25, use_in design", "a0 = 1, use_in design", "a0 = 0.5, use_in both"),
  a0 = c(0.5, 0.25, 1, 0.5), use_in = c("design", "design", "design", "both"), stringsAsFactors = FALSE)
cat("Conventions: one-sided 95%, lower semi-widths 0.10 (Se) and 0.05 (Sp), joint, target 0.80, total N, exact; N_range 2:1450\n")
res <- vector("list", nrow(CURVES))
for (i in seq_len(nrow(CURVES))) {
  t0 <- Sys.time()
  cv <- do.call(bam_assurance_curve, c(list(priors = pri(CURVES$a0[i], CURVES$use_in[i]), N_range = NR), COMMON))
  res[[i]] <- cv
  cat(sprintf("curve %-38s N_first = %s  (%.1f s)\n", CURVES$curve[i], as.character(cv$N_first),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
curves_df <- do.call(rbind, lapply(seq_len(nrow(CURVES)), function(i)
  data.frame(curve = CURVES$curve[i], a0 = CURVES$a0[i], use_in = CURVES$use_in[i],
             N = res[[i]]$curve$N, assurance = res[[i]]$curve$assurance, stringsAsFactors = FALSE)))
utils::write.csv(curves_df, file.path(out_dir, "R1_Figure3_applied_curves.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("WROTE R1_Figure3_applied_curves.csv (%d rows)\n", nrow(curves_df)))
# A(N) is not monotone in N: count its decreases (A(N) < A(N - 1)) and report the largest one, per curve.
cat("Decreases of A(N) along each curve (N_range 2:1450):\n")
for (i in seq_len(nrow(CURVES))) {
  ci <- curves_df[curves_df$curve == CURVES$curve[i], ]
  ci <- ci[order(ci$N), ]
  dA <- diff(ci$assurance)
  dn <- which(dA < 0)
  if (length(dn) > 0L) {
    k <- dn[which.min(dA[dn])]
    cat(sprintf("  curve %-38s decreases = %d; largest decrease = %.3e between N = %d and N = %d; N of the decreases from %d to %d\n",
                CURVES$curve[i], length(dn), -dA[k], ci$N[k], ci$N[k + 1L], ci$N[min(dn)], ci$N[max(dn) + 1L]))
  } else {
    cat(sprintf("  curve %-38s decreases = 0\n", CURVES$curve[i]))
  }
}

N_first <- vapply(res, function(cv) as.numeric(cv$N_first), numeric(1))
anchor("FIG3-01", "first crossing of 0.80: main; a0 0.25; a0 1; both (T5-06, T5-10, T5-07)", N_first, c(497, 520, 480, 258))
sc <- function(s) sens_tab$N_total[sens_tab$scenario == s]
anchor("FIG3-01b", "crossings equal to the applied_example.R tables (main, a0_0.25_design, a0_1_design, both_onesided_joint)",
       N_first, c(main_tab$N_total[main_tab$scenario == "design_onesided_joint"], sc("a0_0.25_design"), sc("a0_1_design"),
                  main_tab$N_total[main_tab$scenario == "both_onesided_joint"]))
c1 <- curves_df[curves_df$curve == CURVES$curve[1], ]
A497 <- c1$assurance[c1$N == 497]
anchor("T5-06", "A(497) of the main curve equals A_N of the main table (<= 1e-10)", A497,
       main_tab$A_N[main_tab$scenario == "design_onesided_joint"], tol = 1e-10)
for (i in seq_len(nrow(CURVES))) {
  ci <- curves_df[curves_df$curve == CURVES$curve[i], ]
  anchor(sprintf("FIG3-02-%d", i), sprintf("curve %d: A at its N_first >= 0.80 and A(N_first - 1) < 0.80", i),
         c(ci$assurance[ci$N == N_first[i]] >= 0.80 - 1e-12, ci$assurance[ci$N == N_first[i] - 1] < 0.80), c(1, 1))
}

## ---- figure ---------------------------------------------------------------
N_PLAN <- 927; N_ENR <- 1369; N_MAIN <- 497
COL <- c("#000000", "#0072B2", "#D55E00", "#009E73")   # Okabe-Ito
LTY <- c(1, 2, 3, 4); LWD <- c(2.2, 1.8, 1.8, 1.8)
W_IN <- 170 / 25.4; H_IN <- W_IN * 5 / 7
draw <- function() {
  par(mar = c(4.3, 4.6, 0.8, 0.8), cex.axis = 0.8, cex.lab = 0.85, las = 1, mgp = c(2.9, 0.7, 0))
  plot(NA, xlim = c(0, 1450), ylim = c(0, 1), xlab = "Total sample size (N)", ylab = "Exact assurance, A(N)", xaxs = "i", yaxs = "i")
  abline(h = 0.80, col = "grey45", lty = 5, lwd = 1)
  abline(v = c(N_PLAN, N_ENR), col = "grey35", lty = 3, lwd = 1)
  for (i in rev(seq_len(nrow(CURVES)))) {
    ci <- curves_df[curves_df$curve == CURVES$curve[i], ]
    lines(ci$N, ci$assurance, col = COL[i], lty = LTY[i], lwd = LWD[i])
  }
  points(N_MAIN, A497, pch = 19, cex = 1.1)
  text(N_MAIN + 20, 0.74, sprintf("N = %d", N_MAIN), pos = 4, cex = 0.78)
  text(N_PLAN - 15, 0.45, "planned (frequentist)
N = 927", adj = c(1, 0.5), cex = 0.72)
  text(N_ENR - 15, 0.45, "enrolled
N = 1369", adj = c(1, 0.5), cex = 0.72)
  # the label and the legend sit away from the dotted lines at N = 927 and N = 1369
  text(15, 0.82, "target 0.80", adj = c(0, 0), cex = 0.72, col = "grey30")
  lg <- legend("bottomright", legend = CURVES$curve, col = COL, lty = LTY, lwd = LWD, cex = 0.72, bty = "n", seg.len = 3, plot = FALSE)
  legend(x = N_PLAN - 15 - lg$rect$w, y = 0.03 + lg$rect$h, legend = CURVES$curve, col = COL, lty = LTY, lwd = LWD, cex = 0.72, bty = "n", seg.len = 3)
}
emit <- function(name) {
  grDevices::tiff(file.path(out_dir, paste0(name, ".tiff")), width = W_IN, height = H_IN, units = "in", res = 600, compression = "lzw")
  draw(); grDevices::dev.off()
  ok <- tryCatch({ grDevices::cairo_pdf(file.path(out_dir, paste0(name, ".pdf")), width = W_IN, height = H_IN); TRUE },
                 error = function(e) { grDevices::pdf(file.path(out_dir, paste0(name, ".pdf")), width = W_IN, height = H_IN); FALSE })
  draw(); grDevices::dev.off()
  cat(sprintf("WROTE %s.pdf and %s.tiff (%.1f mm x %.1f mm, 600 dpi; cairo_pdf = %s)\n", name, name, 170, H_IN * 25.4, ok))
}
emit(FIG_NAME)

cat(sprintf("\nelapsed: %.1f s\n", as.numeric(difftime(Sys.time(), T_START, units = "secs"))))
if (length(ANCH$fail) == 0L) {
  cat(sprintf("SUMMARY: %d anchors, all PASS\n", ANCH$n))
} else {
  cat(sprintf("SUMMARY: %d anchors, %d PASS, %d **CHECK**: %s\n", ANCH$n, ANCH$pass, length(ANCH$fail), paste(ANCH$fail, collapse = ", ")))
}
sink(); close(logcon)
if (length(ANCH$fail)) quit(save = "no", status = 1)
