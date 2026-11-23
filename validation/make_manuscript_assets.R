# Generate the manuscript figures and tables for dtasamplesize
# ---------------------------------------------------------------------------
# This script regenerates every figure and every data-derived table used in
# the accompanying article, from the installed package and nothing else: no
# number here is copied by hand from a previous run. Figures are written as
# a vector PDF (BMC's preferred format for line art) and a 600 dpi TIFF;
# tables are written as CSV. EVERYTHING is written to ONE output directory:
# the first command-line argument or, if none is given, validation/output/
# inside the package tree (created if it does not exist). The script never
# derives a path from folders above the package, so it cannot write into
# the manuscript folders. Re-running it regenerates all of it, because
# every quantity is either a closed-form / exact calculation or a Monte
# Carlo simulation under a fixed seed.
#
# SCOPE. The article describes bam_sample_size(),
# bam_assurance_curve(), mc_validate_buderer() and joint_sample_size();
# ss_unified(), the AUC gate inside it, ss_net_benefit(), ss_imperfect_ref(),
# ss_adaptive_prevalence() and ss_time_dependent_roc() remain in the package
# but are not described in the article, and this script builds nothing for
# them. It writes:
#   Figure 1 (R1): EXACT curve of P(full CI width <= 0.14) against the
#     number of diseased participants n (Se = 0.85, d = 0.07, n = 40..200),
#     Wald (truncated) and Wilson, from mc_validate_buderer(method =
#     "exact") evaluated n by n. It replaces the earlier Monte Carlo figure.
#   Figure 2 (R1): EXACT joint assurance A(N) of the worked example, from
#     bam_assurance_curve().
#   Table 1 (the estimators described in the article), the feature
#   matrix of v0.6.9 (Table 4 of v0.6.9, read from table3_feature_matrix_source.csv; the R1 matrix is Table S3.3e of Additional file 3, built by make_manuscript_tables_R1.R) and the comparison of the three
#   methods under common assumptions (table4_method_comparison.csv; Table 2 of v0.6.9 and of the R1, the R1 version being built by make_manuscript_tables_R1.R).
# The cross-validation table (Table 3 of v0.6.9; Table S2.5 of Additional file 2 in the R1) is written separately by
# cross_validation.R -- see validation/README.md.
#
# ANCHOR CHECKS. Both figures are guarded by anchor checks (values quoted
# by the article). If any fails, the script stops before drawing: it never
# adjusts an anchor.
#
# Usage, from the package root (installed package, any working directory):
#   Rscript validation/make_manuscript_assets.R [output_dir]
# Default output_dir: validation/output/ (inside the package tree).

suppressMessages(library(dtasamplesize))
options(width = 120, dtasamplesize.warn_small_B = FALSE)

## ---------------------------------------------------------------------
## Locate the package tree from this script's own path (Rscript --file=),
## or from the working directory when sourced. The output directory is the
## first command-line argument; the default lives INSIDE the package tree.
## ---------------------------------------------------------------------
locate_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  if (basename(getwd()) == "validation") return(normalizePath(getwd()))
  if (dir.exists("validation")) return(normalizePath(file.path(getwd(), "validation")))
  stop(
    "Cannot determine the location of make_manuscript_assets.R. Run it via\n",
    "  Rscript validation/make_manuscript_assets.R [output_dir]\n",
    "from the dtasamplesize package root, or from inside validation/.",
    call. = FALSE
  )
}

script_dir <- locate_script_dir()                  # .../dtasamplesize/validation
pkg_root   <- dirname(script_dir)                  # .../dtasamplesize
cli_args   <- commandArgs(trailingOnly = TRUE)
out_dir    <- if (length(cli_args) >= 1 && nzchar(cli_args[1])) cli_args[1] else
  file.path(script_dir, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_dir    <- normalizePath(out_dir, winslash = "/")
figures_dir <- out_dir   # figures and tables share the one output directory
assets_dir  <- out_dir

cat("R version:          ", R.version.string, "\n")
cat("dtasamplesize:      ", as.character(utils::packageVersion("dtasamplesize")), "\n")
cat("Output directory:   ", out_dir, "\n\n")

## Anchor bookkeeping: an anchor that fails stops the script at once; the
## passed anchors are listed in the summary line printed at the end.
anchor_log <- character(0)
check_anchor <- function(name, pass, detail) {
  cat(sprintf("ANCHOR %-34s %s  (%s)\n", name,
              if (isTRUE(pass)) "PASS" else "FAIL", detail))
  if (!isTRUE(pass)) {
    stop("Anchor check failed: ", name, " (", detail, "). Refusing to ",
         "continue: the anchor is not adjusted.", call. = FALSE)
  }
  anchor_log[[length(anchor_log) + 1L]] <<- paste0(name, "=PASS")
  invisible(TRUE)
}

## ---------------------------------------------------------------------
## Figure geometry (H-08): 170 mm output width, keeping the 7:5 (w:h)
## aspect ratio previously hard-coded as 7 in x 5 in (177.8 mm x 127.0 mm).
## 170 / 177.8 = 0.9561, so this is a ~4.4% linear shrink, not a re-layout;
## base-graphics text is set in POINTS via cex (see draw_fig1()/draw_fig2(),
## cex = 0.78-0.80 against the device default pointsize = 12), which is an
## absolute unit independent of the device's physical width in inches, so
## the smallest figure text stays ~9.4-9.6 pt at the new size, above the
## 8 pt minimum, unchanged from before.
## ---------------------------------------------------------------------
FIG_WIDTH_MM  <- 170
FIG_WIDTH_IN  <- FIG_WIDTH_MM / 25.4
FIG_HEIGHT_IN <- FIG_WIDTH_IN * 5 / 7
cat(sprintf("Figure size: %.4f in x %.4f in (%.1f mm x %.1f mm)\n\n",
            FIG_WIDTH_IN, FIG_HEIGHT_IN, FIG_WIDTH_MM, FIG_HEIGHT_IN * 25.4))

## ---------------------------------------------------------------------
## emit(): write one figure as a 600 dpi TIFF and a vector PDF, with no
## embedded title (the caption belongs in the manuscript, not the image).
## cairo_pdf keeps line art vector and embeds fonts; if cairo support is
## unavailable on the host R build, fall back to the base pdf() device
## (still vector, just without guaranteed font embedding).
## ---------------------------------------------------------------------
emit <- function(name, draw, w_in, h_in) {
  grDevices::tiff(file.path(figures_dir, paste0(name, ".tiff")),
                  width = w_in, height = h_in, units = "in", res = 600,
                  compression = "lzw")
  draw(); grDevices::dev.off()

  ok <- tryCatch({
    grDevices::cairo_pdf(file.path(figures_dir, paste0(name, ".pdf")),
                          width = w_in, height = h_in)
    TRUE
  }, error = function(e) {
    grDevices::pdf(file.path(figures_dir, paste0(name, ".pdf")),
                    width = w_in, height = h_in)
    FALSE
  })
  draw(); grDevices::dev.off()
  invisible(ok)
}

## ---------------------------------------------------------------------
## run_checked(): evaluate a sample-size search and fail loudly if it hit
## the "did not converge within N_range" condition -- the function
## returned silently once the range ran out, at a value that is not a
## real crossing of the target assurance. An advisory warning (e.g. "the
## joint assurance is close to its own Monte Carlo noise") is not that
## failure mode: it is printed, not raised, because it is informative
## about how tight the margin is, not a sign that the search gave up.
## ---------------------------------------------------------------------
run_checked <- function(label, expr) {
  msgs <- character(0)
  result <- withCallingHandlers(
    expr,
    warning = function(w) {
      msgs[[length(msgs) + 1]] <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  )
  non_convergence <- grep(
    "No N in N_range|No n in n_range|Consider expanding",
    msgs, value = TRUE
  )
  if (length(non_convergence)) {
    stop(sprintf(
      "%s: search did not converge within the supplied range -- %s",
      label, paste(non_convergence, collapse = " | ")
    ), call. = FALSE)
  }
  if (length(msgs)) {
    cat(sprintf("  [%s] advisory warning(s) (not a convergence failure):\n", label))
    for (m in msgs) cat("    - ", m, "\n", sep = "")
  }
  result
}

## Shared design parameters, used throughout (see the article's Methods):
## a moderately informative Se prior, a prevalence prior centred at 0.20,
## and the half-width precision targets Se +/- 0.07 and Sp +/- 0.05 at 80%
## assurance. PRIOR_SP is the HARMONIZED, informative specificity prior
## (Beta(18, 2), mean 0.90) used only in the comparison of methods (table4_method_comparison.csv), where it is compared
## against the classical and Monte Carlo methods under matched point
## assumptions (Se = 0.85, Sp = 0.90, prevalence = 0.20); it is NOT the
## package's own default Sp prior. Figure 2 instead uses the package's
## actual default, vague Sp prior, Beta(2, 2) -- see PRIOR_SP_VAGUE below
## and the article's Results/Limitations for why the two differ (N = 678
## under the vague default vs. N = 672 under the harmonized informative
## prior).
PRIOR_SE   <- c(17, 3)
PRIOR_SP   <- c(18, 2)                 # harmonized (informative), method comparison only
PRIOR_PREV <- c(4, 16)                 # Beta(4, 16): mean 0.20
DELTA_SE   <- 0.07                     # half-width, as used by the classical
DELTA_SP   <- 0.05                     # and Monte Carlo methods in the method comparison
TARGET_ASSURANCE <- 0.80
B_MAIN     <- 20000
SEED       <- 2026

E_SE   <- PRIOR_SE[1]   / sum(PRIOR_SE)     # 0.85
E_SP   <- PRIOR_SP[1]   / sum(PRIOR_SP)     # 0.90 (harmonized prior mean)
E_PREV <- PRIOR_PREV[1] / sum(PRIOR_PREV)   # 0.20

## ===========================================================================
## FIGURE 1 -- EXACT assurance of the classical Buderer sample size
## ===========================================================================
## R1 version (replaces the Monte Carlo figure): P(full CI width <= 2 d =
## 0.14) for the sensitivity interval as a function of the number of
## diseased participants n, for Se = 0.85 and d = 0.07, n = 40..200, computed
## EXACTLY (binomial enumeration, no simulation) with
## mc_validate_buderer(method = "exact") evaluated for each n. Two interval
## types: Wald truncated to [0, 1] and Wilson score; both two-sided 95%.
## The curve is NOT monotone (sawtooth, from the discreteness of x/n): it
## is plotted as computed, without smoothing.
FIG1_SE <- 0.85; FIG1_D <- 0.07
ns_fig1 <- 40:200
exact_p <- function(n) {
  r <- mc_validate_buderer(Se = FIG1_SE, d = FIG1_D, n_diseased = n,
                           ci_method = "both", method = "exact")$results
  stats::setNames(r$P_width_target, r$ci_method)[c("wald", "wilson")]
}
P1 <- t(vapply(ns_fig1, exact_p, numeric(2)))
colnames(P1) <- c("wald", "wilson")
p_wald_fig1   <- P1[, "wald"]
p_wilson_fig1 <- P1[, "wilson"]
n_bud_fig1 <- buderer_n(FIG1_SE, FIG1_D)

first_ge <- function(p, thr = 0.80) ns_fig1[which(p >= thr)[1]]
# smallest n such that P(n') >= thr for every n' from n up to 200
stay_ge <- function(p, thr = 0.80) {
  below <- which(p < thr)
  if (!length(below)) ns_fig1[1] else ns_fig1[max(below) + 1L]
}
first_wald   <- first_ge(p_wald_fig1);   first_wilson <- first_ge(p_wilson_fig1)
stay_wald    <- stay_ge(p_wald_fig1);    stay_wilson  <- stay_ge(p_wilson_fig1)
p100_wald    <- p_wald_fig1[ns_fig1 == 100]
p100_wilson  <- p_wilson_fig1[ns_fig1 == 100]

check_anchor("F1 Buderer n = 100", n_bud_fig1 == 100,
             sprintf("buderer_n = %d", n_bud_fig1))
check_anchor("F1 P(100), Wald = 0.5683", isTRUE(round(p100_wald, 4) == 0.5683),
             sprintf("%.10f", p100_wald))
check_anchor("F1 P(100), Wilson = 0.5683", isTRUE(round(p100_wilson, 4) == 0.5683),
             sprintf("%.10f", p100_wilson))
check_anchor("F1 first n P>=0.80, Wald = 114", isTRUE(first_wald == 114L),
             sprintf("first n = %s", first_wald))
check_anchor("F1 first n P>=0.80, Wilson = 114", isTRUE(first_wilson == 114L),
             sprintf("first n = %s", first_wilson))
check_anchor("F1 stays >=0.80 to 200, Wald = 117", isTRUE(stay_wald == 117L),
             sprintf("n = %s", stay_wald))
check_anchor("F1 stays >=0.80 to 200, Wilson = 114", isTRUE(stay_wilson == 114L),
             sprintf("n = %s", stay_wilson))

fig1_data <- data.frame(n_diseased = ns_fig1, P_wald = p_wald_fig1,
                        P_wilson = p_wilson_fig1)
write.csv(fig1_data, file.path(out_dir, "figure1_exact_curve.csv"),
          row.names = FALSE)

COL_WALD <- "#1f4e79"; COL_WILSON <- "#d9822b"
draw_fig1 <- function() {
  par(mar = c(4.3, 4.6, 1.2, 1))
  plot(ns_fig1, p_wald_fig1, type = "l", lwd = 2, col = COL_WALD,
       xlab = "Number of diseased participants (n)",
       ylab = "P(full CI width \u2264 0.14)", ylim = c(0, 1), las = 1)
  lines(ns_fig1, p_wilson_fig1, lwd = 2, col = COL_WILSON, lty = 5)
  abline(h = 0.80, lty = 3, col = "grey40")
  abline(v = n_bud_fig1, lty = 2, col = "#c00000")
  points(n_bud_fig1, p100_wald, pch = 19, col = "#c00000")
  text(n_bud_fig1 + 3, 0.25,
       sprintf("Buderer n=%d\nP=%.4f", n_bud_fig1, p100_wald),
       col = "#c00000", cex = 0.78, pos = 4)
  legend("bottomright", bty = "n", cex = 0.8,
         legend = c("Wald (truncated), exact", "Wilson, exact",
                    "0.80 target", "Buderer n"),
         col = c(COL_WALD, COL_WILSON, "grey40", "#c00000"),
         lty = c(1, 5, 3, 2), lwd = c(2, 2, 1, 1))
}
emit("Figure_1", draw_fig1, FIG_WIDTH_IN, FIG_HEIGHT_IN)
cat(sprintf(paste0("Figure 1: Buderer n=%d -> P=%.4f (Wald), %.4f (Wilson); ",
                   "first n with P>=0.80: %d (Wald), %d (Wilson); smallest n ",
                   "with P>=0.80 up to 200: %d (Wald), %d (Wilson)\n"),
            n_bud_fig1, p100_wald, p100_wilson, first_wald, first_wilson,
            stay_wald, stay_wilson))

## ===========================================================================
## FIGURE 2 -- EXACT joint assurance A(N) of the worked example
## ===========================================================================
## Se ~ Beta(17, 3), Sp ~ Beta(2, 2) (the package's own vague default, NOT
## the harmonized Beta(18, 2) of the method comparison), prevalence ~ Beta(4, 16) (NOT the
## package default c(6, 14), so it is passed explicitly), total CI widths
## 0.14 (Se) and 0.10 (Sp), level 0.95, target assurance 0.80. The curve is
## computed with bam_assurance_curve(), the exported function that gives
## A(N) at every N with the same calculation as bam_sample_size(method =
## "exact") and no early stop.
PRIOR_SP_VAGUE <- c(2, 2)   # bam_sample_size() default prior_sp
curve_fig2 <- bam_assurance_curve(
  prior_se = PRIOR_SE, prior_sp = PRIOR_SP_VAGUE, prior_prev = PRIOR_PREV,
  delta_se = 2 * DELTA_SE, delta_sp = 2 * DELTA_SP,
  target_assurance = TARGET_ASSURANCE, N_range = 300:900)
N_grid_fig2   <- curve_fig2$curve$N
assur_fig2    <- curve_fig2$curve$assurance
n_cross_fig2  <- curve_fig2$N_first
assur_at_cross_fig2 <- assur_fig2[N_grid_fig2 == n_cross_fig2]
assur_at_prev_fig2  <- assur_fig2[N_grid_fig2 == (n_cross_fig2 - 1L)]

check_anchor("F2 N_first = 678", isTRUE(n_cross_fig2 == 678L),
             sprintf("N_first = %s", n_cross_fig2))
check_anchor("F2 A(678) = 0.8003489948",
             isTRUE(abs(assur_at_cross_fig2 - 0.8003489948) < 1e-10),
             sprintf("%.10f", assur_at_cross_fig2))
check_anchor("F2 A(677) = 0.7996848824",
             isTRUE(abs(assur_at_prev_fig2 - 0.7996848824) < 1e-10),
             sprintf("%.10f", assur_at_prev_fig2))

write.csv(curve_fig2$curve, file.path(out_dir, "figure2_exact_curve.csv"),
          row.names = FALSE)

draw_fig2 <- function() {
  par(mar = c(4.3, 4.6, 1.2, 1))
  plot(N_grid_fig2, assur_fig2, type = "l", lwd = 2, col = COL_WALD,
       xlab = "Total sample size (N)",
       ylab = "Joint assurance for Se and Sp", ylim = c(0, 1), las = 1)
  abline(h = TARGET_ASSURANCE, lty = 3, col = "grey40")
  abline(v = n_cross_fig2, lty = 2, col = "#c00000")
  points(n_cross_fig2, assur_at_cross_fig2, pch = 19, col = "#c00000")
  text(n_cross_fig2 + 14, 0.28,
       sprintf("N=%d\nassurance=%.4f", n_cross_fig2, assur_at_cross_fig2),
       col = "#c00000", cex = 0.78, pos = 4)
  legend("topleft", bty = "n", cex = 0.8,
         legend = c("Exact joint assurance", "0.80 target",
                    sprintf("N=%d (first N reaching 0.80)", n_cross_fig2)),
         col = c(COL_WALD, "grey40", "#c00000"),
         lty = c(1, 3, 2), lwd = c(2, 1, 1))
}
emit("Figure_2", draw_fig2, FIG_WIDTH_IN, FIG_HEIGHT_IN)
cat(sprintf(
  "Figure 2: crossing at N=%d, assurance=%.10f (N-1=%d, assurance=%.10f)\n",
  n_cross_fig2, assur_at_cross_fig2, n_cross_fig2 - 1L, assur_at_prev_fig2))

## ===========================================================================
## TABLE 4 -- the three surviving methods for the joint precision of Se
##            and Sp, under a common set of point assumptions
## ===========================================================================
## Reduced from the previous release's five-method comparison: only the
## three methods the article still describes (Buderer classical, BAM exact
## mode with harmonized priors, and joint_sample_size()) remain. The
## imperfect-reference and unified-framework rows are removed along with
## the v0.6.9 Figures 3-4 above (see the SCOPE note at the top of this file). Each
## method contributes the N it actually produces under the article's common
## assumptions (Se = 0.85, Sp = 0.90, prevalence = 0.20; for the Bayesian
## method, the harmonized prior means Se ~ Beta(17, 3), Sp ~ Beta(18, 2),
## prevalence ~ Beta(4, 16)) -- not a value read off a figure, since neither
## remaining figure sweeps this comparison.

## Method 1: Buderer (classical, deterministic) -- per-arm requirement
## combined across Se and Sp via the common prevalence.
se_arm_n <- buderer_n(0.85, 0.07)
sp_arm_n <- buderer_n(0.90, 0.05)
N_buderer <- ceiling(max(se_arm_n / 0.20, sp_arm_n / 0.80))
h_buderer <- list(N = N_buderer, assurance = NA_real_,
                   assurance_type = "deterministic (normal-approximation target width; not a simulated assurance)")

## Method 2: BAM, exact joint Beta-Binomial search, full-width deltas,
## HARMONIZED (informative) Sp prior -- see the shared-parameters note
## above for why this differs from Figure 2's vague default.
bam_res <- run_checked("Method comparison: BAM (exact, harmonized priors)",
  bam_sample_size(prior_se = PRIOR_SE, prior_sp = PRIOR_SP, prior_prev = PRIOR_PREV,
                   delta_se = 2 * DELTA_SE, delta_sp = 2 * DELTA_SP,
                   target_assurance = TARGET_ASSURANCE, method = "exact",
                   B = B_MAIN, seed = SEED))
h_bam <- list(N = bam_res$N_total, assurance = bam_res$joint_assurance,
               assurance_type = "exact Bayesian (Beta-Binomial), closed form, no Monte Carlo error")

## Method 3: joint Se + Sp precision by direct Monte Carlo, prospective-
## cohort design (joint_sample_size()'s default). AUC and delta_auc are
## left at the function's own defaults (AUC = 0.90, delta_auc = 0.05):
## joint_sample_size() always evaluates a deterministic AUC gate
## internally (see ?joint_sample_size, @details), but the article no
## longer describes that gate, and it is not part of the joint_prob_se_sp
## estimand compared here (see the "notes" column below).
joint_res <- run_checked("Method comparison: joint Se/Sp (joint_sample_size)",
  joint_sample_size(Se = E_SE, Sp = E_SP, delta_se = DELTA_SE, delta_sp = DELTA_SP,
                     prev = E_PREV, design = "cohort",
                     target_prob = TARGET_ASSURANCE, B = B_MAIN, seed = SEED))
h_joint <- list(N = joint_res$n_total, assurance = joint_res$joint_prob_se_sp,
                 assurance_type = "Monte Carlo (cohort design) for Se & Sp jointly")

table4 <- data.frame(
  method = c("Buderer (classical)", "BAM (exact mode, harmonized priors)",
             "Joint Se/Sp (joint_sample_size)"),
  N = c(h_buderer$N, h_bam$N, h_joint$N),
  assurance_achieved = c(h_buderer$assurance, h_bam$assurance, h_joint$assurance),
  assurance_type = c(h_buderer$assurance_type, h_bam$assurance_type, h_joint$assurance_type),
  Se = c(0.85, E_SE, E_SE),
  Sp = c(0.90, E_SP, E_SP),
  prevalence = c(0.20, E_PREV, E_PREV),
  prior_se = c(NA, "17,3", NA),
  prior_sp = c(NA, "18,2", NA),
  prior_prev = c(NA, "4,16", NA),
  delta_se = c(DELTA_SE, 2 * DELTA_SE, DELTA_SE),
  delta_sp = c(DELTA_SP, 2 * DELTA_SP, DELTA_SP),
  B = c(NA, B_MAIN, B_MAIN),
  seed = c(NA, SEED, SEED),
  notes = c(
    "N = ceiling(max(buderer_n(0.85,0.07)/0.20, buderer_n(0.90,0.05)/0.80))",
    "delta_se/delta_sp given as full CI widths (2x the half-width used elsewhere); prior_sp harmonized to Beta(18,2) (package/article-default prior_sp for Figure 2 is the vague Beta(2,2))",
    sprintf(paste("AUC and delta_auc left at joint_sample_size()'s own defaults",
                   "(%.2f and %.2f); the function evaluates a deterministic AUC",
                   "gate internally but the article no longer describes it, and",
                   "it is not part of the joint_prob_se_sp estimand reported here"),
            joint_res$AUC, 0.05)
  ),
  stringsAsFactors = FALSE
)
write.csv(table4, file.path(assets_dir, "table4_method_comparison.csv"), row.names = FALSE)
cat("\nMethod comparison (table4_method_comparison.csv; Table 2 of v0.6.9 and of the R1, the R1 version being built by make_manuscript_tables_R1.R):\n")
print(table4[, c("method", "N", "assurance_achieved")], row.names = FALSE)

## A stale copy of the previous (five-method) version of this table may be
## left over from an earlier release under its old filename; remove it so
## manuscript_assets/ does not carry two versions of the same table.
old_table4_path <- file.path(assets_dir, "table_harmonized_comparison.csv")
if (file.exists(old_table4_path)) {
  file.remove(old_table4_path)
  cat("Removed stale", old_table4_path, "(superseded by table4_method_comparison.csv)\n")
}

## ===========================================================================
## TABLE 1 -- the estimators described in this article
## ===========================================================================
## The article's Table 1 caption is explicit that this table lists only the
## three estimators actually DESCRIBED in this article (mc_validate_buderer,
## bam_sample_size, joint_sample_size); the package's other exported
## helpers and its further, currently experimental modules are named in the
## caption text itself, not as additional table rows (see manuscript.md).
## This differs from the previous release, which listed every function
## exported in NAMESPACE. The three purposes below are the article's own
## Table 1 wording verbatim, not the functions' \title{} (those are
## one-line Rd titles written for a different audience and do not match
## the article's phrasing -- e.g. bam_sample_size()'s \title{} is
## "Bayesian Assurance Method for DTA Sample Size", not the fuller
## description the article's table gives).
##
## The guard below still checks against NAMESPACE, not to derive the
## purposes, but to catch silent drift in the other direction: if one of
## these three functions were ever removed from NAMESPACE, this script
## must fail loudly rather than publish a table for a function the
## installed package no longer exports.
ns_lines <- readLines(file.path(pkg_root, "NAMESPACE"), warn = FALSE)
export_lines <- grep("^export\\(", ns_lines, value = TRUE)
fn_names <- sub("^export\\(([^)]+)\\)$", "\\1", export_lines)

table1_functions <- c("mc_validate_buderer", "bam_sample_size", "joint_sample_size")
table1_purposes <- c(
  mc_validate_buderer = "Monte Carlo validation of the classical formula (diagnostic, not a planner)",
  bam_sample_size     = "Bayesian assurance sample size for the joint precision of Se and Sp, under uncertain prevalence (closed-form exact mode available)",
  joint_sample_size   = "Joint Se and Sp precision by direct Monte Carlo simulation, cohort design"
)
missing_from_namespace <- setdiff(table1_functions, fn_names)
if (length(missing_from_namespace)) {
  stop("Table 1: the following functions are described in the article's ",
       "Table 1 but are no longer exported by the installed package's ",
       "NAMESPACE: ", paste(missing_from_namespace, collapse = ", "),
       ". Check that R-package/dtasamplesize/NAMESPACE and the article ",
       "are still in sync.", call. = FALSE)
}

table1 <- data.frame(Function = table1_functions,
                      Purpose = unname(table1_purposes[table1_functions]),
                      stringsAsFactors = FALSE)
write.csv(table1, file.path(assets_dir, "table1_functions.csv"), row.names = FALSE)
cat("\nTable 1: ", nrow(table1), " estimators described in the article\n", sep = "")

## ===========================================================================
## TABLE 3 -- feature comparison with existing R packages
## ===========================================================================
## Unchanged: this matrix is curated by hand -- it records a judgment about
## what each comparator package's public functions actually do, not
## something that can be derived mechanically -- and lives as a separate,
## versioned CSV (table3_feature_matrix_source.csv) precisely so that the
## curation is visible and diffable on its own, independent of this script.
## The block below only re-emits it, together with a provenance note. The
## capability content (which comparator does what) has not changed in this
## revision; only the cell-text wording shown in the manuscript's rendered
## table is more compact than this source CSV's fuller annotations.
table3_source_path <- file.path(script_dir, "table3_feature_matrix_source.csv")
table3 <- read.csv(table3_source_path, stringsAsFactors = FALSE, check.names = FALSE)

## Guard against a silent CSV-parsing corruption: a stray unquoted comma
## inside a field (e.g. "No (foo, bar)") makes a data row parse to more
## fields than the header has, which makes read.csv() promote the first
## column to row names instead of a data column -- every column after it
## then shifts left by one, and write.csv(..., row.names = FALSE) discards
## the row names, so the shift produces no error or warning at all. Check
## for that failure mode explicitly, right after reading the source, so a
## malformed source CSV stops the build instead of quietly shipping a
## mis-attributed comparison table.
expected_cols <- c("Capability", "presize", "pROC", "MKpower", "epiR",
                    "MKmisc", "SampleSizeDiagnostics", "dtasamplesize")
if (ncol(table3) != length(expected_cols)) {
  stop("Feature-matrix source has ", ncol(table3), " columns, expected ",
       length(expected_cols), " (", paste(expected_cols, collapse = ", "),
       "). This usually means a data row in ", table3_source_path,
       " has an unescaped comma inside an unquoted field, which shifts",
       " read.csv()'s column count.", call. = FALSE)
}
if (!identical(colnames(table3)[1], "Capability")) {
  stop("Feature-matrix source: first column is named '", colnames(table3)[1],
       "', not 'Capability'. read.csv() likely promoted the Capability",
       " column to row names because a data row in ", table3_source_path,
       " parsed to more fields than the header -- check for an unquoted",
       " comma inside a cell.", call. = FALSE)
}
if (!identical(rownames(table3), as.character(seq_len(nrow(table3))))) {
  stop("Feature-matrix source: row.names(table3) are not the default 1..",
       nrow(table3), " -- read.csv() promoted a data column to row names,",
       " which means some row in ", table3_source_path, " has more fields",
       " than the header. Check every row for an unescaped comma inside",
       " an unquoted field and quote it.", call. = FALSE)
}
if (any(is.na(table3)) || any(table3 == "")) {
  bad <- which(is.na(table3) | table3 == "", arr.ind = TRUE)
  stop("Feature-matrix source has ", nrow(bad), " empty/NA cell(s) after parsing ",
       "(row, col): ", paste(sprintf("(%d,%d)", bad[, 1], bad[, 2]),
       collapse = "; "), ". A well-formed matrix should have no empty",
       " cells -- this is the signature of a row that parsed to the wrong",
       " number of fields. Check ", table3_source_path,
       " for an unescaped comma inside an unquoted field.", call. = FALSE)
}

write.csv(table3, file.path(assets_dir, "table3_feature_matrix.csv"), row.names = FALSE)

## Provenance: the version of each comparator package that was actually
## read (its exported function list and help pages) while curating the
## matrix above, plus whatever version (if any) is installed in the
## environment that is running this script right now, so a later reader
## can tell whether the matrix might be stale relative to their own
## installation. Versions are never invented: a package absent from the
## running environment is reported as such, not silently skipped.
compare_pkgs <- c("presize", "pROC", "MKpower", "epiR", "MKmisc")
reviewed_version <- c(presize = "0.3.11", pROC = "1.19.0.1", MKpower = "1.1",
                       epiR = "2.0.95", MKmisc = "2.0")
reviewed_date <- "2026-08-24"

installed_version <- vapply(compare_pkgs, function(p) {
  if (requireNamespace(p, quietly = TRUE)) {
    as.character(utils::packageVersion(p))
  } else {
    NA_character_
  }
}, character(1))

provenance <- data.frame(
  package = compare_pkgs,
  version_reviewed_for_this_matrix = unname(reviewed_version[compare_pkgs]),
  reviewed_on = reviewed_date,
  currently_installed_version = unname(installed_version),
  status = ifelse(
    is.na(installed_version),
    "NOT INSTALLED in the environment that generated this file -- capability entries above are based on the CRAN documentation reviewed on the date above, not re-verified here",
    ifelse(installed_version == unname(reviewed_version[compare_pkgs]),
           "installed; matches the version reviewed",
           "installed; DIFFERENT version than reviewed -- re-check the matrix before trusting it")
  ),
  stringsAsFactors = FALSE
)
write.csv(provenance, file.path(assets_dir, "table3_provenance.csv"), row.names = FALSE)
cat("Feature matrix (table3_feature_matrix.csv; Table 4 of v0.6.9, Table S3.3e of Additional file 3 in the R1): ", nrow(table3), " capability rows x ", ncol(table3) - 1,
    " comparator packages; provenance written for ", nrow(provenance), " packages\n", sep = "")

## ---------------------------------------------------------------------
cat("\nFiles written to", out_dir, ":\n")
cat(sort(list.files(out_dir)), sep = "\n")
cat("\nSUMMARY: ", length(anchor_log), " anchors, all PASS: ",
    paste(anchor_log, collapse = "; "), "\n", sep = "")
cat("\n=== make_manuscript_assets.R done ===\n")
