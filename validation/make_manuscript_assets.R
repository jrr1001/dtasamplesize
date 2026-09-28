# Generate the manuscript figures and tables for dtasamplesize
# ---------------------------------------------------------------------------
# This script regenerates every figure and every data-derived table used in
# the accompanying article, from the installed package and nothing else: no
# number here is copied by hand from a previous run. Figures are written as
# a vector PDF (BMC's preferred format for line art) and a 600 dpi TIFF, to
# submission-BMC-MRM/figures/. Tables are written as CSV to
# R-package/manuscript_assets/. Re-running this script regenerates all of it
# byte-for-byte, because every quantity is either a closed-form calculation
# or a Monte Carlo simulation under a fixed seed.
#
# SCOPE (article now reduced to two figures and four tables). The article
# describes only bam_sample_size(), mc_validate_buderer() and
# joint_sample_size(); ss_unified(), the AUC gate inside it, ss_net_benefit(),
# ss_imperfect_ref(), ss_adaptive_prevalence() and ss_time_dependent_roc()
# remain in the package but are no longer described in this article, and this
# script no longer builds figures or tables for them. This script writes
# Figure 1 (unchanged), Figure 2 (the joint-assurance-vs-N curve, new in this
# revision), Table 1 (the three estimators described in the article), Table 3
# (the feature matrix, read from table3_feature_matrix_source.csv) and
# Table 4 (the three surviving methods compared under common assumptions).
# Table 2 (cross-validation of core formulas) is written separately by
# cross_validation.R, not by this script -- see validation/README.md.
#
# Run after installing the package, e.g. from the package root:
#   Rscript validation/make_manuscript_assets.R
# (the script locates the repository layout from its own path, so it does
# not depend on the current working directory).

suppressMessages(library(dtasamplesize))
options(width = 120, dtasamplesize.warn_small_B = FALSE)

## ---------------------------------------------------------------------
## Locate the repository layout from this script's own path, so the
## script works regardless of the caller's working directory. Falls back
## to a couple of common invocation patterns (running from the package
## root, or from inside validation/) if --file= is not available (e.g.
## when the script is sourced interactively).
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
    "  Rscript validation/make_manuscript_assets.R\n",
    "from the dtasamplesize package root, or from inside validation/.",
    call. = FALSE
  )
}

script_dir      <- locate_script_dir()                 # .../dtasamplesize/validation
pkg_root        <- dirname(script_dir)                 # .../dtasamplesize
r_package_root  <- dirname(pkg_root)                    # .../R-package
repo_root       <- dirname(r_package_root)               # repository root

figures_dir <- file.path(repo_root, "submission-BMC-MRM", "figures")
assets_dir  <- file.path(r_package_root, "manuscript_assets")
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(assets_dir, showWarnings = FALSE, recursive = TRUE)

cat("Package version:", as.character(utils::packageVersion("dtasamplesize")), "\n")
cat("Figures ->", figures_dir, "\n")
cat("Tables   ->", assets_dir, "\n\n")

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
## (Beta(18, 2), mean 0.90) used only in Table 4, where it is compared
## against the classical and Monte Carlo methods under matched point
## assumptions (Se = 0.85, Sp = 0.90, prevalence = 0.20); it is NOT the
## package's own default Sp prior. Figure 2 instead uses the package's
## actual default, vague Sp prior, Beta(2, 2) -- see PRIOR_SP_VAGUE below
## and the article's Results/Limitations for why the two differ (N = 678
## under the vague default vs. N = 672 under the harmonized informative
## prior).
PRIOR_SE   <- c(17, 3)
PRIOR_SP   <- c(18, 2)                 # harmonized (informative), Table 4 only
PRIOR_PREV <- c(4, 16)                 # Beta(4, 16): mean 0.20
DELTA_SE   <- 0.07                     # half-width, as used by the classical
DELTA_SP   <- 0.05                     # and Monte Carlo methods in Table 4
TARGET_ASSURANCE <- 0.80
B_MAIN     <- 20000
SEED       <- 2026

E_SE   <- PRIOR_SE[1]   / sum(PRIOR_SE)     # 0.85
E_SP   <- PRIOR_SP[1]   / sum(PRIOR_SP)     # 0.90 (harmonized prior mean)
E_PREV <- PRIOR_PREV[1] / sum(PRIOR_PREV)   # 0.20

## ===========================================================================
## FIGURE 1 -- Assurance of the classical Buderer sample size
## ===========================================================================
## Unchanged from the previous release: the probability that the realised
## 95% CI width for Se is within the target, as a function of the number
## of diseased subjects, showing that the Buderer n itself clears the
## target width in only a slim majority of studies.
ns <- seq(40, 320, by = 10)
assur_fig1 <- sapply(ns, function(n) {
  mc_validate_buderer(Se = 0.85, d = 0.07, n_diseased = n, B = 4000,
                      ci_method = "wald", seed = SEED)$results$P_width_target[1]
})
n_bud_fig1 <- buderer_n(0.85, 0.07)
n_80_fig1  <- ns[which(assur_fig1 >= 0.80)[1]]

draw_fig1 <- function() {
  par(mar = c(4.3, 4.6, 1.2, 1))
  plot(ns, assur_fig1, type = "l", lwd = 2, col = "#1f4e79",
       xlab = "Number of diseased subjects (n)",
       ylab = "P(CI width <= target)", ylim = c(0, 1), las = 1)
  abline(h = 0.80, lty = 3, col = "grey40")
  abline(v = n_bud_fig1, lty = 2, col = "#c00000")
  points(n_bud_fig1, assur_fig1[ns == n_bud_fig1], pch = 19, col = "#c00000")
  text(n_bud_fig1 - 4, 0.36,
       sprintf("Buderer n=%d\nassurance=%.2f", n_bud_fig1, assur_fig1[ns == n_bud_fig1]),
       col = "#c00000", cex = 0.78, pos = 2)
  if (!is.na(n_80_fig1)) {
    abline(v = n_80_fig1, lty = 2, col = "#2e7d32")
    text(n_80_fig1 + 6, 0.72, sprintf("n=%d for\n0.80 assurance", n_80_fig1),
         col = "#2e7d32", cex = 0.78, pos = 4)
  }
  legend("bottomright", bty = "n", cex = 0.8,
         legend = c("Assurance curve", "0.80 target", "Buderer n", "n for 0.80"),
         col = c("#1f4e79", "grey40", "#c00000", "#2e7d32"),
         lty = c(1, 3, 2, 2), lwd = c(2, 1, 1, 1))
}
emit("Figure_1", draw_fig1, FIG_WIDTH_IN, FIG_HEIGHT_IN)
cat(sprintf("Figure 1: Buderer n=%d -> assurance %.4f; n for 0.80 assurance = %s\n",
            n_bud_fig1, assur_fig1[ns == n_bud_fig1], n_80_fig1))

## ===========================================================================
## FIGURE 2 -- Joint assurance for sensitivity and specificity as a
##             function of total sample size N (exact Beta-Binomial)
## ===========================================================================
## NEW in this revision. REPLACES the previous release's imperfect-reference
## multiplier heat map (Figure 2), the net-benefit-by-threshold plot
## (Figure 3) and the nested-sequence bar chart (Figure 4) -- none of ss_
## unified(), ss_net_benefit(), ss_imperfect_ref() and the AUC gate is
## described by the article any more (see the SCOPE note at the top of this
## file). This is the single figure the article's Results and Figure 2
## legend refer to: bam_sample_size(method = "exact")'s joint assurance for
## Se and Sp, as a function of the total sample size N, under the article's
## worked-example priors -- Se ~ Beta(17, 3) and Sp ~ Beta(2, 2) (vague; the
## package defaults; see PRIOR_SP_VAGUE below, deliberately NOT the
## harmonized informative Sp prior used in Table 4), with prevalence ~
## Beta(4, 16) (NOT the package default prior_prev, which is c(6, 14)) --
## and full-width targets delta_se = 0.14, delta_sp = 0.10 (2x the half-widths used
## elsewhere in this script for the classical/Monte Carlo comparators),
## crossing the 0.80 target at N = 678 (article Results: "the required
## total sample size is N = 678, at which the joint assurance is
## 0.8003489948 ... at N = 677 the assurance is 0.7996848824").
##
## The curve is built by calling the SAME internal exact-calculation
## helpers that bam_sample_size(method = "exact") itself calls
## (.bam_exact_width_prob() / .bam_exact_joint_assurance(); see
## ?bam_sample_size, @details, for the closed-form Beta-Binomial derivation
## these implement) rather than calling the exported bam_sample_size()
## once per candidate N: that search function stops at the FIRST N whose
## assurance reaches target_assurance, so it cannot supply the points
## beyond the crossing needed to draw a curve. Calling the internal helpers
## directly evaluates the identical formula at every N, with no early stop,
## and builds the O(N^2) per-arm qbeta cache once for the whole curve
## rather than once per point -- the same caching bam_sample_size() itself
## relies on to stay fast.
PRIOR_SP_VAGUE <- c(2, 2)   # bam_sample_size()'s own default prior_sp
ALPHA_CI_FIG2  <- 0.95      # bam_sample_size()'s own default alpha_ci

N_grid_fig2 <- seq(350, 850, by = 1)
N_max_fig2  <- max(N_grid_fig2)
ci_lower_q_fig2 <- (1 - ALPHA_CI_FIG2) / 2
ci_upper_q_fig2 <- 1 - ci_lower_q_fig2

P_se_fig2 <- dtasamplesize:::.bam_exact_width_prob(
  N_max_fig2, PRIOR_SE[1], PRIOR_SE[2], 2 * DELTA_SE,
  ci_lower_q_fig2, ci_upper_q_fig2)
P_sp_fig2 <- dtasamplesize:::.bam_exact_width_prob(
  N_max_fig2, PRIOR_SP_VAGUE[1], PRIOR_SP_VAGUE[2], 2 * DELTA_SP,
  ci_lower_q_fig2, ci_upper_q_fig2)
## Degenerate-arm convention (n_d = 0 or n_nd = 0 always counts as a
## failure): the same convention bam_sample_size() enforces internally,
## both under method = "exact" and method = "monte_carlo" -- see
## ?bam_sample_size, @details.
P_se_fig2[1] <- 0
P_sp_fig2[1] <- 0

assur_fig2 <- vapply(N_grid_fig2, function(N) {
  dtasamplesize:::.bam_exact_joint_assurance(
    N, PRIOR_PREV[1], PRIOR_PREV[2], P_se_fig2, P_sp_fig2)
}, numeric(1))

n_cross_fig2 <- N_grid_fig2[which(assur_fig2 >= TARGET_ASSURANCE)[1]]
assur_at_cross_fig2 <- assur_fig2[N_grid_fig2 == n_cross_fig2]
assur_at_prev_fig2  <- assur_fig2[N_grid_fig2 == (n_cross_fig2 - 1L)]

## Anchor check, same pattern as the rest of this script (see Figure 1's
## companion check in reproduce_manuscript.R): refuse to draw the figure if
## the published crossing does not reproduce, rather than silently plotting
## a curve nobody checked. The comparison values carry full precision
## because they are computed exactly, with no Monte Carlo error to round
## away.
if (!isTRUE(all.equal(n_cross_fig2, 678L)) ||
    !isTRUE(all.equal(assur_at_cross_fig2, 0.8003489948, tolerance = 1e-8)) ||
    !isTRUE(all.equal(assur_at_prev_fig2, 0.7996848824, tolerance = 1e-8))) {
  stop(sprintf(
    paste("Figure 2 anchor check failed: crossing at N = %s (assurance",
          "%.10f), N-1 assurance %.10f; expected crossing at N = 678,",
          "assurance 0.8003489948 at N = 678 and 0.7996848824 at N = 677.",
          "Refusing to draw a figure whose verified crossing does not",
          "reproduce."),
    n_cross_fig2, assur_at_cross_fig2, assur_at_prev_fig2
  ), call. = FALSE)
}

draw_fig2 <- function() {
  par(mar = c(4.3, 4.6, 1.2, 1))
  plot(N_grid_fig2, assur_fig2, type = "l", lwd = 2, col = "#1f4e79",
       xlab = "Total sample size (N)",
       ylab = "Joint assurance for Se and Sp", ylim = c(0, 1), las = 1)
  abline(h = TARGET_ASSURANCE, lty = 3, col = "grey40")
  abline(v = n_cross_fig2, lty = 2, col = "#c00000")
  points(n_cross_fig2, assur_at_cross_fig2, pch = 19, col = "#c00000")
  text(n_cross_fig2 + 14, 0.28,
       sprintf("N=%d\nassurance=%.4f", n_cross_fig2, assur_at_cross_fig2),
       col = "#c00000", cex = 0.78, pos = 4)
  legend("bottomright", bty = "n", cex = 0.8,
         legend = c("Exact joint assurance", "0.80 target",
                    sprintf("N=%d (selected)", n_cross_fig2)),
         col = c("#1f4e79", "grey40", "#c00000"),
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
## Figures 3-4 above (see the SCOPE note at the top of this file). Each
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
bam_res <- run_checked("Table 4: BAM (exact, harmonized priors)",
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
joint_res <- run_checked("Table 4: joint Se/Sp (joint_sample_size)",
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
cat("\nTable 4 (method comparison):\n")
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
  stop("Table 3 source has ", ncol(table3), " columns, expected ",
       length(expected_cols), " (", paste(expected_cols, collapse = ", "),
       "). This usually means a data row in ", table3_source_path,
       " has an unescaped comma inside an unquoted field, which shifts",
       " read.csv()'s column count.", call. = FALSE)
}
if (!identical(colnames(table3)[1], "Capability")) {
  stop("Table 3 source: first column is named '", colnames(table3)[1],
       "', not 'Capability'. read.csv() likely promoted the Capability",
       " column to row names because a data row in ", table3_source_path,
       " parsed to more fields than the header -- check for an unquoted",
       " comma inside a cell.", call. = FALSE)
}
if (!identical(rownames(table3), as.character(seq_len(nrow(table3))))) {
  stop("Table 3 source: row.names(table3) are not the default 1..",
       nrow(table3), " -- read.csv() promoted a data column to row names,",
       " which means some row in ", table3_source_path, " has more fields",
       " than the header. Check every row for an unescaped comma inside",
       " an unquoted field and quote it.", call. = FALSE)
}
if (any(is.na(table3)) || any(table3 == "")) {
  bad <- which(is.na(table3) | table3 == "", arr.ind = TRUE)
  stop("Table 3 source has ", nrow(bad), " empty/NA cell(s) after parsing ",
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
cat("Table 3: ", nrow(table3), " capability rows x ", ncol(table3) - 1,
    " comparator packages; provenance written for ", nrow(provenance), " packages\n", sep = "")

## ---------------------------------------------------------------------
cat("\nFiles written to", figures_dir, ":\n")
cat(sort(list.files(figures_dir, pattern = "^Figure_")), sep = "\n")
cat("\nFiles written to", assets_dir, ":\n")
cat(sort(list.files(assets_dir)), sep = "\n")
cat("\n=== make_manuscript_assets.R done ===\n")
