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
## a moderately informative Se prior, a fairly tight Sp prior, a
## prevalence prior centred at 0.20, and the half-width precision targets
## Se +/- 0.07 and Sp +/- 0.05 at 80% assurance. loss_rate = 0 keeps the
## reported N as "enrolled = analysed", so it is directly comparable
## across methods that do and do not model attrition.
PRIOR_SE   <- c(17, 3)
PRIOR_SP   <- c(18, 2)
PRIOR_PREV <- c(4, 16)                 # Beta(4, 16): mean 0.20
DELTA_SE   <- 0.07
DELTA_SP   <- 0.05
TARGET_ASSURANCE <- 0.80
LOSS_RATE  <- 0
B_MAIN     <- 20000
SEED       <- 2026

E_SE   <- PRIOR_SE[1]   / sum(PRIOR_SE)     # 0.85
E_SP   <- PRIOR_SP[1]   / sum(PRIOR_SP)     # 0.90
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
       ylab = "P(CI width <= target)", ylim = c(0, 1))
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
emit("Figure_1", draw_fig1, 7, 5)
cat(sprintf("Figure 1: Buderer n=%d -> assurance %.4f; n for 0.80 assurance = %s\n",
            n_bud_fig1, assur_fig1[ns == n_bud_fig1], n_80_fig1))

## ===========================================================================
## FIGURE 2 -- Required Se-variance multiplier for an imperfect reference
##             standard, at the estimand a study actually plans for
## ===========================================================================
## REPLACES the previous release's heat map of the Rogan-Gladen prevalence
## factor 1/(Se_ref + Sp_ref - 1)^2 relabelled as a sample-size multiplier.
## That relabelling was structurally invalid: as documented in
## ss_imperfect_ref()'s Details ("The retired VIF"), the multiplier that
## actually governs the corrected-Se sample size is not a function of
## (Se_ref, Sp_ref) alone -- it also depends on prevalence, Se and Sp -- so
## no function of (Se_ref, Sp_ref) alone can represent it, and no amount of
## relabelling the old heat map fixes that. This figure instead plots the
## quantity ss_imperfect_ref() itself reports for exactly this purpose,
## multiplier_se_corrected, over the same (Se_ref, Sp_ref) grid as before,
## with Se, Sp and prevalence held fixed at the article's harmonized
## planning values (E_SE = 0.85, E_SP = 0.90, E_PREV = 0.20, defined
## above). B = 0 and sensitivity_table = FALSE because multiplier_se_corrected
## is a closed-form delta-method ratio (see ss_imperfect_ref()'s Details) --
## there is no Monte Carlo replicate to validate here, and the function's
## own internal sensitivity_table uses a coarser grid built for a different
## purpose.
se_ref_grid <- seq(0.80, 0.99, by = 0.01)
sp_ref_grid <- seq(0.80, 0.99, by = 0.01)

## ss_imperfect_ref() enforces both Se_ref + Sp_ref > 1 and the min_youden
## guard (default 0.5, i.e. Se_ref + Sp_ref >= 1.5) as hard errors, not as
## NAs in a vectorised return. A grid point that fails either is caught
## here and recorded as excluded (NA in the matrix, left uncoloured by
## image() and skipped by contour()) instead of aborting the whole figure.
mult_se_corrected <- function(se_ref, sp_ref) {
  res <- tryCatch(
    ss_imperfect_ref(Se = E_SE, Sp = E_SP, prev = E_PREV,
                      d_se = DELTA_SE, d_sp = DELTA_SP,
                      Se_ref = se_ref, Sp_ref = sp_ref,
                      loss_rate = LOSS_RATE, B = 0,
                      sensitivity_table = FALSE, seed = SEED),
    error = function(e) NULL
  )
  if (is.null(res)) NA_real_ else res$multiplier_se_corrected
}
mult_fig2 <- outer(se_ref_grid, sp_ref_grid, Vectorize(mult_se_corrected))
n_excluded_fig2 <- sum(is.na(mult_fig2))
if (n_excluded_fig2 > 0) {
  cat(sprintf(
    "  [Figure 2] %d/%d grid points excluded (reference Youden index below the min_youden guard)\n",
    n_excluded_fig2, length(mult_fig2)))
}

## Anchor check: at Se_ref = 0.90, Sp_ref = 0.95 (Se = 0.85, Sp = 0.90,
## prev = 0.20) the required multiplier was independently verified to be
## 2.470, against 1.384 for the retired Rogan-Gladen factor at the same
## point. Refuse to draw the figure if this stops reproducing, rather than
## silently plotting a value nobody checked.
anchor_fig2 <- mult_se_corrected(0.90, 0.95)
if (!isTRUE(all.equal(anchor_fig2, 2.470, tolerance = 0.001))) {
  stop(sprintf(
    paste("Figure 2 anchor check failed: multiplier_se_corrected(Se_ref=0.90,",
          "Sp_ref=0.95 | Se=%.2f, Sp=%.2f, prev=%.2f) = %.6f, expected",
          "approximately 2.470. Refusing to draw a figure whose verified",
          "anchor point does not reproduce."),
    E_SE, E_SP, E_PREV, anchor_fig2
  ), call. = FALSE)
}

## Colour scale: LINEAR, not logarithmic. Over this grid the multiplier
## spans about an 8.6-fold range (observed range printed below) -- wide,
## but not the multiple-orders-of-magnitude spread that would make a
## linear ramp unreadable, and the intended reading ("this reference
## standard needs roughly this many times the classical N") is direct on
## a linear scale and requires translation on a log one. zlim is anchored
## at 1 (no inflation, i.e. a perfect reference, same convention as the
## retired figure) up to the grid's observed maximum, rounded up; unlike
## the retired figure, nothing is pmin()-capped, because the real range
## here does not need it to stay legible.
zlim_fig2 <- c(1, ceiling(max(mult_fig2, na.rm = TRUE)))

draw_fig2 <- function() {
  layout(matrix(c(1, 2), nrow = 1), widths = c(6, 1))
  cols <- hcl.colors(50, "YlOrRd", rev = TRUE)
  par(mar = c(4.3, 4.6, 1.2, 1))
  image(se_ref_grid, sp_ref_grid, mult_fig2, col = cols, zlim = zlim_fig2,
        xlab = "Reference sensitivity (Se_ref)",
        ylab = "Reference specificity (Sp_ref)")
  ## Levels 2/4/6/8 verified against mult_fig2 to fall strictly inside the
  ## observed range [1.22, 10.54] with a substantial area of the grid on
  ## both sides of each -- none of them touches an edge of the domain or
  ## of the data range.
  contour(se_ref_grid, sp_ref_grid, mult_fig2, add = TRUE,
          levels = c(2, 4, 6, 8), labcex = 0.7, col = "grey20")
  par(mar = c(4.3, 0.5, 1.2, 3.6))
  zseq <- seq(zlim_fig2[1], zlim_fig2[2], length.out = 50)
  image(1, zseq, matrix(zseq, nrow = 1), col = cols, zlim = zlim_fig2,
        axes = FALSE, xlab = "", ylab = "")
  axis(4, at = c(1, 2, 4, 6, 8, zlim_fig2[2]), las = 1)
  mtext("Se multiplier", side = 4, line = 2.3, cex = 0.9); box()
  layout(1)
}
emit("Figure_2", draw_fig2, 8, 5)
cat(sprintf(
  "Figure 2: multiplier_se_corrected(Se_ref=0.90,Sp_ref=0.95)=%.3f (Rogan-Gladen VIF=%.3f at same point) ; range over grid = [%.3f, %.3f]\n",
  anchor_fig2, 1 / (0.90 + 0.95 - 1)^2,
  min(mult_fig2, na.rm = TRUE), max(mult_fig2, na.rm = TRUE)))

## ===========================================================================
## FIGURE 3 -- Sample size for conclusive net benefit
## ===========================================================================
## Unchanged: required total N, by decision threshold, for the lower 95%
## confidence limit of net benefit to exceed both default strategies,
## under the prospective-cohort (design = "cohort", the default) variance.
nb_fig3 <- run_checked("Figure 3 (net benefit)",
                        ss_net_benefit(B = 4000, seed = SEED))
dnb_fig3 <- nb_fig3$N_by_pt

draw_fig3 <- function() {
  par(mar = c(4.3, 4.6, 1.2, 1))
  plot(dnb_fig3$pt, dnb_fig3$N_required, type = "b", pch = 19, lwd = 2,
       col = "#6a1b9a",
       xlab = "Threshold probability (p_t)", ylab = "Required total N")
  grid(col = "grey85")
}
emit("Figure_3", draw_fig3, 7, 5)
cat("Figure 3 N by threshold:", paste(dnb_fig3$N_required, collapse = ", "),
    " (conservative N =", nb_fig3$N_conservative, ")\n")

## ===========================================================================
## FIGURE 4 -- Required N as a NESTED sequence of added uncertainty sources
## ===========================================================================
## Five bars, each adding exactly one source of uncertainty on top of the
## previous bar, all else held fixed. This replaces the earlier "method
## comparison" framing (five unrelated methods, each with its own default
## assumptions) with a design in which each step is a strict superset of
## the one before it, so that the increase from bar to bar is attributable
## to the single feature that step turns on.
##
## Step 1 (Buderer, deterministic) has no Monte Carlo component. Steps 2-5
## all go through ss_unified(): step 2 turns on joint Se/Sp precision by
## simulating estimation uncertainty against a PERFECT reference standard
## (Se_ref = Sp_ref = 1, so the "apparent" accuracy in that call is the
## true accuracy); step 3 additionally requires the AUC half-width target;
## step 4 replaces the perfect reference with the imperfect one actually
## assumed elsewhere in the article (Se_ref = 0.90, Sp_ref = 0.95); step 5
## additionally requires a conclusive (CI-based) net benefit at every
## threshold in the default pt_range.
##
## IMPORTANT: ss_unified()'s default N_range (seq(200, 1000, by = 20)) is
## too narrow for steps 4-5 -- the search would exhaust the range, return
## max(N_range) without a real crossing, and silently look like a smaller
## (wrong) answer instead of failing loudly. Each call below is given an
## explicit N_range verified to bracket the true crossing comfortably
## inside its interior, and run_checked() aborts if any call nonetheless
## fails to converge.

## Step 1: classical Buderer, applied per-arm and then combined across the
## Se and Sp requirements via the common prevalence -- no simulation.
se_arm_n <- buderer_n(0.85, 0.07)
sp_arm_n <- buderer_n(0.90, 0.05)
N_step1  <- ceiling(max(se_arm_n / 0.20, sp_arm_n / 0.80))

## Step 2: + joint Se/Sp precision (perfect reference, no AUC, no NB).
u_step2 <- run_checked("Figure 4 step 2 (joint Se/Sp precision)",
  ss_unified(prior_se = PRIOR_SE, prior_sp = PRIOR_SP, prior_prev = PRIOR_PREV,
             Se_ref = 1, Sp_ref = 1, loss_rate = LOSS_RATE,
             delta_se = DELTA_SE, delta_sp = DELTA_SP, delta_auc = 0,
             check_nb = FALSE, target_assurance = TARGET_ASSURANCE,
             N_range = seq(750, 1000, by = 10), B = B_MAIN, seed = SEED))
N_step2 <- u_step2$n_total

## Step 3: + AUC gate (perfect reference still; delta_auc turned on).
u_step3 <- run_checked("Figure 4 step 3 (+ AUC gate)",
  ss_unified(prior_se = PRIOR_SE, prior_sp = PRIOR_SP, prior_prev = PRIOR_PREV,
             Se_ref = 1, Sp_ref = 1, loss_rate = LOSS_RATE,
             delta_se = DELTA_SE, delta_sp = DELTA_SP, delta_auc = 0.06,
             check_nb = FALSE, target_assurance = TARGET_ASSURANCE,
             N_range = seq(750, 1000, by = 10), B = B_MAIN, seed = SEED))
N_step3 <- u_step3$n_total

## Step 4: + imperfect reference standard (Se_ref/Sp_ref no longer perfect).
u_step4 <- run_checked("Figure 4 step 4 (+ imperfect reference)",
  ss_unified(prior_se = PRIOR_SE, prior_sp = PRIOR_SP, prior_prev = PRIOR_PREV,
             Se_ref = 0.90, Sp_ref = 0.95, loss_rate = LOSS_RATE,
             delta_se = DELTA_SE, delta_sp = DELTA_SP, delta_auc = 0.06,
             check_nb = FALSE, target_assurance = TARGET_ASSURANCE,
             N_range = seq(1050, 1400, by = 10), B = B_MAIN, seed = SEED))
N_step4 <- u_step4$n_total

## Step 5: + net benefit (inference-based, both default comparisons).
u_step5 <- run_checked("Figure 4 step 5 (+ net benefit)",
  ss_unified(prior_se = PRIOR_SE, prior_sp = PRIOR_SP, prior_prev = PRIOR_PREV,
             Se_ref = 0.90, Sp_ref = 0.95, loss_rate = LOSS_RATE,
             delta_se = DELTA_SE, delta_sp = DELTA_SP, delta_auc = 0.06,
             check_nb = TRUE, target_assurance = TARGET_ASSURANCE,
             N_range = seq(1950, 2600, by = 10), B = B_MAIN, seed = SEED))
N_step5 <- u_step5$n_total

fig4_vals <- c(N_step1, N_step2, N_step3, N_step4, N_step5)
fig4_labels <- c(
  "1. Buderer\n(deterministic)",
  "2. + joint Se/Sp\nprecision",
  "3. + AUC\ngate",
  "4. + imperfect\nreference",
  "5. + net\nbenefit"
)
fig4_deltas <- diff(fig4_vals)

draw_fig4 <- function() {
  ## Sequential, single-hue ramp: distinguishable by lightness alone, so it
  ## reads correctly in grayscale print and under the common forms of
  ## color-vision deficiency.
  pal <- grDevices::colorRampPalette(c("#deebf7", "#08306b"))(5)
  par(mar = c(6.4, 4.8, 1.4, 1))
  bp <- barplot(fig4_vals, col = pal, border = "grey20",
                names.arg = fig4_labels, las = 1, cex.names = 0.72,
                ylab = "Total required sample size (N)",
                ylim = c(0, max(fig4_vals) * 1.20))
  text(bp, fig4_vals, labels = fig4_vals, pos = 3, cex = 0.88, font = 2)
  ## Delta labels between consecutive bars, so the flat step (2 -> 3) is
  ## shown exactly as it is: the AUC gate adds nothing at this operating
  ## point, and that is a real result, not an omission.
  for (i in seq_along(fig4_deltas)) {
    xmid <- mean(bp[i:(i + 1)])
    ymid <- max(fig4_vals[i], fig4_vals[i + 1]) + 0.06 * max(fig4_vals)
    is_flat <- fig4_deltas[i] == 0
    text(xmid, ymid, sprintf("%+d", fig4_deltas[i]), cex = 0.76,
         col = if (is_flat) "grey45" else "grey15",
         font = if (is_flat) 3 else 1)
  }
}
emit("Figure_4", draw_fig4, 8, 5.5)
cat("Figure 4 nested-sequence values:",
    paste(sprintf("step%d=%d", seq_along(fig4_vals), fig4_vals), collapse = ", "), "\n")
cat("Figure 4 deltas between consecutive steps:",
    paste(sprintf("%+d", fig4_deltas), collapse = ", "), "\n")

## ===========================================================================
## TABLE -- harmonized comparison of the five planning approaches
## ===========================================================================
## The five methods evaluated under the same common assumptions (Methods),
## each contributing the N it actually produces under those assumptions --
## not a value read off Figure 4, since three of the five methods (BAM,
## the joint Se/Sp+AUC gate, and the imperfect-reference correction) are
## not part of the nested Figure 4 sequence at all.

## Method 1: Buderer (classical, deterministic) -- identical to Figure 4
## step 1, reused here rather than recomputed.
h_buderer <- list(N = N_step1, assurance = NA_real_,
                   assurance_type = "deterministic (normal-approximation target width; not a simulated assurance)")

## Method 2: BAM, exact joint Beta-Binomial search, full-width deltas.
bam_res <- run_checked("Harmonized table: BAM (exact)",
  bam_sample_size(prior_se = PRIOR_SE, prior_sp = PRIOR_SP, prior_prev = PRIOR_PREV,
                   delta_se = 2 * DELTA_SE, delta_sp = 2 * DELTA_SP,
                   target_assurance = TARGET_ASSURANCE, method = "exact",
                   B = B_MAIN, seed = SEED))
h_bam <- list(N = bam_res$N_total, assurance = bam_res$joint_assurance,
               assurance_type = "exact Bayesian (Beta-Binomial), closed form, no Monte Carlo error")

## Method 3: joint Se + Sp precision with a deterministic AUC gate,
## prospective-cohort variance for the Se/Sp Monte Carlo (package defaults
## for AUC = 0.90 and delta_auc = 0.05, which are not part of the five
## common parameters above because this method, unlike ss_unified(), has
## no prior on Se/Sp/prevalence to harmonize).
joint_res <- run_checked("Harmonized table: joint Se/Sp + AUC",
  joint_sample_size(Se = E_SE, Sp = E_SP, delta_se = DELTA_SE, delta_sp = DELTA_SP,
                     prev = E_PREV, design = "cohort",
                     target_prob = TARGET_ASSURANCE, B = B_MAIN, seed = SEED))
h_joint <- list(N = joint_res$n_total, assurance = joint_res$joint_prob_se_sp,
                 assurance_type = "Monte Carlo (cohort design) for Se & Sp jointly; AUC via a deterministic Hanley-McNeil gate")

## Method 4: imperfect reference, apparent-sensitivity estimand --
## deterministic given Se, Sp, prev, Se_ref, Sp_ref; B = 0 skips the
## optional Monte Carlo bias check (see ss_imperfect_ref()'s
## mc_validation), which is not needed for the sample size itself.
imp_res <- ss_imperfect_ref(Se = E_SE, Sp = E_SP, d_se = DELTA_SE, d_sp = DELTA_SP,
                             prev = E_PREV, Se_ref = 0.90, Sp_ref = 0.95,
                             loss_rate = LOSS_RATE, B = 0, seed = SEED)
h_imperfect <- list(N = imp_res$n_total, assurance = NA_real_,
                      assurance_type = paste(
                        "deterministic closed form for the apparent",
                        "(reference-biased) sensitivity; no variance",
                        "inflation factor is applied; not simulated"))

## Method 5: unified framework with every source of uncertainty active --
## identical to Figure 4 step 5, reused here rather than recomputed.
h_unified <- list(N = N_step5, assurance = u_step5$joint_assurance,
                    assurance_type = "Monte Carlo, joint over Se, Sp, AUC and net benefit simultaneously; lower-bound decision rule")

## joint_sample_size() does not echo delta_auc back in its return value;
## the call above left it at the function's own default (0.05), so that
## default is recorded here as a literal, not re-derived at runtime.
JOINT_DELTA_AUC_DEFAULT <- 0.05

harmonized <- data.frame(
  method = c("Buderer (classical)", "BAM (exact mode)",
             "Joint Se/Sp + AUC", "Imperfect reference", "Unified (all sources active)"),
  N = c(h_buderer$N, h_bam$N, h_joint$N, h_imperfect$N, h_unified$N),
  assurance_achieved = c(h_buderer$assurance, h_bam$assurance, h_joint$assurance,
                          h_imperfect$assurance, h_unified$assurance),
  assurance_type = c(h_buderer$assurance_type, h_bam$assurance_type, h_joint$assurance_type,
                      h_imperfect$assurance_type, h_unified$assurance_type),
  Se = c(0.85, E_SE, E_SE, E_SE, NA),
  Sp = c(0.90, E_SP, E_SP, E_SP, NA),
  prevalence = c(0.20, E_PREV, E_PREV, E_PREV, NA),
  prior_se = c(NA, "17,3", NA, NA, "17,3"),
  prior_sp = c(NA, "18,2", NA, NA, "18,2"),
  prior_prev = c(NA, "4,16", NA, NA, "4,16"),
  delta_se = c(DELTA_SE, 2 * DELTA_SE, DELTA_SE, DELTA_SE, DELTA_SE),
  delta_sp = c(DELTA_SP, 2 * DELTA_SP, DELTA_SP, DELTA_SP, DELTA_SP),
  AUC = c(NA, NA, joint_res$AUC, NA, NA),
  delta_auc = c(NA, NA, JOINT_DELTA_AUC_DEFAULT, NA, 0.06),
  Se_ref = c(NA, NA, NA, 0.90, 0.90),
  Sp_ref = c(NA, NA, NA, 0.95, 0.95),
  loss_rate = c(NA, NA, NA, LOSS_RATE, LOSS_RATE),
  check_net_benefit = c(NA, NA, NA, NA, TRUE),
  B = c(NA, B_MAIN, B_MAIN, 0, B_MAIN),
  seed = c(NA, SEED, SEED, SEED, SEED),
  notes = c(
    "N = ceiling(max(buderer_n(0.85,0.07)/0.20, buderer_n(0.90,0.05)/0.80))",
    "delta_se/delta_sp given as full CI widths (2x the half-width used elsewhere); prior_prev harmonized to Beta(4,16) (package default is Beta(6,14))",
    "AUC and delta_auc left at joint_sample_size()'s own defaults (0.90 and 0.05); not part of the five harmonized parameters, since this method has no Se/Sp/prevalence prior to harmonize",
    "Se_ref/Sp_ref harmonized to the values used in Figure 4 steps 4-5 (package defaults are 0.92/0.95)",
    "identical configuration to Figure 4 step 5"
  ),
  stringsAsFactors = FALSE
)
write.csv(harmonized, file.path(assets_dir, "table_harmonized_comparison.csv"), row.names = FALSE)
cat("\nHarmonized comparison table:\n")
print(harmonized[, c("method", "N", "assurance_achieved")], row.names = FALSE)

## ===========================================================================
## TABLE 5 -- feasibility ceiling of the net-benefit criterion, by
##            decision-threshold range
## ===========================================================================
## ss_unified()'s nb_ceiling is the largest joint assurance the check_nb
## criterion could ever reach as N -> Inf (see ?ss_unified, "The CI-based
## net-benefit criterion has a hard ceiling"). It is computed BEFORE the
## grid search, from the priors, Se_ref, Sp_ref, pt_range and seed alone --
## not from N_range or from the search's own B -- so the cheapest possible
## call that still returns it uses a single-value N_range and a small B.
## Whether that one-point search itself "converges" is irrelevant here
## (only nb_ceiling is read off the return value), so its warnings --
## either "criterion unreachable at this target_assurance" for a ceiling
## below TARGET_ASSURANCE, or "no N reached target assurance" for the
## single grid point -- are expected for several of the ranges below and
## are deliberately suppressed rather than treated as failures.
table5_pt_ranges <- list(
  c(0.15, 0.40),   # used throughout the article (Figure 4 step 5, Table 4)
  c(0.20, 0.30),
  c(0.15, 0.25),
  c(0.10, 0.40),
  c(0.10, 0.30),
  c(0.05, 0.50)
)

nb_ceiling_for <- function(pt_range) {
  suppressWarnings(ss_unified(
    prior_se = PRIOR_SE, prior_sp = PRIOR_SP, prior_prev = PRIOR_PREV,
    Se_ref = 0.90, Sp_ref = 0.95, loss_rate = LOSS_RATE,
    delta_se = DELTA_SE, delta_sp = DELTA_SP, delta_auc = 0.06,
    check_nb = TRUE, pt_range = pt_range,
    target_assurance = TARGET_ASSURANCE,
    N_range = 500, B = 100, seed = SEED
  ))$nb_ceiling
}
table5_ceiling <- vapply(table5_pt_ranges, nb_ceiling_for, numeric(1))

## Anchor check: the range actually used throughout the article (0.15-0.40)
## is checked against 0.8796, NOT against a value read off this package's
## own Monte Carlo output. 0.8796 was derived independently of
## nb_assurance_ceiling() entirely, by tensor-product Gauss-Legendre
## quadrature over the (prev, Se, Sp) priors (n = 300 and n = 600 nodes
## per dimension agree to six decimals: 0.879606), and cross-checked
## against a 5e6-draw Monte Carlo run. Anchoring against a number the
## package computed itself would only confirm that the code reproduces
## its own prior output, not that the output is correct -- exactly the
## failure this check exists to rule out. The tolerance (2e-4) is set
## against nb_assurance_ceiling()'s own documented Monte Carlo error
## (below 1e-4 at its default B_ceiling; see ?ss_unified's internal
## nb_assurance_ceiling), with headroom for ordinary run-to-run varia-
## tion. Refuse to write the table if this stops reproducing, rather
## than silently publishing a value nobody checked.
if (!isTRUE(all.equal(table5_ceiling[1], 0.8796, tolerance = 2e-4))) {
  stop(sprintf(
    paste("Table 5 anchor check failed: nb_ceiling(pt_range = [0.15, 0.40])",
          "= %.6f, expected approximately 0.8796 (independently derived by",
          "Gauss-Legendre quadrature, not by this package). Refusing to",
          "write a table whose verified anchor point does not reproduce."),
    table5_ceiling[1]
  ), call. = FALSE)
}

table5 <- data.frame(
  pt_range = vapply(table5_pt_ranges,
                     function(pr) sprintf("%.2f-%.2f", pr[1], pr[2]),
                     character(1)),
  ceiling = table5_ceiling,
  reaches_target_0.80 = ifelse(table5_ceiling >= TARGET_ASSURANCE, "Yes", "No"),
  stringsAsFactors = FALSE
)
write.csv(table5, file.path(assets_dir, "table5_nb_ceiling.csv"), row.names = FALSE)
cat("\nTable 5 (net-benefit feasibility ceiling by threshold range):\n")
print(table5, row.names = FALSE)

## ===========================================================================
## TABLE 1 -- functions exported by dtasamplesize, derived from NAMESPACE
## ===========================================================================
## Read the actual NAMESPACE rather than maintaining a hand-written list,
## so the table cannot silently drift out of sync with the package's
## exports. The one-line description of each function is its Rd \title{},
## read from the installed help database -- again, not retyped by hand.
ns_lines <- readLines(file.path(pkg_root, "NAMESPACE"), warn = FALSE)
export_lines <- grep("^export\\(", ns_lines, value = TRUE)
fn_names <- sort(sub("^export\\(([^)]+)\\)$", "\\1", export_lines))

rd_title <- function(fn) {
  rd_path <- file.path(pkg_root, "man", paste0(fn, ".Rd"))
  if (!file.exists(rd_path)) return(NA_character_)
  rd_lines <- readLines(rd_path, warn = FALSE)
  title_line <- grep("\\\\title\\{", rd_lines, value = TRUE)[1]
  if (is.na(title_line)) return(NA_character_)
  sub(".*\\\\title\\{(.*)\\}.*", "\\1", title_line)
}
purposes <- vapply(fn_names, rd_title, character(1))
if (anyNA(purposes)) {
  stop("Table 1: no \\title{} found for: ",
       paste(fn_names[is.na(purposes)], collapse = ", "),
       ". Check that man/ is in sync with NAMESPACE.", call. = FALSE)
}

table1 <- data.frame(Function = fn_names, Purpose = unname(purposes),
                      stringsAsFactors = FALSE)
write.csv(table1, file.path(assets_dir, "table1_functions.csv"), row.names = FALSE)
cat("\nTable 1: ", nrow(table1), " exported functions (from NAMESPACE)\n", sep = "")

## ===========================================================================
## TABLE 3 -- feature comparison with existing R packages
## ===========================================================================
## This matrix is curated by hand -- it records a judgment about what each
## comparator package's public functions actually do, not something that
## can be derived mechanically -- and lives as a separate, versioned CSV
## (table3_feature_matrix_source.csv) precisely so that the curation is
## visible and diffable on its own, independent of this script. The block
## below only re-emits it, together with a provenance note.
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
