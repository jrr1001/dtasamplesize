# Reproduce every headline number reported in the dtasamplesize manuscript
# ------------------------------------------------------------------------
# This script regenerates, from scratch, each quantity quoted in the article,
# using the same parameters and the same random seed. Every number the article
# reports is produced here; nothing is read from a stored dataset, because the
# study analyses no empirical data: all data are simulated inside the package's
# functions from the stated parameters and are fully determined by the seed.
#
# WHAT THIS SCRIPT DOES NOT DO. Every check below is a REPRODUCIBILITY
# check: it confirms that the installed package, run today, prints the same
# numbers the article prints. It is NOT, and cannot be, a VALIDITY check.
# A defect that biases a result reproduces just as cleanly as a result that
# does not -- a PASS below only means the code and the article agree with
# each other, not that either is statistically correct. This distinction is
# not academic: the numbers this script currently checks against were
# themselves revised twice (0.3.0, then 0.5.0 and 0.6.0) after earlier
# releases passed their own reproducibility checks against defective
# numbers. Validity is established separately, release by release, in
# NEWS.md and by independent review of the methods -- never by this script
# agreeing with itself.
#
# Run after installing the package:  R -f validation/reproduce_manuscript.R
#
# Optional: set DUMP_DATA <- TRUE below to also write the underlying simulated
# samples to CSV files in validation/simulated_data/ (see "Simulated data" at
# the end). The CSVs are a convenience: they are regenerated exactly by this
# script, so the script itself is the authoritative record of the data.

library(dtasamplesize)
options(width = 110, digits = 6)

# The article reports the numbers produced by version 0.6.0. Earlier releases
# contain the validity errors that 0.3.0, 0.5.0 and 0.6.0 correct in turn
# (see NEWS.md) -- notably a net-benefit sample size of 110 instead of 240,
# a BAM search that reported the marginal per-arm assurance instead of the
# joint one, a Figure 4 comparison confounded with three different
# prevalences across methods, an adaptive-prevalence design that discarded
# its own stage-1 pilot and understated recruitment by 21-36%, an
# imperfect-reference sample size inflated by a formula that belongs to a
# different estimand (prevalence, not the index test's Se/Sp), and a
# ss_unified() stopping rule whose selected N reached the declared
# assurance only 84% of the time across seeds -- so earlier installations
# return different values. Stop early rather than let a stale installation
# look like a discrepancy in the article.
if (utils::packageVersion("dtasamplesize") < "0.6.0") {
  stop("dtasamplesize ", utils::packageVersion("dtasamplesize"), " is installed, ",
       "but this script reproduces the numbers of version 0.6.0 or later.\n",
       "  Install the current source first, e.g.\n",
       "    remotes::install_github('jrr1001/dtasamplesize')\n",
       "  or, from a local copy of the repository:\n",
       "    install.packages('.', repos = NULL, type = 'source')",
       call. = FALSE)
}

DUMP_DATA <- FALSE          # set to TRUE to also export the simulated samples
SEED      <- 2026           # the seed cited throughout the manuscript

ok  <- function(x) if (isTRUE(x)) "PASS" else "**CHECK**"
sep <- function(t) cat("\n==========", t, "==========\n")
near <- function(a, b, tol) isTRUE(abs(a - b) <= tol)

## ---- 1. The motivating result (Results; Figure 1) ----
sep("1. Assurance gap of the classical Buderer sample size")
n_bud <- buderer_n(0.85, 0.07)
cat(sprintf("buderer_n(Se=0.85, d=0.07) = %d   (article: 100)   %s\n",
            n_bud, ok(n_bud == 100)))
cat(sprintf("buderer_n(Se=0.90, d=0.05) = %d   (article: 139)   %s\n",
            buderer_n(0.90, 0.05), ok(buderer_n(0.90, 0.05) == 139)))

mv <- mc_validate_buderer(Se = 0.85, d = 0.07, n_diseased = 100,
                          B = 4000, ci_method = "wald", seed = SEED)$results
cat(sprintf("P(CI width <= target) = %.4f   (article: 0.5730)   %s\n",
            mv$P_width_target[1], ok(near(mv$P_width_target[1], 0.5730, 5e-4))))
cat(sprintf("Wald coverage         = %.4f   (article: 0.9380)   %s\n",
            mv$coverage[1], ok(near(mv$coverage[1], 0.9380, 5e-4))))

## ---- 2. Net benefit under a prospective cohort (Results; Figure 3) ----
sep("2. Net-benefit sample size, cohort design")
nb <- suppressWarnings(ss_net_benefit(B = 4000, seed = SEED))
cat(sprintf("Conservative required N = %d   (article: 240)   %s\n",
            nb$N_conservative, ok(nb$N_conservative == 240)))
cat("N by threshold p_t:", paste(nb$N_by_pt$N_required, collapse = ", "), "\n")
cat("  (article, p_t = 0.10 to 0.50: 140, 60, 50, 60, 70, 80, 110, 160, 240)\n")

## ---- 3. Why the cohort variance matters (Results) ----
sep("3. Fixed-margin vs cohort sampling standard deviation")
Se <- 0.85; Sp <- 0.90; prev <- 0.20; pt <- 0.20; N <- 500
w  <- pt / (1 - pt)
n_d <- floor(N * prev); n_nd <- N - n_d
se_fixed <- sqrt((n_d/N)^2 * Se*(1-Se)/n_d + (n_nd/N)^2 * w^2 * Sp*(1-Sp)/n_nd)
set.seed(SEED)
nd_r <- rbinom(2e5, N, prev)
NB_r <- (rbinom(2e5, nd_r, Se) - w * rbinom(2e5, N - nd_r, 1 - Sp)) / N
sd_cohort <- sd(NB_r)
cat(sprintf("SE used under fixed margins   = %.4f   (article: 0.0077)   %s\n",
            se_fixed, ok(near(se_fixed, 0.0077, 5e-4))))
cat(sprintf("True cohort sampling SD       = %.4f   (article: 0.0174)   %s\n",
            sd_cohort, ok(near(sd_cohort, 0.0174, 1e-3))))
cat(sprintf("Ratio (times too small)       = %.2f    (article: 2.25)     %s\n",
            sd_cohort / se_fixed, ok(near(sd_cohort/se_fixed, 2.25, 0.05))))

## ---- 4. Cohort variance agrees with Monte Carlo (Table 2) ----
sep("4. Closed-form cohort variance vs Monte Carlo")
P1 <- prev * Se; P2 <- (1 - prev) * (1 - Sp)
analytic <- sqrt((P1*(1-P1) + w^2*P2*(1-P2) + 2*w*P1*P2) / N)
cat(sprintf("analytic = %.5f  MC = %.5f  ratio = %.4f   (article: ~1.002, ~0.2%%)   %s\n",
            analytic, sd_cohort, analytic / sd_cohort,
            ok(near(analytic / sd_cohort, 1, 0.01))))

## ---- 5. Apparent sensitivity against an imperfect reference (Results) ----
sep("5. Imperfect reference: apparent sensitivity and its bias")
ir <- ss_imperfect_ref(B = 6000, seed = SEED, sensitivity_table = FALSE)
app <- ir$mc_validation$se_apparent[ir$mc_validation$scenario == "adjusted"]
cat(sprintf("apparent Se = %.3f   bias = %.3f   (article: 0.764, -0.086)   %s\n",
            app, app - 0.85, ok(near(app, 0.764, 0.01))))

## ---- 6. Required N as a nested sequence of uncertainty sources (Figure 4) ----
sep("6. Figure 4: five bars, each adding one source of uncertainty")
# Shared design: Se ~ Beta(17,3), Sp ~ Beta(18,2), prevalence ~ Beta(4,16)
# (mean 0.20), delta_se = 0.07, delta_sp = 0.05, target assurance 0.80,
# loss_rate = 0, B = 20000. Step 1 is the deterministic Buderer formula;
# steps 2-5 all go through ss_unified(), each turning on exactly one more
# constraint than the step before it. See make_manuscript_assets.R for the
# full derivation, including why each N_range below is wide enough that the
# search converges well inside it rather than exhausting the range.
PRIOR_SE <- c(17, 3); PRIOR_SP <- c(18, 2); PRIOR_PREV <- c(4, 16)
DELTA_SE <- 0.07; DELTA_SP <- 0.05; B_MAIN <- 20000

N_f4_1 <- ceiling(max(buderer_n(0.85, 0.07) / 0.20,
                       buderer_n(0.90, 0.05) / 0.80))
N_f4_2 <- suppressWarnings(ss_unified(
  prior_se = PRIOR_SE, prior_sp = PRIOR_SP, prior_prev = PRIOR_PREV,
  Se_ref = 1, Sp_ref = 1, loss_rate = 0, delta_se = DELTA_SE, delta_sp = DELTA_SP,
  delta_auc = 0, check_nb = FALSE, N_range = seq(750, 1000, by = 10),
  B = B_MAIN, seed = SEED))$n_total
N_f4_3 <- suppressWarnings(ss_unified(
  prior_se = PRIOR_SE, prior_sp = PRIOR_SP, prior_prev = PRIOR_PREV,
  Se_ref = 1, Sp_ref = 1, loss_rate = 0, delta_se = DELTA_SE, delta_sp = DELTA_SP,
  delta_auc = 0.06, check_nb = FALSE, N_range = seq(750, 1000, by = 10),
  B = B_MAIN, seed = SEED))$n_total
N_f4_4 <- suppressWarnings(ss_unified(
  prior_se = PRIOR_SE, prior_sp = PRIOR_SP, prior_prev = PRIOR_PREV,
  Se_ref = 0.90, Sp_ref = 0.95, loss_rate = 0, delta_se = DELTA_SE, delta_sp = DELTA_SP,
  delta_auc = 0.06, check_nb = FALSE, N_range = seq(1050, 1400, by = 10),
  B = B_MAIN, seed = SEED))$n_total
N_f4_5_res <- suppressWarnings(ss_unified(
  prior_se = PRIOR_SE, prior_sp = PRIOR_SP, prior_prev = PRIOR_PREV,
  Se_ref = 0.90, Sp_ref = 0.95, loss_rate = 0, delta_se = DELTA_SE, delta_sp = DELTA_SP,
  delta_auc = 0.06, check_nb = TRUE, N_range = seq(1950, 2600, by = 10),
  B = B_MAIN, seed = SEED))
N_f4_5 <- N_f4_5_res$n_total

got_f4 <- c(N_f4_1, N_f4_2, N_f4_3, N_f4_4, N_f4_5)
exp_f4 <- c(500, 844, 844, 1155, 2351)
lab_f4 <- c("1. Buderer (deterministic)", "2. + joint Se/Sp precision",
            "3. + AUC gate", "4. + imperfect reference", "5. + net benefit")
for (i in seq_along(got_f4))
  cat(sprintf("%-28s = %5d   (article: %4d)   %s\n",
      lab_f4[i], got_f4[i], exp_f4[i], ok(got_f4[i] == exp_f4[i])))
cat("Note: step 2 -> step 3 is expected to be FLAT (+0): the AUC gate does\n",
    "not exclude any additional replicate at this operating point.\n", sep = "")

# Step 5 also reports the achievable ceiling of the check_nb criterion: the
# largest joint assurance the CI-based net-benefit criterion could ever
# reach as N -> Inf under this step's priors, Se_ref, Sp_ref and pt_range
# (default c(0.15, 0.40), not overridden above). This is a new diagnostic
# in 0.6.0 (see NEWS.md), not a search result, so it does not move with N --
# it depends only on the priors, Se_ref, Sp_ref, pt_range and the seed, not
# on N_range or on the search's own B (see ?ss_unified). That is also why
# the same ceiling can be read off the full Step-5 search above for the
# article's own pt_range, and off the cheapest possible call -- a
# single-value N_range and a small B -- for the other five threshold
# ranges reported in Table 5, checked below.
#
# Tolerance: nb_assurance_ceiling()'s default B_ceiling = 30000000 draws
# (raised from 200000 after the ceiling values were found to be
# systematically off by ~0.001 and correlated across pt_range -- see
# NEWS.md) is documented (?ss_unified, internal nb_assurance_ceiling) to
# keep the Monte Carlo error of the ceiling below about 1e-4 for ceilings
# in the typical 0.7-0.95 range; since both the seed and B_ceiling are
# fixed, the value is in fact exactly reproducible run to run, but the
# declared tolerance below (2e-3) also absorbs the +/-0.0005 from rounding
# the article's published figures to three decimals.
NB_CEILING_TOL <- 2e-3
nb_ceiling_for <- function(pt_range) {
  suppressWarnings(ss_unified(
    prior_se = PRIOR_SE, prior_sp = PRIOR_SP, prior_prev = PRIOR_PREV,
    Se_ref = 0.90, Sp_ref = 0.95, loss_rate = 0,
    delta_se = DELTA_SE, delta_sp = DELTA_SP, delta_auc = 0.06,
    check_nb = TRUE, pt_range = pt_range, target_assurance = 0.80,
    N_range = 500, B = 100, seed = SEED
  ))$nb_ceiling
}

table5_pt_ranges <- list(c(0.15, 0.40), c(0.20, 0.30), c(0.15, 0.25),
                          c(0.10, 0.40), c(0.10, 0.30), c(0.05, 0.50))
table5_labels <- c("0.15-0.40 (used in this article)", "0.20-0.30",
                    "0.15-0.25", "0.10-0.40", "0.10-0.30", "0.05-0.50")
table5_expected <- c(0.880, 0.969, 0.956, 0.695, 0.753, 0.000)
# The first range reuses the ceiling already returned by the full Step-5
# search above instead of recomputing it, since the two are identical by
# construction (see note above).
table5_got <- c(N_f4_5_res$nb_ceiling,
                 vapply(table5_pt_ranges[-1], nb_ceiling_for, numeric(1)))
for (i in seq_along(table5_got))
  cat(sprintf("Table 5 ceiling, pt_range = %-33s = %.3f   (article: %.3f)   %s\n",
      table5_labels[i], table5_got[i], table5_expected[i],
      ok(near(table5_got[i], table5_expected[i], NB_CEILING_TOL))))

## ---- 7. Harmonized comparison of five planning approaches (new table) ----
sep("7. Harmonized comparison table: five methods, common assumptions")
# Unlike Figure 4 (a nested sequence, all but one bar computed by
# ss_unified()), this table calls each method's own native function once,
# under the same five common parameters where each method has the
# corresponding argument.
E_SE <- PRIOR_SE[1] / sum(PRIOR_SE); E_SP <- PRIOR_SP[1] / sum(PRIOR_SP)
E_PREV <- PRIOR_PREV[1] / sum(PRIOR_PREV)

N_h_bud <- N_f4_1  # identical configuration to Figure 4 step 1
N_h_bam <- suppressWarnings(bam_sample_size(
  prior_se = PRIOR_SE, prior_sp = PRIOR_SP, prior_prev = PRIOR_PREV,
  delta_se = 2 * DELTA_SE, delta_sp = 2 * DELTA_SP, target_assurance = 0.80,
  method = "exact", B = B_MAIN, seed = SEED))$N_total
N_h_joi <- suppressWarnings(joint_sample_size(
  Se = E_SE, Sp = E_SP, delta_se = DELTA_SE, delta_sp = DELTA_SP, prev = E_PREV,
  design = "cohort", target_prob = 0.80, B = B_MAIN, seed = SEED))$n_total
N_h_imp_res <- ss_imperfect_ref(
  Se = E_SE, Sp = E_SP, d_se = DELTA_SE, d_sp = DELTA_SP, prev = E_PREV,
  Se_ref = 0.90, Sp_ref = 0.95, loss_rate = 0, B = 0, seed = SEED)
N_h_imp <- N_h_imp_res$n_total  # estimand = "apparent" (default): the
                                 # naive-analysis N, not the corrected one
N_h_uni <- N_f4_5  # identical configuration to Figure 4 step 5

got_h <- c(N_h_bud, N_h_bam, N_h_joi, N_h_imp, N_h_uni)
exp_h <- c(500, 672, 580, 732, 2351)
lab_h <- c("Buderer (classical)", "BAM (exact mode)", "Joint Se/Sp + AUC",
           "Imperfect reference (apparent)", "Unified (all sources active)")
for (i in seq_along(got_h))
  cat(sprintf("%-28s = %5d   (article: %4d)   %s\n",
      lab_h[i], got_h[i], exp_h[i], ok(got_h[i] == exp_h[i])))

## ---- 8. Imperfect reference: the corrected estimand (new in 0.6.0) ----
sep("8. Imperfect reference: apparent vs corrected estimand")
# ss_imperfect_ref() always computes and returns BOTH estimands regardless
# of which one `estimand` selects for the generic n_total slot checked in
# section 7 above. This checks the OTHER one -- the misclassification-
# corrected Se/Sp -- from the same call already made above (N_h_imp_res),
# under the identical common configuration used throughout this script.
cat(sprintf("N_corrected (imperfect ref) = %5d   (article: %4d)   %s\n",
            N_h_imp_res$N_corrected, 1235,
            ok(N_h_imp_res$N_corrected == 1235)))

## ---- Simulated data (optional export) ----
# The article analyses no empirical data. Each quantity above is computed from
# samples drawn inside the package from the stated parameters under seed 2026,
# so the parameters plus the seed fully determine the data. If a tangible file
# is required, the block below writes the underlying samples to CSV; rerunning
# this script reproduces them byte for byte.
if (isTRUE(DUMP_DATA)) {
  sep("Exporting the simulated samples to CSV")
  dir.create("validation/simulated_data", showWarnings = FALSE, recursive = TRUE)

  # (a) Buderer validation: B binomial samples of successes among 100 diseased
  set.seed(SEED)
  buderer_draws <- rbinom(4000, 100, 0.85)
  write.csv(data.frame(replicate = seq_along(buderer_draws),
                       n_diseased = 100, true_Se = 0.85,
                       observed_positives = buderer_draws),
            "validation/simulated_data/buderer_validation_samples.csv",
            row.names = FALSE)

  # (b) Cohort net benefit: the 2x2 cell counts per replicate at N = 500
  set.seed(SEED)
  nd <- rbinom(2e5, 500, 0.20)
  TP <- rbinom(2e5, nd, 0.85); FP <- rbinom(2e5, 500 - nd, 0.10)
  write.csv(data.frame(replicate = seq_len(2e5), N = 500,
                       n_diseased = nd, n_nondiseased = 500 - nd,
                       true_positives = TP, false_positives = FP,
                       false_negatives = nd - TP,
                       true_negatives = (500 - nd) - FP),
            "validation/simulated_data/cohort_netbenefit_samples.csv",
            row.names = FALSE)

  cat("Written to validation/simulated_data/:\n",
      " - buderer_validation_samples.csv  (4,000 replicates)\n",
      " - cohort_netbenefit_samples.csv   (200,000 replicates)\n", sep = "")
}

cat("\n=== reproduce_manuscript done ===\n")
