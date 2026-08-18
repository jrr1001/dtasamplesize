# Reproduce every headline number reported in the dtasamplesize manuscript
# ------------------------------------------------------------------------
# This script regenerates, from scratch, each quantity quoted in the article,
# using the same parameters and the same random seed. Every number the article
# reports is produced here; nothing is read from a stored dataset, because the
# study analyses no empirical data: all data are simulated inside the package's
# functions from the stated parameters and are fully determined by the seed.
#
# Run after installing the package:  R -f validation/reproduce_manuscript.R
#
# Optional: set DUMP_DATA <- TRUE below to also write the underlying simulated
# samples to CSV files in validation/simulated_data/ (see "Simulated data" at
# the end). The CSVs are a convenience: they are regenerated exactly by this
# script, so the script itself is the authoritative record of the data.

library(dtasamplesize)
options(width = 110, digits = 6)

# The article reports the numbers produced by version 0.4.0. Earlier releases
# contain the validity errors that 0.3.0 corrects (see NEWS.md), so they return
# different values -- notably a net-benefit sample size of 110 instead of 240,
# and a joint sample size computed from a different AUC default. Stop early
# rather than let a stale installation look like a discrepancy in the article.
if (utils::packageVersion("dtasamplesize") < "0.3.0") {
  stop("dtasamplesize ", utils::packageVersion("dtasamplesize"), " is installed, ",
       "but this script reproduces the numbers of version 0.3.0 or later.\n",
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

## ---- 6. Required N by planning approach (Figure 4) ----
sep("6. Five planning approaches at a common prevalence of 0.20")
PREV <- 0.20; PRIOR_PREV <- c(4, 16)          # Beta(4,16): mean 0.20
N_bud <- ceiling(max(buderer_n(0.85, 0.07)/PREV, buderer_n(0.90, 0.05)/(1-PREV)))
N_bam <- suppressWarnings(bam_sample_size(B = 2000, seed = SEED,
             prior_prev = PRIOR_PREV, n_range = 20:800))$N_total_median
N_joi <- suppressWarnings(joint_sample_size(B = 2000, seed = SEED,
             prev = PREV, N_range = seq(100, 900, by = 20)))$n_total
N_imp <- ss_imperfect_ref(B = 0, prev = PREV)$N_adjusted_loss
N_uni <- suppressWarnings(ss_unified(B = 1200, seed = SEED,
             prior_prev = PRIOR_PREV, N_range = seq(200, 1600, by = 20),
             delta_auc = 0, check_nb = FALSE))$n_total
got <- c(N_bud, N_bam, N_joi, N_imp, N_uni)
exp <- c(500, 591, 580, 773, 1000)
for (i in seq_along(got))
  cat(sprintf("%-22s = %5d   (article: %4d)   %s\n",
      c("Buderer (classical)","Bayesian assurance","Joint Se/Sp + AUC gate",
        "Imperfect reference","Unified framework")[i],
      got[i], exp[i], ok(got[i] == exp[i])))

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
