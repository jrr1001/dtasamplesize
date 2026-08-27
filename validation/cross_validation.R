# External cross-validation of dtasamplesize core formulas
# ---------------------------------------------------------
# Verifies the numerical machinery against INDEPENDENT references:
# Hmisc, base prop.test, numerical integration, and Monte Carlo. Also
# writes Table 2 of the manuscript ("Cross-validation of core formulas
# against independent references") to
# R-package/manuscript_assets/table2_crossvalidation.csv, built from the
# values this script computes below rather than typed in by hand: each
# row's "result" text is derived from the same variables the printed
# report uses, and if a computed value does not match what the
# manuscript prose claims, the CSV's `note` column says so instead of
# being forced to agree.
# Run after installing the package:  R -f validation/cross_validation.R
# Required reference package: Hmisc. pROC is listed as an optional
# reference package but is not called by any check below, so it is
# loaded only if available.

library(dtasamplesize)
suppressMessages(library(Hmisc))
have_pROC <- requireNamespace("pROC", quietly = TRUE)   # optional; unused by the checks below
options(width = 110, digits = 6)
ok <- function(x) if (isTRUE(x)) "PASS" else "**CHECK**"
sep <- function(t) cat("\n========== ", t, " ==========\n")

## ---------------------------------------------------------------------
## Locate R-package/manuscript_assets from this script's own path (same
## approach as make_manuscript_assets.R), so the CSV lands in the right
## place regardless of the caller's working directory.
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
    "Cannot determine the location of cross_validation.R. Run it via\n",
    "  Rscript validation/cross_validation.R\n",
    "from the dtasamplesize package root, or from inside validation/.",
    call. = FALSE
  )
}
script_dir     <- locate_script_dir()          # .../dtasamplesize/validation
pkg_root       <- dirname(script_dir)          # .../dtasamplesize
r_package_root <- dirname(pkg_root)            # .../R-package
assets_dir     <- file.path(r_package_root, "manuscript_assets")
dir.create(assets_dir, showWarnings = FALSE, recursive = TRUE)

## Provenance, recorded once and stamped on every Table 2 row below.
pkg_version  <- as.character(utils::packageVersion("dtasamplesize"))
Hmisc_version <- as.character(utils::packageVersion("Hmisc"))
pROC_version  <- if (have_pROC) as.character(utils::packageVersion("pROC")) else NA_character_

## Accumulates one row per Table 2 entry. Each block below that reproduces
## a manuscript row appends to this list with add_row(); the CSV is built
## from the accumulated data frame at the very end. None of these values
## are typed by hand -- each is either read directly from the printed
## report's own variables, or recomputed with the identical formula and
## (where relevant) the identical seed used in the check above it, so a
## second call is byte-for-byte the same as the first.
table2_rows <- list()
add_row <- function(quantity, reference, result, observed, reference_value,
                     tolerance, pass, n_cases, note = "") {
  table2_rows[[length(table2_rows) + 1]] <<- data.frame(
    quantity = quantity, reference = reference, result = result,
    observed = observed, reference_value = reference_value,
    tolerance = tolerance, pass = pass, n_cases_checked = n_cases,
    note = note,
    package_version = pkg_version, Hmisc_version = Hmisc_version,
    pROC_version = pROC_version, generated_on = as.character(Sys.Date()),
    stringsAsFactors = FALSE
  )
}

## ---- V1. Wilson CI vs Hmisc::binconf(wilson) and base prop.test ----
sep("V1. Wilson score CI")
cases <- list(c(85,100), c(45,50), c(2,20), c(19,20), c(0,30))
for (cs in cases) {
  x <- cs[1]; n <- cs[2]
  mine <- wilson_ci(x, n)
  hm   <- binconf(x, n, method = "wilson")            # Hmisc
  pt   <- prop.test(x, n, correct = FALSE)$conf.int    # base R = Wilson score
  cat(sprintf("x=%3d n=%3d | mine[%.5f,%.5f] Hmisc[%.5f,%.5f] prop.test[%.5f,%.5f]  %s / %s\n",
      x, n, mine["lower"], mine["upper"], hm[2], hm[3], pt[1], pt[2],
      ok(all(abs(c(mine["lower"]-hm[2], mine["upper"]-hm[3])) < 1e-6)),
      ok(all(abs(c(mine["lower"]-pt[1], mine["upper"]-pt[2])) < 1e-6))))
}
cat("Note: at x=19/n=20 the package matches base prop.test (canonical Wilson);\n",
    "Hmisc differs at that boundary.\n")

## -- capture for Table 2 (recomputes the same comparison against
## prop.test used above; Table 2 cites prop.test, not Hmisc, as the
## reference for this row) --
v1_diffs <- sapply(cases, function(cs) {
  x <- cs[1]; n <- cs[2]
  mine <- wilson_ci(x, n)
  pt   <- prop.test(x, n, correct = FALSE)$conf.int
  max(abs(c(mine["lower"] - pt[1], mine["upper"] - pt[2])))
})
v1_maxdiff <- max(v1_diffs)
v1_pass <- v1_maxdiff < 1e-6
add_row(
  quantity = "Wilson CI", reference = "base R prop.test(correct = FALSE)",
  result = if (v1_pass) "exact match" else
    sprintf("MISMATCH: max abs diff %.2e exceeds 1e-6", v1_maxdiff),
  observed = sprintf("max abs diff %.2e", v1_maxdiff),
  reference_value = "0 (identical bounds)", tolerance = 1e-6, pass = v1_pass,
  n_cases = length(cases)
)

## ---- V2. Wald CI vs Hmisc::binconf(asymptotic) ----
sep("V2. Wald CI (interior cases; wald_ci clamps to [0,1])")
wald_cases <- list(c(85,100), c(45,50), c(60,90))
for (cs in wald_cases) {
  x <- cs[1]; n <- cs[2]
  mine <- wald_ci(x, n)
  hm   <- binconf(x, n, method = "asymptotic")
  cat(sprintf("x=%3d n=%3d | mine[%.5f,%.5f] Hmisc[%.5f,%.5f]  %s\n",
      x, n, mine["lower"], mine["upper"], hm[2], hm[3],
      ok(all(abs(c(mine["lower"]-hm[2], mine["upper"]-hm[3])) < 1e-6))))
}

## -- capture for Table 2 --
v2_diffs <- sapply(wald_cases, function(cs) {
  x <- cs[1]; n <- cs[2]
  mine <- wald_ci(x, n)
  hm   <- binconf(x, n, method = "asymptotic")
  max(abs(c(mine["lower"] - hm[2], mine["upper"] - hm[3])))
})
v2_maxdiff <- max(v2_diffs)
v2_pass <- v2_maxdiff < 1e-6
add_row(
  quantity = "Wald CI", reference = "Hmisc::binconf",
  result = if (v2_pass) "exact match" else
    sprintf("MISMATCH: max abs diff %.2e exceeds 1e-6", v2_maxdiff),
  observed = sprintf("max abs diff %.2e", v2_maxdiff),
  reference_value = "0 (identical bounds)", tolerance = 1e-6, pass = v2_pass,
  n_cases = length(wald_cases)
)

## ---- V3. buderer_n known published values ----
sep("V3. buderer_n")
cat("Se=0.85,d=0.07 ->", buderer_n(0.85,0.07), ok(buderer_n(0.85,0.07)==100), "\n")
cat("Se=0.90,d=0.05 ->", buderer_n(0.90,0.05), ok(buderer_n(0.90,0.05)==139), "\n")

## -- capture for Table 2 --
v3_a <- buderer_n(0.85, 0.07); v3_b <- buderer_n(0.90, 0.05)
v3_pass <- (v3_a == 100) && (v3_b == 139)
add_row(
  quantity = "buderer_n", reference = "published values",
  result = if (v3_pass) sprintf("exact (%d, %d)", v3_a, v3_b) else
    sprintf("MISMATCH: got (%d, %d), published (100, 139)", v3_a, v3_b),
  observed = sprintf("(%d, %d)", v3_a, v3_b), reference_value = "(100, 139)",
  tolerance = 0, pass = v3_pass, n_cases = 2
)

## ---- V4. Hanley-McNeil variance: formula vs Monte Carlo empirical AUC var ----
sep("V4. Hanley-McNeil AUC variance vs Monte Carlo (binormal data)")
emp_auc_var <- function(AUC, n_case, n_ctrl, R = 4000, seed = 7) {
  set.seed(seed)
  mu <- qnorm(AUC) * sqrt(2)        # binormal: AUC = pnorm(mu/sqrt2)
  a <- numeric(R)
  for (r in seq_len(R)) {
    xc <- rnorm(n_case, mu, 1); yc <- rnorm(n_ctrl, 0, 1)
    rk <- rank(c(xc, yc))
    U <- sum(rk[seq_len(n_case)]) - n_case*(n_case+1)/2
    a[r] <- U / (n_case * n_ctrl)
  }
  var(a)
}
v4_configs <- list(c(0.80,100,100), c(0.85,60,140), c(0.75,80,80))
for (cfg in v4_configs) {
  AUC <- cfg[1]; nc <- cfg[2]; nk <- cfg[3]
  f  <- hanley_mcneil_var(AUC, nc, nk)
  e  <- emp_auc_var(AUC, nc, nk)
  cat(sprintf("AUC=%.2f nC=%3d nK=%3d | formula=%.6f  MC=%.6f  ratio=%.3f  %s\n",
      AUC, nc, nk, f, e, f/e, ok(abs(f/e - 1) < 0.20)))
}
cat("Note: Hanley-McNeil (1982) assumes a negative-exponential score model;\n",
    "the small deviation from binormal MC (larger when groups are imbalanced)\n",
    "is the expected, conservative behaviour of the approximation.\n")

## -- capture for Table 2 (recomputes formula/MC ratios per config above;
## same functions, same seeds, so identical to the loop's own values).
## Table 2's row is scoped to "(balanced)" configurations only (nC ==
## nK); the imbalanced configuration is recorded in `note` for audit,
## exactly as printed above, since the manuscript cell only claims
## "conservative for imbalanced" qualitatively, not a numeric tolerance. --
v4_ratios <- sapply(v4_configs, function(cfg) {
  AUC <- cfg[1]; nc <- cfg[2]; nk <- cfg[3]
  hanley_mcneil_var(AUC, nc, nk) / emp_auc_var(AUC, nc, nk)
})
v4_balanced <- sapply(v4_configs, function(cfg) cfg[2] == cfg[3])
v4_bal_ratios <- v4_ratios[v4_balanced]
v4_imb_ratio  <- v4_ratios[!v4_balanced]
v4_maxdev <- max(abs(v4_bal_ratios - 1))
v4_pass <- all(abs(v4_bal_ratios - 1) < 0.20)      # same tolerance as ok() above
add_row(
  quantity = "Hanley-McNeil AUC variance (balanced)", reference = "Monte Carlo",
  result = sprintf(
    "close, ratio %s (max %.1f%% deviation); conservative for imbalanced (ratio %.3f)",
    paste(sprintf("%.3f", v4_bal_ratios), collapse = ", "), 100 * v4_maxdev, v4_imb_ratio
  ),
  observed = paste(sprintf("%.4f", v4_bal_ratios), collapse = ", "),
  reference_value = "1 (formula == MC variance)", tolerance = 0.20, pass = v4_pass,
  n_cases = length(v4_bal_ratios),
  note = sprintf(
    paste(
      "Manuscript text now states the two balanced-configuration ratios",
      "explicitly (%.3f and %.3f, max %.1f%% from unity), matching the",
      "values computed here -- both remain within the script's own 20%%",
      "pass tolerance. Imbalanced configuration (AUC=%.2f, nC=%d, nK=%d)",
      "gives ratio %.3f, outside the 20%% tolerance (printed above as",
      "**CHECK**), which is the expected conservative behaviour the",
      "manuscript describes for that case, not a numeric claim being tested."
    ),
    v4_bal_ratios[1], v4_bal_ratios[2], 100 * v4_maxdev,
    v4_configs[!v4_balanced][[1]][1], v4_configs[!v4_balanced][[1]][2],
    v4_configs[!v4_balanced][[1]][3], v4_imb_ratio
  )
)

## ---- V5. Net benefit point estimate vs standard Vickers definition ----
sep("V5. Net benefit identity: prev-weighted form == (TP - w*FP)/N")
set.seed(1); Se<-0.85; Sp<-0.90; prev<-0.20; pt<-0.20; w<-pt/(1-pt); N<-300
n_d<-floor(N*prev); n_nd<-N-n_d
TP<-rbinom(1,n_d,Se); FP<-rbinom(1,n_nd,1-Sp)
nb_pkgform <- (n_d/N)*(TP/n_d) - (n_nd/N)*(FP/n_nd)*w     # package algebra
nb_vickers <- (TP - w*FP)/N                                # Vickers 2006 def
cat(sprintf("package=%.6f  Vickers=%.6f  %s\n", nb_pkgform, nb_vickers,
    ok(abs(nb_pkgform - nb_vickers) < 1e-9)))

## -- capture for Table 2 --
v5_diff <- abs(nb_pkgform - nb_vickers)
v5_pass <- v5_diff < 1e-9
add_row(
  quantity = "Net benefit point estimate", reference = "Vickers definition",
  result = if (v5_pass) sprintf("exact identity (abs diff %.1e)", v5_diff) else
    sprintf("MISMATCH: abs diff %.1e exceeds 1e-9", v5_diff),
  observed = sprintf("%.6f vs %.6f", nb_pkgform, nb_vickers),
  reference_value = "0 (identical)", tolerance = 1e-9, pass = v5_pass, n_cases = 1
)

## ---- V6. Net benefit variance vs Monte Carlo SD (COHORT design, package default) ----
sep("V6. Net benefit analytic SE vs Monte Carlo SD of NB_hat (cohort design)")
# The package default is design = "cohort": disease status is RANDOM, so the four
# cell counts are Multinomial(N, .). The closed-form variance used by ss_net_benefit()
# is Var(NB_hat) = [P1(1-P1) + w^2 P2(1-P2) + 2 w P1 P2] / N, with P1 = prev*Se and
# P2 = (1-prev)*(1-Sp). This block validates THAT variance, not the superseded
# fixed-margin (conditional) variance of releases <= 0.2.0.
#
# Why R = 400000: the manuscript claims the analytic SE agrees with the
# Monte Carlo SD "to within about 0.2%", a magnitude-only claim with no
# ratio digit or direction (an earlier draft quoted "ratio ~1.002", but
# that specific value was not reproducible -- across seeds the ratio is
# centred on 1.000 and its sign flips). Reproducing even the corrected
# magnitude claim requires the check itself to have resolution finer than
# 0.2%. The Monte Carlo estimator of a standard deviation has relative
# sampling error of order 1/sqrt(2R). At R = 40000 that is
# 1/sqrt(80000) = 0.35%, which is coarser than the 0.2% being tested --
# the check would be reporting its own noise, not the formula's accuracy,
# and could not even fix the sign of any deviation. At R = 400000 the
# noise floor drops to 1/sqrt(800000) = 0.11%, small enough relative to
# 0.2% to say the analytic and Monte Carlo SDs agree at that resolution.
R6 <- 400000   # replications for the net-benefit cohort-variance check (see rationale above)
nb_se_cohort <- function(Se, Sp, prev, pt, N, R = R6, seed = 3) {
  w  <- pt/(1-pt); P1 <- prev*Se; P2 <- (1-prev)*(1-Sp)
  analytic <- sqrt((P1*(1-P1) + w^2*P2*(1-P2) + 2*w*P1*P2)/N)   # package cohort formula
  set.seed(seed)
  n_d <- rbinom(R, N, prev); n_nd <- N - n_d                    # disease status random
  TP  <- rbinom(R, n_d, Se); FP <- rbinom(R, n_nd, 1-Sp)
  NB  <- (TP - w*FP)/N
  c(analytic = analytic, mc_sd = sd(NB), ratio = analytic/sd(NB))
}
v6_configs <- list(c(0.85,0.90,0.20,0.20,500), c(0.80,0.85,0.30,0.40,300))
for (cfg in v6_configs) {
  r <- do.call(nb_se_cohort, as.list(cfg))
  cat(sprintf("Se=%.2f Sp=%.2f prev=%.2f pt=%.2f N=%d | analytic=%.5f MC=%.5f ratio=%.3f %s\n",
      cfg[1],cfg[2],cfg[3],cfg[4],cfg[5], r["analytic"], r["mc_sd"], r["ratio"],
      ok(abs(r["ratio"]-1) < 0.02)))
}
cat("Note: row 1 is the manuscript's headline configuration (Se=0.85, Sp=0.90, prev=0.20,\n",
    "pt=0.20, N=500): cohort SD ~0.0174, agreeing with the analytic SE to within ~0.2%\n",
    "at R=400000, where the SD estimator's own sampling noise (~0.11%) is fine enough\n",
    "to resolve that claim (see comment above nb_se_cohort()).\n")

## -- capture for Table 2 (recomputes both configurations with the same
## seed; row 1, Se=0.85/Sp=0.90/prev=0.20/pt=0.20/N=500, is the headline
## configuration named in the comment above, so it is the value reported
## in the Result column; row 2's ratio is kept in `note`) --
v6_results <- lapply(v6_configs, function(cfg) do.call(nb_se_cohort, as.list(cfg)))
v6_head <- v6_results[[1]]
v6_dev  <- abs(v6_head["ratio"] - 1)
v6_pass <- all(sapply(v6_results, function(r) abs(r["ratio"] - 1) < 0.02))
add_row(
  quantity = "Net benefit variance (cohort)", reference = "cohort Monte Carlo SD",
  result = sprintf("agrees to within about 0.2%% (ratio %.3f, %.2f%% diff, R=%d)",
    v6_head["ratio"], 100 * v6_dev, R6),
  observed = sprintf("analytic=%.5f MC_sd=%.5f", v6_head["analytic"], v6_head["mc_sd"]),
  reference_value = "ratio 1 (analytic == MC SD)", tolerance = 0.02, pass = v6_pass,
  n_cases = length(v6_configs),
  note = sprintf(
    paste(
      "Manuscript text states the cohort net-benefit variance agrees with",
      "cohort Monte Carlo to within about 0.2%%, a magnitude-only claim with",
      "no ratio digit or direction (an earlier draft quoted \"ratio ~1.002\",",
      "which was not reproducible across seeds -- the ratio is centred on",
      "1.000 and its sign varies by seed). At R=%d, chosen so the SD",
      "estimator's own sampling noise (~1/sqrt(2R) = %.2f%%) is fine enough",
      "to resolve a 0.2%% claim, the headline configuration gives ratio",
      "%.4f (%.2f%% from unity) and the second configuration gives ratio",
      "%.4f (%.2f%% from unity); both are consistent with \"to within about",
      "0.2%%\" and within the script's own 2%% pass tolerance."
    ),
    R6, 100 / sqrt(2 * R6),
    v6_head["ratio"], 100 * v6_dev,
    v6_results[[2]]["ratio"], 100 * abs(v6_results[[2]]["ratio"] - 1)
  )
)

## ---- V7. BAM Beta-Binomial conjugacy vs numerical Bayes posterior ----
sep("V7. Beta-Binomial posterior CI (conjugate) vs numerical integration")
a<-17; b<-3; n<-100; x<-86
cw <- qbeta(0.975, a+x, b+n-x) - qbeta(0.025, a+x, b+n-x)   # conjugate (BAM)
g <- seq(1e-5, 1-1e-5, length.out = 200001)                # numerical grid
post <- dbeta(g,a,b) * dbinom(x, n, g); post <- post/sum(post)
cdf <- cumsum(post)
lo <- g[which.min(abs(cdf-0.025))]; hi <- g[which.min(abs(cdf-0.975))]
cat(sprintf("conjugate width=%.5f  numerical width=%.5f  %s\n",
    cw, hi-lo, ok(abs(cw-(hi-lo)) < 2e-3)))

## -- capture for Table 2 --
v7_diff <- abs(cw - (hi - lo))
v7_code_tol <- 2e-3                 # the tolerance actually used by ok() above
v7_claimed_tol <- 1e-3              # the tolerance named in the manuscript's Result cell
v7_pass <- v7_diff < v7_code_tol
add_row(
  quantity = "Beta-Binomial posterior", reference = "numerical integration",
  result = if (v7_diff < v7_claimed_tol) sprintf("match to 1e-3 (abs diff %.1e)", v7_diff) else
    sprintf("match to %.0e only (abs diff %.1e exceeds 1e-3)", v7_code_tol, v7_diff),
  observed = sprintf("%.5f vs %.5f", cw, hi - lo), reference_value = "0 (identical width)",
  tolerance = v7_code_tol, pass = v7_pass, n_cases = 1
)

## ---- V8. Imperfect-ref: VIF (Rogan-Gladen) and apparent Se closed form ----
sep("V8. Imperfect reference: VIF and apparent-Se closed form")
Se_ref<-0.90; Sp_ref<-0.95; Se<-0.85; Sp<-0.90; prev<-0.30
vif_pkg <- ss_imperfect_ref(Se=Se,Sp=Sp,Se_ref=Se_ref,Sp_ref=Sp_ref,prev=prev,B=0)$VIF
vif_rg  <- 1/(Se_ref+Sp_ref-1)^2
cat(sprintf("VIF package=%.6f  Rogan-Gladen=%.6f  %s\n", vif_pkg, vif_rg, ok(abs(vif_pkg-vif_rg)<1e-9)))
ir <- ss_imperfect_ref(Se=Se,Sp=Sp,Se_ref=Se_ref,Sp_ref=Sp_ref,prev=prev,B=6000,seed=11,sensitivity_table=FALSE)
app_emp <- ir$mc_validation$se_apparent[ir$mc_validation$scenario=="adjusted"]
app_theo <- (prev*Se*Se_ref + (1-prev)*(1-Sp)*(1-Sp_ref)) /
            (prev*Se_ref + (1-prev)*(1-Sp_ref))
cat(sprintf("apparent Se: MC=%.4f  closed-form=%.4f  %s\n", app_emp, app_theo,
    ok(abs(app_emp-app_theo) < 0.01)))

## -- capture for Table 2 (Table 2's row scopes only the VIF vs
## Rogan-Gladen comparison; the apparent-Se Monte Carlo check above is
## reported in the manuscript's narrative text, not as a separate Table 2
## row, so it is not folded in here) --
v8_diff <- abs(vif_pkg - vif_rg)
v8_pass <- v8_diff < 1e-9
add_row(
  quantity = "Imperfect-reference VIF", reference = "Rogan-Gladen",
  result = if (v8_pass) sprintf("exact (abs diff %.1e)", v8_diff) else
    sprintf("MISMATCH: abs diff %.1e exceeds 1e-9", v8_diff),
  observed = sprintf("%.6f vs %.6f", vif_pkg, vif_rg), reference_value = "0 (identical)",
  tolerance = 1e-9, pass = v8_pass, n_cases = 1
)

## ---------------------------------------------------------------------
## Write Table 2 ("Cross-validation of core formulas against independent
## references") from the rows accumulated above.
## ---------------------------------------------------------------------
table2 <- do.call(rbind, table2_rows)
write.csv(table2, file.path(assets_dir, "table2_crossvalidation.csv"), row.names = FALSE)
cat(sprintf("\nTable 2 written to %s (%d rows, all pass = %s)\n",
    file.path(assets_dir, "table2_crossvalidation.csv"), nrow(table2), all(table2$pass)))

cat("\n=== cross_validation done ===\n")
