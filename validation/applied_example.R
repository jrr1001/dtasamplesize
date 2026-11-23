#!/usr/bin/env Rscript
# =============================================================================
# applied_example.R
#
# Worked example with published data: planning a validation study of a rapid
# antigen test (Panbio COVID-19 Ag, primary care, Spain) from a published
# pilot, with the exact Bayesian assurance method of dtasamplesize 0.7.0.
#
# Purpose
#   Full planning flow: initial prior and its compatibility with
#   the pilot; design prior from the pilot (power prior, fixed weight a0);
#   analysis prior (use_in); exact sample size (one-sided and two-sided,
#   joint and marginal); sensitivity (a0, width targets, assurance target,
#   prevalence prior, local Hellinger contour); frequentist complement with
#   fixed values; post-study prior-data conflict check with the validation
#   data (computed here because the vignette shows it; misspecified_priors.R studies it
#   systematically).
#
# Usage
#   Rscript applied_example.R <output_dir>
#   Writes ONLY to <output_dir> (created if missing): log, CSV tables, PNG.
#   Libraries: if DTA_PKG_LIB and/or DTA_EXTRA_LIB are set they are prepended
#   to .libPaths() (in that order).
#   Single core, no parallel code. Everything is exact (no Monte Carlo).
#
# Data (read from the full texts of the two articles in PubMed Central)
#   Pilot: Albert E, et al. Clin Microbiol Infect 2021;27(3):472.e7-472.e10,
#     doi:10.1016/j.cmi.2020.11.004, PMC7662075. Abstract and Results, first
#     paragraph: 412 patients, 43 positive by RT-PCR and antigen test, 358
#     negative by both, 11 RT-PCR+/antigen-. FP = 0 is derived
#     (43 + 358 + 11 = 412; specificity 100% stated).
#     TP = 43, FN = 11, FP = 0, TN = 358. Single-gate cohort; prevalence 54/412.
#   Validation: Bulilete O, et al. J Infect 2021;82(3):391-398,
#     doi:10.1016/j.jinf.2021.02.014, PMC7881288. Results, "Enrollment and
#     characteristics" (1369 enrolled, 140 RT-PCR positive, 7 RT-PCR
#     excluded) and "False-negative and false-positive results" (40 of 140
#     false negatives; 2 of 1222 false positives); Table 2 of Bulilete et al., global row
#     (N = 1362, Se 71.4%, Sp 99.8%). TP = 100, FN = 40, FP = 2, TN = 1220.
#     Published frequentist sample size: 927 (Methods, Statistical analysis:
#     prevalence 15%, error 5%, Se 90%); enrolled 1369; analysed 1362.
#     Manufacturer claims quoted in Background: Se 93.3%, Sp 99.4%; product
#     profile thresholds quoted there (Se >= 80%; a 97% threshold that the
#     text also calls "sensitivity", an evident erratum for specificity).
#   Targets of the one-sided example: semi-widths 0.10 (Se) and 0.05 (Sp),
#     95% one-sided, assurance 0.80, joint (the MHRA-based values used by
#     Williamson et al. 2023, Diagn Progn Res 7:17, in their illustrative
#     example).
#
# Conventions (named explicitly; each is printed in the log)
#   * Width criterion "lower_half_width": posterior median minus the
#     one-sided lower bound at the (1 - alpha_ci) quantile, alpha_ci = 0.95.
#     delta_se = 0.10, delta_sp = 0.05.
#   * Width criterion "full_width" (two-sided): full width of the equal-tailed
#     95% credible interval. "Equivalent" two-sided target declared here:
#     full width = 2 x the one-sided semi-width (0.20 for Se, 0.10 for Sp).
#     It is a convention, not an equivalence of the two intervals.
#   * Objective: joint (Se and Sp widths both met in the same study) unless
#     "marginal" (Se only, assurance_for = "se").
#   * Success rule: width <= delta, compared exactly by the package; the
#     script compares with TOL = 1e-12 (width/target and N/A comparisons).
#   * Sample size N is the TOTAL N (diseased + non-diseased), smallest N
#     whose exact assurance A(N) >= 0.80; A(N - 1) is reported.
#   * Design prior: pilot-updated power prior, base Beta(1,1) for Se, Sp and
#     prevalence, cohort pilot (prevalence updated). Base a0 = 0.5.
#   * Analysis prior: principal = use_in "design" (analysis prior = base
#     Beta(1,1); the new study is analysed without borrowing the pilot);
#     alternative = use_in "both" (the updated prior is also the analysis
#     prior). See the interpretation section of the log for the reason.
#   * Initial-prior check with the manufacturer-centred prior: ILLUSTRATIVE
#     effective size NEFF = 10 (Se Beta(9.33, 0.67), Sp Beta(9.94, 0.06));
#     the author must confirm it. 20 is shown as a variant.
#   * Conflict rule: tail probabilities of the exact Beta-Binomial prior
#     predictive, flag when p <= 0.05 (one tail) or the two-sided value
#     <= 0.05 (bam_prior_data_conflict defaults).
#   * Hellinger: H^2 = 1 - Bhattacharyya coefficient, epsilon = 0.01 (declared
#     here; Wilson et al. used 0.00354 around informative priors), 24
#     equally spaced angles in [-pi, pi) (n_phi = 24, the package default is
#     72; 24 keeps the run time acceptable), n_ref = N of the principal
#     scenario, one initial prior perturbed at a time.
#   * Frequentist complement: Buderer, n = ceil(z^2 p (1-p) / d^2) per arm,
#     N = max(n_diseased / prev, n_nondiseased / (1 - prev)), rounded up.
#     One-sided: z = qnorm(0.95). Se and prevalence from the pilot
#     (43/54, 54/412). Sp from the pilot is 1 (FP = 0) and the Wald-type
#     formula degenerates (n = 0), so a fixed plausible Sp is declared: the
#     lower limit of the two-sided 95% Wilson interval of the pilot (358/358);
#     variants 0.97 and 0.994.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: Rscript applied_example.R <output_dir>")
outdir <- args[1L]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

libs <- c(Sys.getenv("DTA_PKG_LIB"), Sys.getenv("DTA_EXTRA_LIB"))
libs <- libs[nzchar(libs)]
if (length(libs)) .libPaths(c(libs, .libPaths()))
suppressPackageStartupMessages(library(dtasamplesize))

TOL <- 1e-12
t_start <- Sys.time()
logcon <- file(file.path(outdir, "applied_example.log"), open = "wt", encoding = "UTF-8")
sink(logcon, split = TRUE)
on.exit({ sink(); close(logcon) }, add = TRUE)

hdr <- function(txt) cat("\n", strrep("=", 78), "\n", txt, "\n", strrep("=", 78), "\n", sep = "")
say <- function(...) cat(sprintf(...), "\n", sep = "")
# run an expression, log its warnings instead of letting them go to stderr
logw <- function(expr, tag = "") {
  withCallingHandlers(expr, warning = function(w) {
    cat("WARNING", if (nzchar(tag)) paste0(" [", tag, "]") else "", ": ",
        conditionMessage(w), "\n", sep = "")
    invokeRestart("muffleWarning")
  })
}
wcsv <- function(x, name) write.csv(x, file.path(outdir, name), row.names = FALSE)

hdr("applied_example.R: planning a Panbio validation study from a published pilot")
say("date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("%s", R.version.string)
say("packageVersion(dtasamplesize): %s", as.character(packageVersion("dtasamplesize")))
say("find.package(dtasamplesize): %s", find.package("dtasamplesize"))
say("TOL = %g", TOL)

# ---------------------------------------------------------------- data ------
hdr("0. Data and published figures")
pilot <- list(tp = 43, fn = 11, tn = 358, fp = 0)
valid <- list(tp = 100, fn = 40, tn = 1220, fp = 2)
stopifnot(sum(unlist(pilot)) == 412, sum(unlist(valid)) == 1362)
say("Pilot (Albert 2021): TP=%d FN=%d FP=%d TN=%d n=%d; Se=%.4f Sp=%.4f prev=%.4f",
    pilot$tp, pilot$fn, pilot$fp, pilot$tn, sum(unlist(pilot)),
    pilot$tp / (pilot$tp + pilot$fn), pilot$tn / (pilot$tn + pilot$fp),
    (pilot$tp + pilot$fn) / sum(unlist(pilot)))
say("Validation (Bulilete 2021): TP=%d FN=%d FP=%d TN=%d n=%d; Se=%.4f Sp=%.4f prev=%.4f",
    valid$tp, valid$fn, valid$fp, valid$tn, sum(unlist(valid)),
    valid$tp / (valid$tp + valid$fn), valid$tn / (valid$tn + valid$fp),
    (valid$tp + valid$fn) / sum(unlist(valid)))
N_PUBLISHED <- 927L; N_ENROLLED <- 1369L; N_ANALYSED <- 1362L
say("Published: frequentist N = %d (Bulilete Methods); enrolled %d; analysed %d",
    N_PUBLISHED, N_ENROLLED, N_ANALYSED)

flat <- c(1, 1)
pri <- function(a0, use_in, base_se = flat, base_sp = flat, base_prev = flat)
  do.call(bam_prior_from_pilot, c(pilot, list(
    base_prior_se = base_se, base_prior_sp = base_sp,
    base_prior_prev = base_prev, a0 = a0, use_in = use_in)))

# --------------------------------------------- 1. initial prior and pilot ----
hdr("1. Initial prior and its compatibility with the pilot")
cat("Base case: flat Beta(1,1) for Se, Sp and prevalence. The compatibility\n",
    "check is OMITTED for a flat prior: with Beta(1,1) the predictive is uniform\n",
    "and the smallest attainable tail is 1/(n+1), so it cannot flag anything\n",
    "informative. The function is called once below only to show\n",
    "its warning.\n", sep = "")
invisible(logw(bam_prior_data_conflict(tp = pilot$tp, fn = pilot$fn, tn = pilot$tn,
                                       fp = pilot$fp, prior_se = flat, prior_sp = flat),
               "flat initial prior"))

cat("\nSecondary illustration: initial prior centred on the manufacturer claim\n",
    "(Se 0.933, Sp 0.994; Bulilete Background), ILLUSTRATIVE effective size\n",
    "NEFF declared as an assumption (the author must confirm it).\n", sep = "")
manuf <- function(neff) list(se = c(0.933, 1 - 0.933) * neff, sp = c(0.994, 1 - 0.994) * neff)
init_chk <- do.call(rbind, lapply(c(10, 20), function(neff) {
  m <- manuf(neff)
  r <- logw(bam_prior_data_conflict(tp = pilot$tp, fn = pilot$fn, tn = pilot$tn,
                                    fp = pilot$fp, prior_se = m$se, prior_sp = m$sp),
            paste0("manufacturer prior, NEFF=", neff))
  r <- as.data.frame(r)
  data.frame(neff = neff, r[, intersect(c("quantity", "reference", "x_obs", "n", "a", "b",
                                           "p_lower", "p_upper", "p_two_sided",
                                           "flag_lower", "flag_upper", "flag_two_sided"),
                                         names(r))])
}))
print(init_chk, digits = 4, row.names = FALSE)
wcsv(init_chk, "applied_example_initial_prior_check.csv")
cat("Reading: the check reviews the INITIAL prior against the pilot; it says\n",
    "nothing about a0 (heterogeneity between pilot and study is not observed\n",
    "when planning, and a conflict does not say which source fails).\n", sep = "")

# ------------------------------------------------ 2. design prior (pilot) ----
hdr("2. Design prior from the pilot (power prior, fixed weight a0)")
cat("a0 base = 0.5, grid 0.25, 0.5, 1. Rationale (similar populations): the\n",
    "pilot enrolled only symptomatic patients, children included (85 of 412);\n",
    "the validation enrolled adults only, with 54.8% close contacts and 36.7%\n",
    "referred for symptoms (Bulilete 2021), where Se was lower in asymptomatic\n",
    "people. Same test, sample type, operators, reference test and country,\n",
    "contiguous period. Partial similarity: a0 < 1 (0.5 as the base case).\n", sep = "")
P_design <- pri(0.5, "design")
P_both   <- pri(0.5, "both")
print(P_design)
cat("\n")
print(P_both)

cat("\nAnalysis prior: principal = use_in 'design' (analysis prior = base\n",
    "Beta(1,1), design prior = pilot-updated). Alternative = 'both'.\n", sep = "")

# ------------------------------------------------ 3. exact sample size -------
hdr("3. Exact sample size (principal table)")
COMMON <- list(delta_se = 0.10, delta_sp = 0.05, width_type = "lower_half_width",
               alpha_ci = 0.95, target_assurance = 0.80)
main_scen <- list(
  design_onesided_joint = list(priors = P_design),
  both_onesided_joint   = list(priors = P_both),
  design_twosided_joint = list(priors = P_design, width_type = "full_width",
                               delta_se = 0.20, delta_sp = 0.10),
  both_twosided_joint   = list(priors = P_both, width_type = "full_width",
                               delta_se = 0.20, delta_sp = 0.10),
  design_onesided_se    = list(priors = P_design, assurance_for = "se"),
  design_twosided_se    = list(priors = P_design, assurance_for = "se",
                               width_type = "full_width", delta_se = 0.20))
main_desc <- c(
  design_onesided_joint = "PRINCIPAL: one-sided 95%, semi-widths 0.10/0.05, joint, analysis prior = base (use_in design)",
  both_onesided_joint   = "as principal, analysis prior = updated prior (use_in both)",
  design_twosided_joint = "two-sided 95%, full widths 0.20/0.10 (= 2 x semi-width), joint, use_in design",
  both_twosided_joint   = "two-sided 95%, full widths 0.20/0.10, joint, use_in both",
  design_onesided_se    = "one-sided 95%, Se only (marginal), semi-width 0.10, use_in design",
  design_twosided_se    = "two-sided 95%, Se only (marginal), full width 0.20, use_in design")
t0 <- Sys.time()
main_tab <- do.call(bam_sensitivity, c(list(main_scen), COMMON)) |> logw("main table")
say("main table time: %.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs")))
main_tab$description <- unname(main_desc[main_tab$scenario])
stopifnot(all(main_tab$target_reached))
stopifnot(all(main_tab$A_N >= 0.80 - TOL), all(main_tab$A_N_minus_1 < 0.80 + TOL))
print(main_tab[, c("scenario", "width_type", "assurance_for", "delta_se", "delta_sp",
                   "N_total", "A_N", "A_N_minus_1")], digits = 7, row.names = FALSE)
wcsv(main_tab, "applied_example_main_table.csv")
N_main <- main_tab$N_total[main_tab$scenario == "design_onesided_joint"]
say("\nPrincipal N = %d; published frequentist N = %d; enrolled = %d; analysed = %d",
    N_main, N_PUBLISHED, N_ENROLLED, N_ANALYSED)
pubcmp <- data.frame(
  source = c("principal exact N (this script)", "Bulilete 2021 Methods (frequentist)",
             "Bulilete 2021 enrolled", "Bulilete 2021 analysed"),
  N = c(N_main, N_PUBLISHED, N_ENROLLED, N_ANALYSED),
  note = c("one-sided 95%, 0.10/0.05, joint, a0 = 0.5, use_in design",
           "published; prevalence 15%, error 5%, Se 90% (Methods, Statistical analysis)",
           "published (Results, Enrollment)", "published (Table 2 of Bulilete et al., global row)"))
wcsv(pubcmp, "applied_example_published_comparison.csv")

# ------------------------------------------------ 4. sensitivity -------------
hdr("4. Sensitivity (bam_sensitivity)")
base_args <- COMMON
sens <- list(); blk <- character()
add <- function(block, name, sc) { sens[[name]] <<- sc; blk[name] <<- block }
for (a0 in c(0.25, 0.5, 1)) for (u in c("design", "both"))
  add("a0 x use_in", sprintf("a0_%s_%s", a0, u), list(priors = pri(a0, u)))
for (d in c(0.05, 0.075, 0.125, 0.15))
  add("Se semi-width", sprintf("dse_%s", d), list(priors = P_design, delta_se = d))
for (d in c(0.025, 0.075))
  add("Sp semi-width", sprintf("dsp_%s", d), list(priors = P_design, delta_sp = d))
for (tg in c(0.70, 0.90))
  add("assurance target", sprintf("target_%s", tg), list(priors = P_design, target_assurance = tg))
# prevalence prior: (i) alternative INITIAL prior, then updated with the pilot
for (neff in c(20, 100))
  add("prevalence prior", sprintf("prev_init15_n%d_updated", neff),
      list(priors = pri(0.5, "design", base_prev = c(0.15, 0.85) * neff)))
# (ii) design prevalence prior replaced (no pilot update of prevalence)
repl <- function(centre, neff) list(
  prior_se = P_design$prior_se, prior_sp = P_design$prior_sp,
  prior_prev = c(centre, 1 - centre) * neff,
  analysis_prior_se = P_design$analysis_prior_se,
  analysis_prior_sp = P_design$analysis_prior_sp)
add("prevalence prior", "prev_replaced_15_n100", repl(0.15, 100))
add("prevalence prior", "prev_replaced_10_n100_hindsight", repl(0.10, 100))
add("prevalence prior", "prev_replaced_05_n100", repl(0.05, 100))
t0 <- Sys.time()
sens_tab <- do.call(bam_sensitivity, c(list(sens), base_args)) |> logw("sensitivity")
say("sensitivity time: %.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs")))
sens_tab$block <- unname(blk[sens_tab$scenario])
sens_tab$N_change_vs_principal <- sens_tab$N_total - N_main
keep <- c("block", "scenario", "delta_se", "delta_sp", "target_assurance",
          "prior_prev_a", "prior_prev_b", "analysis_prior_se_a", "analysis_prior_sp_a",
          "N_total", "A_N", "A_N_minus_1", "target_reached", "N_change_vs_principal")
print(sens_tab[, keep], digits = 6, row.names = FALSE)
wcsv(sens_tab, "applied_example_sensitivity.csv")
cat("Notes: 'prev_init15_*' = flat prevalence prior replaced by a prior centred on\n",
    "15% (the value assumed by Bulilete et al.) with effective size 20 or 100,\n",
    "THEN updated with the pilot (a0 = 0.5). 'prev_replaced_*' = design prevalence\n",
    "prior set directly, no pilot update of prevalence. The 10% row is hindsight\n",
    "(the validation observed 10.3%) and is shown only as an illustration.\n", sep = "")

# ---- Hellinger local sensitivity -------------------------------------------
hdr("4b. Local sensitivity to the initial prior (Hellinger contour)")
EPS <- 0.01; NPHI <- 24L
say("H2 = 1 - Bhattacharyya; epsilon = %g; n_phi = %d equally spaced angles; n_ref = %d; a0 = 0.5; use_in = design",
    EPS, NPHI, N_main)
hell <- list(); hsum <- NULL
for (wh in c("se", "sp", "prev")) {
  t0 <- Sys.time()
  h <- logw(do.call(bam_sensitivity_local, c(list(
    base_prior_se = flat, base_prior_sp = flat, base_prior_prev = flat, which = wh,
    epsilon = EPS, n_phi = NPHI, pilot = pilot, a0 = 0.5, design = "cohort",
    use_in = "design", n_ref = N_main), COMMON)), paste0("Hellinger ", wh))
  say("which = %s: %.1f s", wh, as.numeric(difftime(Sys.time(), t0, units = "secs")))
  print(h)
  hell[[wh]] <- h
  s <- h$summary
  hsum <- rbind(hsum, data.frame(which = wh, epsilon = EPS, n_phi = NPHI, n_ref = N_main,
                                 A_base = h$base$A_nref, n_star_base = h$base$n_star, s))
}
wcsv(hsum, "applied_example_hellinger_summary.csv")

# ---------------------------------------------- 5. frequentist complement ----
hdr("5. Frequentist complement with fixed values (Buderer)")
se_p <- pilot$tp / (pilot$tp + pilot$fn)
prev_p <- (pilot$tp + pilot$fn) / sum(unlist(pilot))
sp_wl <- unname(wilson_ci(pilot$tn, pilot$tn + pilot$fp, alpha = 0.05)["lower"])
say("Pilot Se = %.4f, prevalence = %.4f; pilot Sp = 1 (FP = 0): Wald-type n = 0 (degenerate).",
    se_p, prev_p)
say("Fixed Sp used: lower limit of the two-sided 95%% Wilson interval of 358/358 = %.4f", sp_wl)
fr <- function(label, se, sp, prev, dse, dsp, alpha) {
  nd <- buderer_n(se, dse, alpha); nn <- buderer_n(sp, dsp, alpha)
  data.frame(case = label, se = se, sp = sp, prev = prev, d_se = dse, d_sp = dsp,
             alpha = alpha, n_diseased = nd, n_nondiseased = nn,
             N_total = max(ceiling(nd / prev), ceiling(nn / (1 - prev))))
}
frq <- rbind(
  fr("one-sided 95%, Sp = Wilson lower limit (principal)", se_p, sp_wl, prev_p, 0.10, 0.05, 0.10),
  fr("one-sided 95%, Sp = 0.97", se_p, 0.97, prev_p, 0.10, 0.05, 0.10),
  fr("one-sided 95%, Sp = 0.994", se_p, 0.994, prev_p, 0.10, 0.05, 0.10),
  fr("two-sided 95%, half-widths 0.10/0.05, Sp = Wilson lower limit", se_p, sp_wl, prev_p, 0.10, 0.05, 0.05),
  fr("one-sided 95%, prevalence 0.15 (Bulilete assumption)", se_p, sp_wl, 0.15, 0.10, 0.05, 0.10),
  fr("REPRODUCTION of published 927: Se 0.90, d 0.05, prev 0.15, two-sided 95%", 0.90, sp_wl, 0.15, 0.05, 0.05, 0.05))
frq$N_total[6] <- ceiling(frq$n_diseased[6] / 0.15)  # Se arm only, as in the published calculation
frq$note <- c(rep("", 5), sprintf("Se arm only; published N = %d; difference %d", N_PUBLISHED,
                                  ceiling(frq$n_diseased[6] / 0.15) - N_PUBLISHED))
print(frq, digits = 4, row.names = FALSE)
wcsv(frq, "applied_example_frequentist.csv")
cat("The 927 is reproduced from the inputs quoted in Bulilete Methods with the\n",
    "Buderer formula (the exact formula they used is [not verified]; the\n",
    "reproduction is a consistency check, not a statement about their method).\n", sep = "")

# --------------------------------------------- 6. post-study conflict --------
hdr("6. Post-study prior-data conflict check with the validation data")
cat("Data: TP=100/140 diseased, TN=1220/1222 non-diseased, 140 diseased of 1362.\n",
    "Reference for inference = ANALYSIS prior; calibration of the assurance =\n",
    "DESIGN prior. With use_in = 'design' the analysis prior is flat Beta(1,1): the\n",
    "analysis-prior row is then a weak check (flat-prior caveat, warning below).\n",
    "Prevalence has no analysis prior; its row refers to the design prior.\n", sep = "")
pc <- NULL
for (a0 in c(0.25, 0.5, 1)) for (u in c("design", "both")) {
  r <- as.data.frame(logw(bam_prior_data_conflict(
    tp = valid$tp, fn = valid$fn, tn = valid$tn, fp = valid$fp, priors = pri(a0, u)),
    sprintf("post-study a0=%s use_in=%s", a0, u)))
  pc <- rbind(pc, data.frame(a0 = a0, use_in = u, r[, intersect(
    c("quantity", "reference", "same_as_design", "x_obs", "n", "a", "b", "p_lower", "p_upper",
      "p_two_sided", "flag_lower", "flag_upper", "flag_two_sided"), names(r))]))
}
print(pc, digits = 4, row.names = FALSE)
wcsv(pc, "applied_example_postconflict.csv")

# ------------------------------------------------ 7. figure ----------------
hdr("7. Assurance curves and figure")
NR <- 2:1450
curves <- list(
  "a0 = 0.5, use_in design (principal)" = pri(0.5, "design"),
  "a0 = 0.25, use_in design" = pri(0.25, "design"),
  "a0 = 1, use_in design" = pri(1, "design"),
  "a0 = 0.5, use_in both" = pri(0.5, "both"))
cv <- lapply(curves, function(p) logw(do.call(bam_assurance_curve, c(
  list(priors = p, N_range = NR), COMMON[c("delta_se", "delta_sp", "width_type",
                                          "alpha_ci", "target_assurance")])), "curve"))
chk <- cv[[1]]$curve$assurance[cv[[1]]$curve$N == N_main]
say("check: A(N=%d) from curve = %.10f; from bam_sensitivity = %.10f; |diff| = %.2e",
    N_main, chk, main_tab$A_N[1], abs(chk - main_tab$A_N[1]))
stopifnot(abs(chk - main_tab$A_N[1]) <= TOL)
for (nm in names(cv)) say("  %-40s N_first = %s", nm, as.character(cv[[nm]]$N_first))
cols <- c("#000000", "#0072B2", "#D55E00", "#009E73"); ltys <- c(1, 2, 3, 4)
png(file.path(outdir, "applied_example_assurance_curve.png"), width = 2400, height = 1500, res = 300)
par(mar = c(4.2, 4.4, 1, 1))
plot(NA, xlim = range(NR), ylim = c(0, 1), xlab = "Total sample size N",
     ylab = "Exact assurance A(N)", las = 1)
abline(h = 0.80, col = "grey50", lty = 5)
for (i in seq_along(cv)) lines(cv[[i]]$curve$N, cv[[i]]$curve$assurance, col = cols[i], lty = ltys[i], lwd = 2)
abline(v = c(N_PUBLISHED, N_ENROLLED), col = "grey30", lty = c(3, 3))
points(N_main, main_tab$A_N[1], pch = 19)
text(N_PUBLISHED, 0.12, sprintf("published\nfrequentist N = %d", N_PUBLISHED), pos = 2, cex = 0.7)
text(N_ENROLLED, 0.12, sprintf("enrolled\n%d", N_ENROLLED), pos = 4, cex = 0.7)
text(N_main, 0.80, sprintf("N = %d", N_main), pos = 1, offset = 1.2, cex = 0.8)
legend("bottomright", legend = c(names(cv), "target 0.80"), col = c(cols, "grey50"),
       lty = c(ltys, 5), lwd = 2, cex = 0.7, bty = "n", inset = c(0.02, 0.2))
dev.off()
say("figure written: applied_example_assurance_curve.png")

# ------------------------------------------------ 8. interpretation ----------
hdr("8. Interpretation")
g <- function(tab, nm, col = "N_total") tab[[col]][tab$scenario == nm]
nb <- g(main_tab, "both_onesided_joint"); n2 <- g(main_tab, "design_twosided_joint")
nse <- g(main_tab, "design_onesided_se")
cat(sprintf("1. The principal flow recommends N = %d total participants (A = %.4f; A(%d) = %.4f).\n",
            N_main, main_tab$A_N[1], N_main - 1L, main_tab$A_N_minus_1[1]),
    sprintf("   The published frequentist plan was %d and %d were enrolled (%d analysed).\n",
            N_PUBLISHED, N_ENROLLED, N_ANALYSED),
    sprintf("   The recommended N is %.0f%% of the published plan.\n", 100 * N_main / N_PUBLISHED), sep = "")
cat(sprintf("2. Analysis prior: with use_in 'both' N = %d (Sp updated prior Beta(%.0f, %.0f) already\n",
            nb, P_both$analysis_prior_sp[1], P_both$analysis_prior_sp[2]),
    "   meets the Sp width before any data are observed, see the warning in section 3),\n",
    "   so the Sp criterion stops constraining N and success is driven by the prior.\n",
    "   'design' is therefore the principal choice: the analysis does not borrow the pilot.\n", sep = "")
cat(sprintf("3. Two-sided full widths 0.20/0.10 give N = %d; marginal Se only gives N = %d (one-sided).\n",
            n2, nse),
    "   Se drives N in every scenario: Sp is close to 1 in the pilot (FP = 0).\n", sep = "")
r_a0 <- sens_tab[sens_tab$block == "a0 x use_in" & grepl("design", sens_tab$scenario), ]
say("4. a0 (use_in design): %s",
    paste(sprintf("a0=%s -> N=%d", sub("a0_([0-9.]+)_.*", "\\1", r_a0$scenario), r_a0$N_total), collapse = "; "))
r_dse <- sens_tab[sens_tab$block == "Se semi-width", ]
say("   Se semi-width: %s", paste(sprintf("%s -> %d", r_dse$delta_se, r_dse$N_total), collapse = "; "))
r_pv <- sens_tab[sens_tab$block == "prevalence prior", ]
say("   prevalence prior: %s", paste(sprintf("%s -> %d", r_pv$scenario, r_pv$N_total), collapse = "; "))
r_tg <- sens_tab[sens_tab$block == "assurance target", ]
say("   assurance target: %s", paste(sprintf("%s -> %d", r_tg$target_assurance, r_tg$N_total), collapse = "; "))
say("5. Hellinger contour (epsilon %g), n* range per perturbed initial prior: %s", EPS,
    paste(sprintf("%s: %s-%s", hsum$which, hsum$n_star_min, hsum$n_star_max), collapse = "; "))
say("6. Frequentist (Buderer, fixed values): principal N = %d (Se 0.796, Sp %.3f, prev 0.131); it ignores\n   prior uncertainty and the probability of reaching the width, so it is not comparable to an assured N.",
    frq$N_total[1], sp_wl)
say("   Reproduction of the published 927 from the Methods inputs: %d.", ceiling(frq$n_diseased[6] / 0.15))
cat("7. Post-study check: see section 6 and the CSV; interpret against the ANALYSIS prior (inferential\n",
    "   reference) and the DESIGN prior (calibration). It does not choose a0.\n", sep = "")

hdr("Done")
say("elapsed: %.1f s", as.numeric(difftime(Sys.time(), t_start, units = "secs")))
