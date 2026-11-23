#!/usr/bin/env Rscript
# =============================================================================
# make_manuscript_tables_R1.R
#
# Generator of the tables of the revised article that derive from the
# validation scripts: Tables 2, 3, 5 and 6 of the main text, and the tables of Additional file 3 (S3.1 to S3.9)
# and two of Additional file 2 (S2.3, S2.4).
#
# It reads the CSV files written by the validation scripts (it does NOT redo
# their calculations, with one exception: the Monte Carlo row of Table 2,
# joint_sample_size(), is recomputed with the installed package) and checks
# every number that the text quotes against the value expected in the article
# ("anchors").
#
# Usage
#   Rscript validation/make_manuscript_tables_R1.R <out_dir> \
#     --assets=<dir> --comparators=<dir> --width=<dir> --negctl=<dir> \
#     --emptyarm=<dir> --wilson=<dir> --williamson=<dir> --bridge=<dir> \
#     --exactmc=<dir> --shortcut=<dir> --applied=<dir> --oc=<dir> --mp=<dir> \
#     --crosscheck=<dir> [--witness=<wilson2022_rule_enumeration_enum.log>]
#   <out_dir> is the first positional argument and the ONLY place where the
#   script writes (default: validation/output/ inside the package tree). The
#   input directories are named, never searched: --assets is the output of
#   make_manuscript_assets.R; --comparators of comparator_calls.R; --width of
#   width_attainment_exact.R; --negctl of negative_control_halfwidth.R;
#   --emptyarm of edge_case_empty_arm.R; --wilson of compare_wilson2022.R;
#   --williamson of compare_williamson2023.R; --bridge of
#   bridge_frequentist_bayesian.R; --exactmc of exact_vs_mc_resolution.R;
#   --shortcut of success_set_shortcut.R; --applied of applied_example.R;
#   --oc of operating_characteristics.R; --mp of misspecified_priors.R;
#   --crosscheck of oc_montecarlo_crosscheck.R.
#   --witness (optional) is the text log with the witness cases of the exact
#   enumeration against the programmed rule of Wilson et al. (S3.1).
#   Libraries: DTA_PKG_LIB and DTA_EXTRA_LIB, if set, are prepended to
#   .libPaths(). Single core. The log lists the path and SHA-256 of every
#   file read.
#
# Conventions of the numbers in the tables (each is also printed in the log)
#   * Widths: the Table 2 inputs are half-widths 0.07 (Se) and 0.05 (Sp), i.e.
#     full widths 0.14 and 0.10; two-sided intervals are equal-tailed 95%.
#   * Searches: the exact rows scan N from 2 and report the smallest N with
#     A(N) >= 0.80 (comparison with TOL = 1e-12).
#   * Selection of rows by key (column key, design/analysis, scenario), never
#     by line number.
#   * Anchors: integers must match exactly; probabilities must match the
#     ledger value rounded to its printed decimals (|difference| <= half a unit
#     of the last printed decimal); a few values are held to 1e-10 or 1e-6 as
#     stated in the line. A failed anchor is printed as **CHECK**, the script
#     finishes writing and exits with status 1.
#   * Display files (_display.csv, .md): probabilities with 4 decimals; values
#     below 0.001 in scientific notation (8.0 x 10^-6); integers as integers;
#     rows with "to" for ranges; empty cell = not applicable.
#
# Outputs (prefix R1_): R1_Table2_planning_questions, R1_Table3_reproductions,
#   R1_Table5_applied_example, R1_Table6_operating_characteristics (each as
#   .csv full precision, _display.csv and .md); R1_AF3_S3.* and R1_AF2_S2.*
#   (a .csv and a .md per table); log make_manuscript_tables_R1.log.
# =============================================================================

options(width = 140, stringsAsFactors = FALSE, dtasamplesize.warn_small_B = FALSE)

## ---- libraries, arguments, paths ------------------------------------------
libs <- c(Sys.getenv("DTA_PKG_LIB"), Sys.getenv("DTA_EXTRA_LIB")); libs <- libs[nzchar(libs)]
if (length(libs)) .libPaths(c(libs, .libPaths()))

locate_script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) == 1) return(dirname(normalizePath(sub("^--file=", "", f[1]), winslash = "/")))
  if (basename(getwd()) == "validation") return(normalizePath(getwd(), winslash = "/"))
  if (dir.exists("validation")) return(normalizePath(file.path(getwd(), "validation"), winslash = "/"))
  stop("Cannot locate the validation/ folder; run with Rscript validation/make_manuscript_tables_R1.R")
}
script_dir <- locate_script_dir()
pkg_root <- dirname(script_dir)

args_all <- commandArgs(trailingOnly = TRUE)
named <- list(); pos <- character(0)
for (a in args_all) {
  if (grepl("^--[A-Za-z0-9_]+=", a)) {
    named[[sub("^--([A-Za-z0-9_]+)=.*$", "\\1", a)]] <- sub("^--[A-Za-z0-9_]+=", "", a)
  } else pos <- c(pos, a)
}
out_dir <- if (length(pos) >= 1 && nzchar(pos[1])) pos[1] else file.path(script_dir, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_dir <- normalizePath(out_dir, winslash = "/")

REQ <- c("assets", "comparators", "width", "negctl", "emptyarm", "wilson", "williamson",
         "bridge", "exactmc", "shortcut", "applied", "oc", "mp", "crosscheck")
miss <- setdiff(REQ, names(named))
if (length(miss)) stop("Missing named arguments: ", paste0("--", miss, "=<dir>", collapse = " "))
for (k in REQ) if (!dir.exists(named[[k]])) stop("Not a directory: --", k, "=", named[[k]])

logcon <- file(file.path(out_dir, "make_manuscript_tables_R1.log"), open = "wt", encoding = "UTF-8")
sink(logcon, split = TRUE)
suppressMessages(library(dtasamplesize))

git_commit <- tryCatch({
  h <- system2("git", c("-C", shQuote(pkg_root), "rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE)
  n <- length(system2("git", c("-C", shQuote(pkg_root), "status", "--porcelain"), stdout = TRUE, stderr = FALSE))
  sprintf("%s%s", h[1], if (n > 0) sprintf(" (+%d uncommitted paths)", n) else "")
}, error = function(e) "not available", warning = function(w) "not available")
cat(sprintf("%s | %s | dtasamplesize %s at %s | tree commit %s | script make_manuscript_tables_R1.R | out %s\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), R.version.string,
            as.character(utils::packageVersion("dtasamplesize")),
            find.package("dtasamplesize"), git_commit, out_dir))
T_START <- Sys.time()

## ---- inputs (with SHA-256) ------------------------------------------------
INPUT_LOG <- data.frame(arg = character(), file = character(), path = character(), sha256 = character())
.cache <- new.env()
rd <- function(argname, file) {
  p <- normalizePath(file.path(named[[argname]], file), winslash = "/", mustWork = FALSE)
  if (!file.exists(p)) stop("Input not found: ", p)
  if (!is.null(.cache[[p]])) return(.cache[[p]])
  sha <- unname(tools::sha256sum(p))
  cat(sprintf("INPUT --%s %s sha256=%s\n", argname, p, sha))
  INPUT_LOG[nrow(INPUT_LOG) + 1L, ] <<- list(argname, file, p, sha)
  x <- utils::read.csv(p, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  .cache[[p]] <- x
  x
}

## ---- anchors --------------------------------------------------------------
ANCH <- new.env(); ANCH$n <- 0L; ANCH$pass <- 0L; ANCH$fail <- character(0)
fmtv <- function(v) paste(ifelse(is.na(v), "NA", if (is.numeric(v)) format(v, digits = 10, trim = TRUE) else as.character(v)), collapse = ", ")
# exp: expected value(s) from the ledger. dec: number of printed decimals (probabilities);
# tol: explicit absolute tolerance; neither: exact equality (integers, text).
anchor <- function(id, desc, obs, exp, dec = NA_integer_, tol = NA_real_) {
  ok <- tryCatch({
    if (length(obs) != length(exp)) FALSE
    else if (is.character(exp)) identical(as.character(obs), exp)
    else {
      o <- suppressWarnings(as.numeric(obs))
      na_ok <- all(is.na(o) == is.na(exp))
      keep <- !is.na(exp)
      d <- abs(o[keep] - exp[keep])
      lim <- if (!is.na(tol)) tol else if (!is.na(dec)) 0.5 * 10^(-dec) + 1e-12 else 0
      na_ok && all(d <= lim)
    }
  }, error = function(e) FALSE)
  ANCH$n <- ANCH$n + 1L
  if (isTRUE(ok)) ANCH$pass <- ANCH$pass + 1L else ANCH$fail <- c(ANCH$fail, id)
  cat(sprintf("ANCHOR %s %s %s (%s)\n", id, desc, if (isTRUE(ok)) "PASS" else "**CHECK**", fmtv(obs)))
  invisible(isTRUE(ok))
}

## ---- formatting and writers -----------------------------------------------
SUP <- c("0" = "⁰", "1" = "¹", "2" = "²", "3" = "³", "4" = "⁴",
         "5" = "⁵", "6" = "⁶", "7" = "⁷", "8" = "⁸", "9" = "⁹")
sci1 <- function(x) {
  e <- floor(log10(abs(x))); m <- round(x / 10^e, 1)
  if (abs(m) >= 10) { m <- m / 10; e <- e + 1 }
  paste0(formatC(m, format = "f", digits = 1), " × 10", if (e < 0) "⁻" else "",
         paste(SUP[strsplit(as.character(abs(e)), "")[[1]]], collapse = ""))
}
fp <- function(x, d = 4) vapply(x, function(v) {
  if (is.na(v)) "" else if (v != 0 && abs(v) < 1e-3) sci1(v) else formatC(v, format = "f", digits = d)
}, character(1))
fi <- function(x) vapply(x, function(v) if (is.na(v)) "" else formatC(v, format = "d", big.mark = ""), character(1))
# automatic display of a data frame column: integers as integers, other numbers as fp()
fauto <- function(x) {
  if (is.numeric(x)) {
    nn <- x[!is.na(x)]
    if (length(nn) && all(abs(nn - round(nn)) < 1e-9 & abs(nn) < 1e15)) fi(round(x)) else fp(x)
  } else if (is.logical(x)) ifelse(is.na(x), "", as.character(x)) else ifelse(is.na(x), "", as.character(x))
}
display_df <- function(df) { out <- df; out[] <- lapply(df, fauto); out }
write_csv <- function(df, name) {
  p <- file.path(out_dir, name)
  utils::write.csv(df, p, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
  cat(sprintf("WROTE %s (%d rows x %d cols)\n", name, nrow(df), ncol(df)))
}
write_md <- function(df, name, caption = NULL) {
  df[] <- lapply(df, function(x) gsub("\n", " ", gsub("|", "\\|", ifelse(is.na(as.character(x)), "", as.character(x)), fixed = TRUE)))
  hdr <- paste0("| ", paste(gsub("|", "\\|", names(df), fixed = TRUE), collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  body <- if (nrow(df)) apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")) else character(0)
  con <- file(file.path(out_dir, name), open = "wt", encoding = "UTF-8")
  writeLines(c(if (!is.null(caption)) c(caption, ""), hdr, sep, body), con); close(con)
  cat(sprintf("WROTE %s\n", name))
}
emit_main <- function(stem, full, display, caption = NULL) {
  write_csv(full, paste0(stem, ".csv"))
  write_csv(display, paste0(stem, "_display.csv"))
  write_md(display, paste0(stem, ".md"), caption)
}
emit_simple <- function(stem, df, caption = NULL) {
  write_csv(df, paste0(stem, ".csv"))
  write_md(display_df(df), paste0(stem, ".md"), caption)
}
one <- function(df, cond) { r <- df[cond, , drop = FALSE]; if (nrow(r) != 1L) stop("expected exactly one row, got ", nrow(r)); r }
num <- function(x) suppressWarnings(as.numeric(x))

# =============================================================================
# Table 2: R1_Table2_planning_questions
# =============================================================================
cat("\n== Table 2: planning questions ==\n")
t4   <- rd("assets", "table4_method_comparison.csv")
br   <- rd("bridge", "bridge_table.csv")
fig2 <- rd("assets", "figure2_exact_curve.csv")
cj   <- rd("comparators", "comparator_calls_results.csv")
TOLW <- 1e-12

bro <- function(d, a) one(br, br$design == d & br$analysis == a)

cat("Recomputing joint_sample_size() (Monte Carlo, B = 20000, seed 2026, default N_range)...\n")
t_js <- Sys.time()
js <- joint_sample_size(Se = 0.85, Sp = 0.90, delta_se = 0.07, delta_sp = 0.05, prev = 0.20,
                        design = "cohort", target_prob = 0.80, B = 20000, seed = 2026)
cat(sprintf("joint_sample_size elapsed %.1f s; fields: %s\n",
            as.numeric(difftime(Sys.time(), t_js, units = "secs")), paste(names(js), collapse = ", ")))
js_N <- js$n_total; js_p <- js$joint_prob_se_sp; js_se <- js$joint_prob_mcse
if (!is.null(js$search_type)) cat("search_type:", js$search_type, "\n")
cat(sprintf("Monte Carlo row: N = %s, joint_prob_se_sp = %.6f, MCSE = %.6f\n", js_N, js_p, js_se))

b_row <- one(t4, t4$method == "Buderer (classical)")
bam_t4 <- one(t4, t4$method == "BAM (exact mode, harmonized priors)")
joint_t4 <- one(t4, t4$method == "Joint Se/Sp (joint_sample_size)")
c4 <- one(cj, cj$id == "c_joint")
# Probability and MCSE of the integer-grid call, read from the note of the row
c4_p  <- as.numeric(sub(".*joint_prob_se_sp = ([0-9.]+).*", "\\1", c4$note))
c4_se <- as.numeric(sub(".*mcse = ([0-9.]+).*", "\\1", c4$note))
if (is.na(c4_p) || is.na(c4_se)) stop("Could not read the probability and MCSE from comparator_calls_results.csv")
f2_ok <- fig2$assurance >= 0.80 - TOLW
f2_N <- min(fig2$N[f2_ok]); f2_A <- fig2$assurance[fig2$N == f2_N]

T2 <- data.frame(
  row = 1:12,
  question = c("Expected plug-in width, each arm separately",
               rep("Probability that both widths are attained", 3),
               rep("Credible interval with a vague prior", 2),
               "Informative analysis prior",
               rep("Prior uncertainty in the design", 3),
               "Bayesian assurance, informative Sp prior",
               "Bayesian assurance, worked example with vague Sp prior"),
  design_assumption = c(rep("Fixed planning values", 7),
                        rep("Beta design priors of the example", 3),
                        "Beta design priors (Se Beta(17, 3), Sp Beta(18, 2), prevalence Beta(4, 16)); analysis = design",
                        "Beta design priors of the worked example (vague Sp prior)"),
  analysis_interval = c("Wald", "Wilson", "Wilson", "Wilson", "Jeffreys", "Uniform",
                        "Beta(17, 3), Beta(18, 2)", "Wilson", "Jeffreys", "Uniform",
                        "Beta(17, 3), Beta(18, 2)", "Equal-tailed credible, analysis = design"),
  computation = c("Closed form", "Exact enumeration", "Monte Carlo, step-10 grid (B = 20000, seed 2026)",
                  "Monte Carlo, integer grid 100:800 (B = 20000, seed 2026; validation/comparator_calls.R)",
                  "Exact enumeration", "Exact enumeration", "Exact enumeration",
                  "Exact enumeration", "Exact enumeration", "Exact enumeration",
                  "Exact enumeration", "Exact enumeration"),
  function_or_script = c("buderer_n()", "validation/bridge_frequentist_bayesian.R", "joint_sample_size()",
                         "joint_sample_size()", "validation/bridge_frequentist_bayesian.R",
                         "validation/bridge_frequentist_bayesian.R", "validation/bridge_frequentist_bayesian.R",
                         "validation/bridge_frequentist_bayesian.R", "validation/bridge_frequentist_bayesian.R",
                         "validation/bridge_frequentist_bayesian.R", "bam_sample_size(method = \"exact\")",
                         "bam_assurance_curve()"),
  N = c(b_row$N,
        bro("point", "wilson")$N_star, js_N, num(c4$value),
        bro("point", "jeffreys")$N_star, bro("point", "uniform")$N_star, bro("point", "informative")$N_star,
        bro("beta", "wilson")$N_star, bro("beta", "jeffreys")$N_star, bro("beta", "uniform")$N_star,
        bro("beta", "informative")$N_star, f2_N),
  probability_at_N = c(NA_real_,
                       bro("point", "wilson")$A_at_N_star, js_p, c4_p,
                       bro("point", "jeffreys")$A_at_N_star, bro("point", "uniform")$A_at_N_star,
                       bro("point", "informative")$A_at_N_star,
                       bro("beta", "wilson")$A_at_N_star, bro("beta", "jeffreys")$A_at_N_star,
                       bro("beta", "uniform")$A_at_N_star, bro("beta", "informative")$A_at_N_star, f2_A),
  mcse = c(NA, NA, js_se, c4_se, rep(NA, 8)),
  ledger_id = paste0("T2-", sprintf("%02d", 1:12)),
  stringsAsFactors = FALSE)
T2$N <- as.integer(T2$N)

anchor("T2-01", "Buderer N", T2$N[1], 500L)
anchor("T2-02", "point/Wilson N, A", c(T2$N[2], T2$probability_at_N[2]), c(578, 0.8012), dec = 4)
anchor("T2-02n", "point/Wilson N exact", T2$N[2], 578L)
anchor("T2-03", "joint_sample_size N", T2$N[3], 580L)
anchor("T2-03p", "joint_sample_size probability", T2$probability_at_N[3], 0.8041, dec = 4)
anchor("T2-03m", "joint_sample_size MCSE (tolerance 1e-6)", T2$mcse[3], 0.002806, tol = 1e-6)
anchor("T2-03x", "recomputed N and probability equal to table4_method_comparison.csv",
       c(T2$N[3], T2$probability_at_N[3]), c(joint_t4$N, joint_t4$assurance_achieved), tol = 1e-12)
anchor("T2-04", "integer grid 100:800 N (comparator_calls_results.csv, c_joint)", T2$N[4], 579L)
anchor("T2-04p", "integer grid 100:800 probability and MCSE (comparator_calls_results.csv, c_joint)",
       c(T2$probability_at_N[4], T2$mcse[4]), c(0.8003, 0.0028), dec = 4)
anchor("T2-05", "point/Jeffreys N, A", c(T2$N[5], T2$probability_at_N[5]), c(576, 0.8014), dec = 4)
anchor("T2-06", "point/uniform N, A", c(T2$N[6], T2$probability_at_N[6]), c(579, 0.8025), dec = 4)
anchor("T2-07", "point/informative N, A", c(T2$N[7], T2$probability_at_N[7]), c(466, 0.8002), dec = 4)
anchor("T2-08", "beta/Wilson N, A", c(T2$N[8], T2$probability_at_N[8]), c(831, 0.8004), dec = 4)
anchor("T2-09", "beta/Jeffreys N, A", c(T2$N[9], T2$probability_at_N[9]), c(824, 0.8000), dec = 4)
anchor("T2-10", "beta/uniform N, A", c(T2$N[10], T2$probability_at_N[10]), c(831, 0.8004), dec = 4)
anchor("T2-11", "beta/informative (BAM) N, A", c(T2$N[11], T2$probability_at_N[11]), c(672, 0.8003), dec = 4)
anchor("T2-11x", "BAM row equal to table4_method_comparison.csv (<= 1e-10)",
       c(T2$N[11], T2$probability_at_N[11]), c(bam_t4$N, bam_t4$assurance_achieved), tol = 1e-10)
anchor("T2-12", "worked example first N with A >= 0.80, A", c(T2$N[12], T2$probability_at_N[12]), c(678, 0.8003), dec = 4)
anchor("T2-12x", "worked example A(678) (<= 1e-10)", T2$probability_at_N[12], 0.8003489948, tol = 1e-10)
anchor("T2-13", "inputs (Se, Sp, prevalence, half-widths) in table4_method_comparison.csv",
       c(b_row$Se, b_row$Sp, b_row$prevalence, b_row$delta_se, b_row$delta_sp), c(0.85, 0.90, 0.20, 0.07, 0.05), tol = 1e-12)
anchor("T2-13b", "Monte Carlo settings B, seed", c(joint_t4$B, joint_t4$seed), c(20000, 2026))
cat("Facts for the note of Table 2 (T2-13): planning values Se 0.85, Sp 0.90, prevalence 0.20; half-widths 0.07 and 0.05 (full widths 0.14 and 0.10); target 0.80; equal-tailed 95% intervals; exact rows scan N from 2; Monte Carlo row B = 20000, seed 2026, MCSE", sprintf("%.6f", js_se), "\n")

T2d <- data.frame(
  Row = T2$row, `Planning question` = T2$question, `Design assumption` = T2$design_assumption,
  `Analysis interval` = T2$analysis_interval, Computation = T2$computation,
  `Function or script` = T2$function_or_script, N = fi(T2$N),
  `Probability at N` = ifelse(is.na(T2$probability_at_N), "", fp(T2$probability_at_N)),
  `MCSE` = ifelse(is.na(T2$mcse), "", fp(T2$mcse)), check.names = FALSE, stringsAsFactors = FALSE)
emit_main("R1_Table2_planning_questions", T2, T2d)

# =============================================================================
# Table 3: R1_Table3_reproductions
# =============================================================================
cat("\n== Table 3: reproductions ==\n")
wt <- rd("wilson", "compare_wilson2022_table.csv")
f32a <- rd("williamson", "f32a_illustrative.csv")
f32b <- rd("williamson", "f32b_figure3.csv")
f32v <- rd("williamson", "f32b_figure3_variants.csv")
f33  <- rd("williamson", "f33_package_vs_app.csv")
wk <- function(key) one(wt, wt$key == key)
wcol <- function(key, col) num(wk(key)[[col]])

L3 <- list()
add3 <- function(block, group, case, quantity, mode, inputs, published, published_source, exact,
                 rule = NA_real_, rule_s2 = NA_real_, ledger, note = "") {
  L3[[length(L3) + 1L]] <<- data.frame(block = block, group = group, case = case, quantity = quantity,
    mode_exercised = mode, inputs = inputs, published = published, published_source = published_source,
    dtasamplesize_exact = exact, programmed_rule = rule, programmed_rule_dataS2 = rule_s2,
    ledger_id = ledger, note = note, stringsAsFactors = FALSE)
}
WMODE <- "Se only (marginal); two-sided 95%; full width; exact enumeration"
addw <- function(group, case, key, quantity, ledger, note = "", mode = WMODE) {
  r <- wk(key)
  add3("Wilson et al. 2022", group, case, quantity, mode, r$inputs, num(r$published), r$published_source,
       num(r$package), num(r$emulation), num(r$dataS2_original), ledger, note)
}
# VAP
addw("VAP", "VAP (Section 5.1)", "VAP|n_star", "n*", "T3-02;T3-03;T3-04")
addw("VAP", "VAP (Section 5.1)", "VAP|A_N104", "A(104)", "T3-03")
addw("VAP", "VAP (Section 5.1)", "VAP|A_N106", "A(106)", "T3-03")
addw("VAP", "VAP (Section 5.1)", "VAP|A_N150", "A(150)", "T3-03")
# B.4
addw("B4flat", "Supplement B.4, Combination (flat analysis prior)", "B4|Combination|n_star", "n*", "T3-05")
addw("B4flat", "Supplement B.4, Hybrid II (flat analysis prior)", "B4|HybridII|n_star", "n*", "T3-05")
addw("B4std", "Supplement B.4, Standard", "B4|Standard|n_star", "n*", "T3-06")
addw("B4std", "Supplement B.4, Hybrid I", "B4|HybridI|n_star", "n*", "T3-06")
# Section 6.1
addw("S610", "Section 6.1, rho = 0.10", "S61|rho0.10|n_star", "n*", "T3-07")
addw("S605", "Section 6.1, rho = 0.05", "S61|rho0.05|n_star", "n*", "T3-08",
     note = "The original code (Data S2) returned no result for this case (empty cell; NaN from N = 1030, ledger T3-17)")
addw("S601", "Section 6.1, rho = 0.01", "S61|rho0.01|n_star", "n*", "T3-09",
     note = "Not reached by N = 4000 with the package")
addw("S601", "Section 6.1, rho = 0.01", "S61|rho0.01|max_A_to_4000", "max A over N = 2 to 4000 (package); A(4000) (programmed rule)", "T3-09")
addw("S601", "Section 6.1, rho = 0.01", "S61|rho0.01|A_N2770", "A(2770)", "T3-09")
# Section 5.3
CM <- "Prior-data conflict; P(X <= x_obs) x 100"
addw("CONFprev", "Section 5.3, prevalence (53 of 150)", "CONF|prev_pctl_lower", "percentile", "T3-10;T3-11", mode = CM)
addw("CONFprev", "Section 5.3, prevalence (53 of 150)", "CONF|prev_tail_upper", "upper tail P(X >= x_obs)", "T3-11", mode = CM)
addw("CONFse", "Section 5.3, Se (51 of 53)", "CONF|se_pctl_lower", "percentile", "T3-10;T3-11", mode = CM)
addw("CONFse", "Section 5.3, Se (51 of 53)", "CONF|se_tail_upper", "upper tail P(X >= x_obs)", "T3-11", mode = CM)
# Table 2 of Wilson et al. (grid)
T2M <- "Se only (marginal); two-sided 95%; full width; A(106) and n* on the 63-direction grid of Table 2 of Wilson et al."
for (q in c("Se", "prev")) for (s in c("A_min", "A_max", "nstar_at_Amin", "nstar_at_Amax")) {
  lab <- c(A_min = "A(106) minimum", A_max = "A(106) maximum", nstar_at_Amin = "n* at the minimum of A",
           nstar_at_Amax = "n* at the maximum of A")[[s]]
  addw(paste0("T2", q), paste0("Table 2, ", if (q == "Se") "Se" else "prevalence"), sprintf("T2|%s|%s", q, s),
       lab, "T3-12;T3-13;T3-14", mode = T2M)
}
# Williamson et al.
WM <- "Joint (Se and Sp); one-sided 95%; half-widths 0.10 (Se) and 0.05 (Sp); assurance 0.80; exact enumeration"
r321 <- one(f32a, f32a$N == 321); r320 <- one(f32a, f32a$N == 320)
add3("Williamson et al. 2023", "WILL", "Illustrative example", "N", WM,
     "Se Beta(25, 7), Sp Beta(30, 2), prevalence Beta(13.56, 122.06)", num(r321$published_N),
     r321$source, num(r321$package_N_first), ledger = "T3-19")
add3("Williamson et al. 2023", "WILL", "Illustrative example", "A(N)", WM, "as above", NA_real_, "",
     num(r321$A_package), ledger = "T3-19")
add3("Williamson et al. 2023", "WILL", "Illustrative example", "A(N - 1)", WM, "as above", NA_real_, "",
     num(r320$A_package), ledger = "T3-19")
FM <- "Se only; one-sided; half-width 0.10; n_c = 10; assurance 0.80; exact enumeration"
for (rho in c(0.3, 0.05)) for (pr in c("Beta(9,3)", "Beta(10,2)", "Beta(11,1)")) {
  r <- one(f32b, f32b$rho_T == rho & f32b$prior_se == pr)
  vv <- f32v[f32v$rho_T == rho & f32v$variant == "V2a", ]
  vcol <- c("Beta(9,3)" = "N_Se_Beta_9_3", "Beta(10,2)" = "N_Se_Beta_10_2", "Beta(11,1)" = "N_Se_Beta_11_1")[[pr]]
  nt <- if (rho == 0.05) sprintf("Prevalence prior derived without rounding, %s; with the prevalence prior rounded to 2 decimals (%s) the package gives %s (variant V2a)",
                                 r$prev_prior, vv$prev_prior, vv[[vcol]]) else ""
  add3("Williamson et al. 2023", sprintf("F3_%s", rho), sprintf("Figure 3, rho_T = %s, Se prior %s", rho, pr),
       "N", FM, sprintf("prevalence prior %s", r$prev_prior), num(r$published_N), r$published_source,
       num(r$package_N), ledger = if (rho == 0.3) "T3-21" else "T3-22;T3-23", note = nt)
}
# Web application
AM <- "Joint (Se and Sp); half-widths; assurance 0.80; web application (Calculator tab, Both)"
app_add <- function(group, case, runs, ledger) {
  r <- f33[f33$run %in% runs, , drop = FALSE]
  stopifnot(nrow(r) == length(runs), length(unique(r$app_N)) == 1, length(unique(r$app_A)) == 1,
            length(unique(r$package_N)) == 1)
  inp <- sprintf("Se %s, Sp %s, prevalence %s; half-widths %s and %s", r$Se_prior[1], r$Sp_prior[1],
                 r$prev_prior[1], r$halfwidth_Se[1], r$halfwidth_Sp[1])
  src <- sprintf("Captures %s (runs %s)", paste(unique(r$capture), collapse = "; "), paste(r$run, collapse = ", "))
  add3("Web application", group, case, "N", AM, inp, r$app_N[1], src, r$package_N[1], ledger = ledger)
  add3("Web application", group, case, "A(N)", AM, inp, r$app_A[1], src,
       r$package_A_at_N[1], ledger = ledger)
}
app_add("APP_A", "Configuration A (illustrative example)", "A-1", "T3-24")
app_add("APP_B", "Configuration B (stable runs)", c("B-3", "B-4"), "T3-25")
app_add("APP_C", "Configuration C (stable runs)", c("C-2", "C-3", "C-4", "C-5"), "T3-26")

T3 <- do.call(rbind, L3)
isdiff <- T3$quantity == "n*" & T3$group %in% c("VAP", "S610", "S605")
T3$published_minus_exact <- ifelse(isdiff, T3$published - T3$dtasamplesize_exact, NA_real_)
T3$relative_to_exact <- ifelse(isdiff, (T3$published - T3$dtasamplesize_exact) / T3$dtasamplesize_exact, NA_real_)
cat("published minus exact (n*) and relative to exact:\n")
print(T3[isdiff, c("case", "published", "dtasamplesize_exact", "published_minus_exact", "relative_to_exact")], row.names = FALSE)
T3 <- T3[, c("block", "group", "case", "quantity", "mode_exercised", "inputs", "published", "published_source",
             "dtasamplesize_exact", "programmed_rule", "programmed_rule_dataS2", "published_minus_exact",
             "relative_to_exact", "ledger_id", "note")]
g3 <- function(group, quantity, col, case = NULL) {
  cond <- T3$group == group & T3$quantity == quantity
  if (!is.null(case)) cond <- cond & T3$case == case
  r <- T3[cond, , drop = FALSE]
  if (nrow(r) != 1L) stop("g3: ", group, "/", quantity, " gave ", nrow(r), " rows")
  r[[col]]
}
g3c <- function(group, case, quantity, col) g3(group, quantity, col, case)

# --- anchors, Wilson et al.
anchor("T3-02", "VAP published n*, A(150)", c(g3("VAP", "n*", "published"), g3("VAP", "A(150)", "published")), c(106, 0.88), dec = 2)
anchor("T3-03", "VAP package n*, A(104), A(106), A(150)",
       c(g3("VAP", "n*", "dtasamplesize_exact"), g3("VAP", "A(104)", "dtasamplesize_exact"),
         g3("VAP", "A(106)", "dtasamplesize_exact"), g3("VAP", "A(150)", "dtasamplesize_exact")),
       c(104, 0.8005, 0.8049, 0.8826), dec = 4)
anchor("T3-04", "VAP programmed rule: emulation n*, Data S2 n*, emulation A(150)",
       c(g3("VAP", "n*", "programmed_rule"), g3("VAP", "n*", "programmed_rule_dataS2"), g3("VAP", "A(150)", "programmed_rule")),
       c(106, 106, 0.8825), dec = 4)
for (grp in c("B4flat", "B4std")) {
  rr <- T3[T3$group == grp, ]
  anchor(if (grp == "B4flat") "T3-05" else "T3-06", sprintf("B.4 %s: published; package; programmed rule (both keys)", grp),
         c(rr$published, rr$dtasamplesize_exact, rr$programmed_rule),
         if (grp == "B4flat") rep(290, 6) else c(106, 106, 104, 104, 106, 106))
}
anchor("T3-07", "6.1 rho 0.10: published; package; rule", c(g3("S610", "n*", "published"), g3("S610", "n*", "dtasamplesize_exact"), g3("S610", "n*", "programmed_rule")), c(681, 661, 681))
anchor("T3-08", "6.1 rho 0.05: published; package; rule", c(g3("S605", "n*", "published"), g3("S605", "n*", "dtasamplesize_exact"), g3("S605", "n*", "programmed_rule")), c(1643, 1408, 1643))
anchor("T3-08b", "6.1 rho 0.05: Data S2 column empty", g3("S605", "n*", "programmed_rule_dataS2"), NA_real_)
anchor("T3-09", "6.1 rho 0.01: published 2770; package n* not reached (NA); max A; A(2770)",
       c(g3("S601", "n*", "published"), g3("S601", "n*", "dtasamplesize_exact"),
         g3("S601", "max A over N = 2 to 4000 (package); A(4000) (programmed rule)", "dtasamplesize_exact"),
         g3("S601", "A(2770)", "dtasamplesize_exact")), c(2770, NA, 0.6917, 0.6465), dec = 4)
anchor("T3-10", "5.3 published percentiles (prevalence; Se)", c(g3("CONFprev", "percentile", "published"), g3("CONFse", "percentile", "published")), c(99, 76))
anchor("T3-11", "5.3 package percentiles x 100 (prevalence; Se)", c(g3("CONFprev", "percentile", "dtasamplesize_exact"), g3("CONFse", "percentile", "dtasamplesize_exact")), c(99.1, 75.9), dec = 1)
anchor("T3-11b", "5.3 package upper tails (prevalence; Se)", c(g3("CONFprev", "upper tail P(X >= x_obs)", "dtasamplesize_exact"), g3("CONFse", "upper tail P(X >= x_obs)", "dtasamplesize_exact")), c(0.012, 0.390), dec = 3)
pub12 <- c(g3("T2Se", "A(106) minimum", "published"), g3("T2Se", "A(106) maximum", "published"),
           g3("T2Se", "n* at the minimum of A", "published"), g3("T2Se", "n* at the maximum of A", "published"),
           g3("T2prev", "A(106) minimum", "published"), g3("T2prev", "A(106) maximum", "published"),
           g3("T2prev", "n* at the minimum of A", "published"), g3("T2prev", "n* at the maximum of A", "published"))
anchor("T3-12", "Table 2 published: Se A range, n*; prevalence A range, n*", pub12, c(0.73, 0.86, 130, 82, 0.80, 0.81, 108, 104), dec = 2)
rule8 <- c(g3("T2Se", "A(106) minimum", "programmed_rule"), g3("T2Se", "A(106) maximum", "programmed_rule"),
           g3("T2Se", "n* at the minimum of A", "programmed_rule"), g3("T2Se", "n* at the maximum of A", "programmed_rule"),
           g3("T2prev", "A(106) minimum", "programmed_rule"), g3("T2prev", "A(106) maximum", "programmed_rule"),
           g3("T2prev", "n* at the minimum of A", "programmed_rule"), g3("T2prev", "n* at the maximum of A", "programmed_rule"))
isA <- c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE)
match8 <- ifelse(isA, round(rule8, 2) == pub12, rule8 == pub12)
anchor("T3-13", "programmed rule reproduces the eight published values of Table 2 (A to 2 decimals, n* exact)", sum(match8), 8L)
pkg14 <- c(g3("T2Se", "A(106) minimum", "dtasamplesize_exact"), g3("T2Se", "A(106) maximum", "dtasamplesize_exact"),
           g3("T2Se", "n* at the minimum of A", "dtasamplesize_exact"), g3("T2Se", "n* at the maximum of A", "dtasamplesize_exact"),
           g3("T2prev", "A(106) minimum", "dtasamplesize_exact"), g3("T2prev", "A(106) maximum", "dtasamplesize_exact"),
           g3("T2prev", "n* at the minimum of A", "dtasamplesize_exact"), g3("T2prev", "n* at the maximum of A", "dtasamplesize_exact"))
anchor("T3-14", "Table 2 package on the authors' grid: Se A min-max, n*; prevalence", pkg14, c(0.741, 0.861, 129, 76, 0.800, 0.809, 106, 102), dec = 3)
# contour of 629 directions (S3.1)
anchor("T3-15", "exact contour (629 directions): Se A range, n* range; prevalence A range, n* range",
       c(wcol("T2|Se|contour629|A_min", "package"), wcol("T2|Se|contour629|A_max", "package"),
         wcol("T2|Se|contour629|nstar_min_over_contour", "package"), wcol("T2|Se|contour629|nstar_max_over_contour", "package"),
         wcol("T2|prev|contour629|A_min", "package"), wcol("T2|prev|contour629|A_max", "package"),
         wcol("T2|prev|contour629|nstar_min_over_contour", "package"), wcol("T2|prev|contour629|nstar_max_over_contour", "package")),
       c(0.742, 0.871, 73, 129, 0.800, 0.809, 102, 106), dec = 3)
# differences (DIS-04)
dd <- c(g3("VAP", "n*", "published_minus_exact"), g3("S610", "n*", "published_minus_exact"), g3("S605", "n*", "published_minus_exact"))
rr <- c(g3("VAP", "n*", "relative_to_exact"), g3("S610", "n*", "relative_to_exact"), g3("S605", "n*", "relative_to_exact"))
anchor("DIS-04a", "published minus exact n* (VAP; rho 0.10; rho 0.05)", dd, c(2, 20, 235))
anchor("DIS-04b", "relative to exact, per cent to 1 decimal", round(100 * rr, 1), c(1.9, 3.0, 16.7), tol = 1e-9)
# --- anchors, Williamson et al.
anchor("T3-19", "illustrative: published N; package N, A(N), A(N - 1)",
       c(g3("WILL", "N", "published"), g3("WILL", "N", "dtasamplesize_exact"), g3("WILL", "A(N)", "dtasamplesize_exact"),
         g3("WILL", "A(N - 1)", "dtasamplesize_exact")), c(321, 321, 0.8008, 0.7986), dec = 4)
f3_03 <- T3[T3$group == "F3_0.3", ]; f3_005 <- T3[T3$group == "F3_0.05", ]
anchor("T3-21", "Figure 3, rho_T = 0.30: published; package", c(f3_03$published, f3_03$dtasamplesize_exact), c(182, 154, 101, 182, 154, 101))
anchor("T3-22", "Figure 3, rho_T = 0.05: published; package (unrounded)", c(f3_005$published, f3_005$dtasamplesize_exact), c(1630, 1295, 783, 1628, 1294, 782))
vv05 <- f32v[f32v$rho_T == 0.05 & f32v$variant == "V2a", ]
anchor("T3-23a", "V2a (prevalence prior rounded to 2 decimals) reproduces the three published values", c(vv05$N_Se_Beta_9_3, vv05$N_Se_Beta_10_2, vv05$N_Se_Beta_11_1), c(1630, 1295, 783))
pp <- function(s) as.numeric(strsplit(gsub("[Beta() ]", "", s), ",")[[1]])
anchor("T3-23b", "derived prevalence priors to 3 decimals (rho_T 0.05; 0.30)",
       c(round(pp(sub("prevalence prior ", "", f3_005$inputs[1])), 3), round(pp(sub("prevalence prior ", "", f3_03$inputs[1])), 3)),
       c(3.433, 65.229, 96.337, 224.786), tol = 1e-9)
anchor("T3-23c", "rounded prior of the variant V2a", pp(vv05$prev_prior), c(3.43, 65.23), tol = 1e-9)
# --- anchors, application
a_A <- T3[T3$group == "APP_A", ]; a_B <- T3[T3$group == "APP_B", ]; a_C <- T3[T3$group == "APP_C", ]
anchor("T3-24", "application A: app N, A; package N, A", c(a_A$published, a_A$dtasamplesize_exact), c(321, 0.8008, 321, 0.8008), dec = 4)
anchor("T3-25", "application B stable: app N, A; package N, A", c(a_B$published, a_B$dtasamplesize_exact), c(551, 0.8005, 551, 0.8005), dec = 4)
anchor("T3-26", "application C stable: app N, A; package N, A", c(a_C$published, a_C$dtasamplesize_exact), c(552, 0.8007, 552, 0.8007), dec = 4)
anchor("T3-24n", "application A, B, C: package N exact", c(a_A$dtasamplesize_exact[1], a_B$dtasamplesize_exact[1], a_C$dtasamplesize_exact[1]), c(321, 551, 552))

T3d <- local({
  d <- list()
  addd <- function(block, case, mode, inputs, pub, src, exact, rule = "", ledger, note = "") {
    d[[length(d) + 1L]] <<- data.frame(Block = block, Case = case, `Mode exercised` = mode, Inputs = inputs,
      Published = pub, `Published source` = src, `dtasamplesize (exact)` = exact, `Programmed rule` = rule,
      `Published minus exact` = "", `Relative to exact` = "", `Ledger ID` = ledger, Note = note,
      check.names = FALSE, stringsAsFactors = FALSE)
  }
  W <- "Wilson et al. 2022"
  src <- function(group, quantity) T3$published_source[T3$group == group & T3$quantity == quantity][1]
  inp <- function(group) T3$inputs[T3$group == group][1]
  addd(W, "VAP, n*", WMODE, inp("VAP"), fi(g3("VAP", "n*", "published")), src("VAP", "n*"),
       fi(g3("VAP", "n*", "dtasamplesize_exact")),
       sprintf("%s (Data S2: %s)", fi(g3("VAP", "n*", "programmed_rule")), fi(g3("VAP", "n*", "programmed_rule_dataS2"))),
       "T3-02;T3-03;T3-04")
  addd(W, "VAP, A(150)", WMODE, inp("VAP"), fp(g3("VAP", "A(150)", "published"), 2), src("VAP", "A(150)"),
       fp(g3("VAP", "A(150)", "dtasamplesize_exact")), fp(g3("VAP", "A(150)", "programmed_rule")), "T3-03;T3-04")
  addd(W, "B.4, flat analysis prior (Combination; Hybrid II)", WMODE, inp("B4flat"),
       fi(g3("B4flat", "n*", "published", "Supplement B.4, Combination (flat analysis prior)")),
       T3$published_source[T3$group == "B4flat"][1], fi(g3("B4flat", "n*", "dtasamplesize_exact", "Supplement B.4, Combination (flat analysis prior)")),
       fi(g3("B4flat", "n*", "programmed_rule", "Supplement B.4, Combination (flat analysis prior)")), "T3-05")
  addd(W, "B.4, Standard; Hybrid I", WMODE, inp("B4std"),
       fi(g3("B4std", "n*", "published", "Supplement B.4, Standard")), T3$published_source[T3$group == "B4std"][1],
       fi(g3("B4std", "n*", "dtasamplesize_exact", "Supplement B.4, Standard")),
       fi(g3("B4std", "n*", "programmed_rule", "Supplement B.4, Standard")), "T3-06")
  addd(W, "Section 6.1, rho = 0.10", WMODE, inp("S610"), fi(g3("S610", "n*", "published")), src("S610", "n*"),
       fi(g3("S610", "n*", "dtasamplesize_exact")), fi(g3("S610", "n*", "programmed_rule")), "T3-07")
  addd(W, "Section 6.1, rho = 0.05", WMODE, inp("S605"), fi(g3("S605", "n*", "published")), src("S605", "n*"),
       fi(g3("S605", "n*", "dtasamplesize_exact")), fi(g3("S605", "n*", "programmed_rule")), "T3-08",
       "Data S2 returned no result for this case")
  addd(W, "Section 6.1, rho = 0.01", WMODE, inp("S601"), fi(g3("S601", "n*", "published")), src("S601", "n*"),
       sprintf("not reached by N = 4000 (maximum A %s; A(2770) %s)",
               fp(g3("S601", "max A over N = 2 to 4000 (package); A(4000) (programmed rule)", "dtasamplesize_exact")),
               fp(g3("S601", "A(2770)", "dtasamplesize_exact"))),
       sprintf("not reached by N = 4000 (A(4000) %s; A(2770) %s)",
               fp(g3("S601", "max A over N = 2 to 4000 (package); A(4000) (programmed rule)", "programmed_rule")),
               fp(g3("S601", "A(2770)", "programmed_rule"))), "T3-09")
  addd(W, "Section 5.3, conflict percentiles (prevalence; Se)", CM, "VAP priors; 53 of 150 with the disease; 51 of 53 with a positive test",
       sprintf("%s; %s", fi(g3("CONFprev", "percentile", "published")), fi(g3("CONFse", "percentile", "published"))),
       src("CONFprev", "percentile"),
       sprintf("%s; %s (upper tails %s; %s)", fp(g3("CONFprev", "percentile", "dtasamplesize_exact"), 1),
               fp(g3("CONFse", "percentile", "dtasamplesize_exact"), 1),
               fp(g3("CONFprev", "upper tail P(X >= x_obs)", "dtasamplesize_exact"), 3),
               fp(g3("CONFse", "upper tail P(X >= x_obs)", "dtasamplesize_exact"), 3)), "", "T3-10;T3-11")
  for (q in c("Se", "prev")) {
    gq <- paste0("T2", q)
    addd(W, sprintf("Table 2, %s", if (q == "Se") "Se" else "prevalence"), T2M, inp(gq),
         sprintf("A(106) %s to %s; n* %s to %s", fp(g3(gq, "A(106) minimum", "published"), 2), fp(g3(gq, "A(106) maximum", "published"), 2),
                 fi(g3(gq, "n* at the maximum of A", "published")), fi(g3(gq, "n* at the minimum of A", "published"))),
         src(gq, "A(106) minimum"),
         sprintf("A(106) %s to %s; n* %s at the maximum and %s at the minimum of A", fp(g3(gq, "A(106) minimum", "dtasamplesize_exact"), 3),
                 fp(g3(gq, "A(106) maximum", "dtasamplesize_exact"), 3), fi(g3(gq, "n* at the maximum of A", "dtasamplesize_exact")),
                 fi(g3(gq, "n* at the minimum of A", "dtasamplesize_exact"))),
         sprintf("A(106) %s to %s; n* %s and %s", fp(g3(gq, "A(106) minimum", "programmed_rule"), 3), fp(g3(gq, "A(106) maximum", "programmed_rule"), 3),
                 fi(g3(gq, "n* at the maximum of A", "programmed_rule")), fi(g3(gq, "n* at the minimum of A", "programmed_rule"))),
         "T3-12;T3-13;T3-14")
  }
  Wm <- "Williamson et al. 2023"
  addd(Wm, "Illustrative example", WM, T3$inputs[T3$group == "WILL"][1], fi(g3("WILL", "N", "published")), T3$published_source[T3$group == "WILL"][1],
       sprintf("%s (A = %s; A(%s) = %s)", fi(g3("WILL", "N", "dtasamplesize_exact")), fp(g3("WILL", "A(N)", "dtasamplesize_exact")),
               fi(g3("WILL", "N", "dtasamplesize_exact") - 1), fp(g3("WILL", "A(N - 1)", "dtasamplesize_exact"))), "", "T3-19")
  addd(Wm, "Figure 3, rho_T = 0.30 (Se priors Beta(9, 3), Beta(10, 2), Beta(11, 1))", FM, f3_03$inputs[1],
       paste(fi(f3_03$published), collapse = ", "), f3_03$published_source[1], paste(fi(f3_03$dtasamplesize_exact), collapse = ", "), "", "T3-21")
  addd(Wm, "Figure 3, rho_T = 0.05 (Se priors Beta(9, 3), Beta(10, 2), Beta(11, 1))", FM, f3_005$inputs[1],
       paste(fi(f3_005$published), collapse = ", "), f3_005$published_source[1],
       sprintf("%s; with the prevalence prior rounded to 2 decimals: %s", paste(fi(f3_005$dtasamplesize_exact), collapse = ", "),
               paste(fi(c(vv05$N_Se_Beta_9_3, vv05$N_Se_Beta_10_2, vv05$N_Se_Beta_11_1)), collapse = ", ")), "", "T3-22;T3-23")
  for (gg in c("APP_A", "APP_B", "APP_C")) {
    z <- T3[T3$group == gg, ]
    addd("Web application", z$case[1], AM, z$inputs[1], sprintf("%s (A = %s)", fi(z$published[1]), fp(z$published[2])), z$published_source[1],
         sprintf("%s (A = %s)", fi(z$dtasamplesize_exact[1]), fp(z$dtasamplesize_exact[2])), "", z$ledger_id[1])
  }
  out <- do.call(rbind, d)
  # published minus exact (CF-06) for the three rows that have it
  fill <- function(case, grp) {
    i <- which(out$Case == case)
    out[i, "Published minus exact"] <<- fi(g3(grp, "n*", "published_minus_exact"))
    out[i, "Relative to exact"] <<- sprintf("%.1f%%", 100 * g3(grp, "n*", "relative_to_exact"))
  }
  fill("VAP, n*", "VAP"); fill("Section 6.1, rho = 0.10", "S610"); fill("Section 6.1, rho = 0.05", "S605")
  out
})
# The Inputs cells come from compare_wilson2022.R, which is not re-run here; the article states widths as full widths.
T3$inputs <- gsub("total width", "full width", T3$inputs, fixed = TRUE)
T3d$Inputs <- gsub("total width", "full width", T3d$Inputs, fixed = TRUE)
emit_main("R1_Table3_reproductions", T3, T3d)

# =============================================================================
# Table 5: R1_Table5_applied_example
# =============================================================================
cat("\n== Table 5: applied example ==\n")
am <- rd("applied", "applied_example_main_table.csv")
as_ <- rd("applied", "applied_example_sensitivity.csv")
af  <- rd("applied", "applied_example_frequentist.csv")
apc <- rd("applied", "applied_example_published_comparison.csv")
# Pilot counts (Albert et al. 2021, Abstract and Results)
PIL <- c(TP = 43, FN = 11, FP = 0, TN = 358)
anchor("T5-01", "pilot counts and prevalence 54/412", c(PIL, sum(PIL), sum(PIL[c("TP", "FN")]), round(100 * sum(PIL[c("TP", "FN")]) / sum(PIL), 1)),
       c(43, 11, 0, 358, 412, 54, 13.1), tol = 1e-9)
A0 <- 0.5
dp <- c(se_a = 1 + A0 * PIL[["TP"]], se_b = 1 + A0 * PIL[["FN"]], sp_a = 1 + A0 * PIL[["TN"]], sp_b = 1 + A0 * PIL[["FP"]],
        pv_a = 1 + A0 * (PIL[["TP"]] + PIL[["FN"]]), pv_b = 1 + A0 * (sum(PIL) - PIL[["TP"]] - PIL[["FN"]]))
mp1 <- one(am, am$scenario == "design_onesided_joint")
csvp <- c(mp1$prior_se_a, mp1$prior_se_b, mp1$prior_sp_a, mp1$prior_sp_b, mp1$prior_prev_a, mp1$prior_prev_b)
anchor("T5-05", "design priors read from CSV (Se, Sp, prevalence) = 1 + 0.5 x counts", csvp, c(22.5, 6.5, 180, 1, 28, 180), tol = 1e-12)
anchor("T5-05b", "arithmetic from the pilot counts equals the CSV", unname(dp), csvp, tol = 1e-12)
mrow <- function(s) one(am, am$scenario == s)
srow <- function(s) one(as_, as_$scenario == s)
scn_main <- c("design_onesided_joint", "both_onesided_joint", "design_twosided_joint", "both_twosided_joint", "design_onesided_se", "design_twosided_se")
anchor("T5-06", "principal N, A, A(N - 1)", c(mrow("design_onesided_joint")$N_total, mrow("design_onesided_joint")$A_N, mrow("design_onesided_joint")$A_N_minus_1), c(497, 0.8007, 0.7985), dec = 4)
anchor("T5-07", "both one-sided joint N, A", c(mrow("both_onesided_joint")$N_total, mrow("both_onesided_joint")$A_N), c(258, 0.8012), dec = 4)
anchor("T5-08", "two-sided joint N (design; both)", c(mrow("design_twosided_joint")$N_total, mrow("both_twosided_joint")$N_total), c(612, 363))
anchor("T5-09", "Se-only N (one-sided; two-sided)", c(mrow("design_onesided_se")$N_total, mrow("design_twosided_se")$N_total), c(497, 612))
sel <- c("a0_0.25_design", "a0_1_design", "dse_0.075", "dse_0.125", "target_0.7", "target_0.9")
anchor("T5-10", "a0 = 0.25 and 1 (design) N", c(srow("a0_0.25_design")$N_total, srow("a0_1_design")$N_total), c(520, 480))
anchor("T5-12", "Se half-width 0.075 and 0.125 N", c(srow("dse_0.075")$N_total, srow("dse_0.125")$N_total), c(859, 325))
anchor("T5-14", "assurance target 0.70 and 0.90 N", c(srow("target_0.7")$N_total, srow("target_0.9")$N_total), c(457, 558))
fq1 <- one(af, startsWith(af$case, "one-sided 95%, Sp = Wilson lower limit (principal)"))
fq4 <- one(af, startsWith(af$case, "two-sided 95%, half-widths 0.10/0.05, Sp = Wilson lower limit"))
anchor("T5-17", "Buderer formula N with fixed values (one-sided; two-sided)", c(fq1$N_total, fq4$N_total), c(336, 481))
pubN <- function(s) apc$N[startsWith(apc$source, s)]
anchor("T5-03", "published frequentist plan", pubN("Bulilete 2021 Methods"), 927)
anchor("T5-02", "enrolled; analysed", c(pubN("Bulilete 2021 enrolled"), pubN("Bulilete 2021 analysed")), c(1369, 1362))
anchor("T5-19", "principal N as a percentage of 927", round(100 * mrow("design_onesided_joint")$N_total / pubN("Bulilete 2021 Methods")), 54)

T5 <- list()
add5 <- function(block, label, criterion = "", objective = "", analysis_prior = "", value = "", N_total = NA, A_N = NA, A_N_minus_1 = NA, ledger = "") {
  T5[[length(T5) + 1L]] <<- data.frame(block = block, label = label, criterion = criterion, objective = objective,
    analysis_prior = analysis_prior, value = value, N_total = N_total, A_N = A_N, A_N_minus_1 = A_N_minus_1,
    ledger_id = ledger, stringsAsFactors = FALSE)
}
add5("A_inputs", "Pilot counts (TP, FN, FP, TN)", value = sprintf("%d, %d, %d, %d (N = %d; prevalence %d/%d = %.1f%%)", PIL[["TP"]], PIL[["FN"]], PIL[["FP"]], PIL[["TN"]], sum(PIL), PIL[["TP"]] + PIL[["FN"]], sum(PIL), 100 * 54 / 412), ledger = "T5-01")
add5("A_inputs", "Base priors and power-prior weight", value = "Beta(1, 1) for Se, Sp and prevalence; a0 = 0.5", ledger = "T5-04")
add5("A_inputs", "Design prior, Se", value = sprintf("Beta(%s, %s)", format(mp1$prior_se_a), format(mp1$prior_se_b)), ledger = "T5-05")
add5("A_inputs", "Design prior, Sp", value = sprintf("Beta(%s, %s)", format(mp1$prior_sp_a), format(mp1$prior_sp_b)), ledger = "T5-05")
add5("A_inputs", "Design prior, prevalence", value = sprintf("Beta(%s, %s)", format(mp1$prior_prev_a), format(mp1$prior_prev_b)), ledger = "T5-05")
for (s in scn_main) {
  r <- mrow(s)
  tw <- r$width_type == "full_width"
  se_only <- r$assurance_for == "se"
  add5("B_main", s, criterion = if (tw) {
         if (se_only) sprintf("two-sided 95%%, full width %s (Se)", format(r$delta_se)) else
           sprintf("two-sided 95%%, full widths %s and %s", format(r$delta_se), format(r$delta_sp))
       } else {
         if (se_only) sprintf("one-sided 95%%, lower half-width %s (Se)", format(r$delta_se)) else
           sprintf("one-sided 95%%, lower half-widths %s and %s", format(r$delta_se), format(r$delta_sp))
       },
       objective = if (r$assurance_for == "joint") "joint (Se and Sp)" else "marginal (Se only)",
       analysis_prior = if (r$analysis_prior_se_a == 1 && r$analysis_prior_sp_a == 1) "flat Beta(1, 1) (use_in = \"design\")" else "updated pilot prior (use_in = \"both\")",
       N_total = r$N_total, A_N = r$A_N, A_N_minus_1 = r$A_N_minus_1,
       ledger = c(design_onesided_joint = "T5-06", both_onesided_joint = "T5-07", design_twosided_joint = "T5-08",
                  both_twosided_joint = "T5-08", design_onesided_se = "T5-09", design_twosided_se = "T5-09")[[s]])
}
for (s in sel) {
  r <- srow(s)
  add5("C_sensitivity", s, criterion = sprintf("one-sided 95%%, lower half-widths %s and %s; target %s", format(r$delta_se), format(r$delta_sp), format(r$target_assurance)),
       objective = "joint (Se and Sp)", analysis_prior = "flat Beta(1, 1) (use_in = \"design\")",
       N_total = r$N_total, A_N = r$A_N, A_N_minus_1 = r$A_N_minus_1,
       ledger = if (startsWith(s, "a0")) "T5-10" else if (startsWith(s, "dse")) "T5-12" else "T5-14")
}
add5("D_comparison", "Buderer's formula (expected plug-in width), fixed values, one-sided 95% (Se 0.796, prevalence 0.131, Sp = Wilson lower limit 0.9894)",
     criterion = "one-sided 95%, half-widths 0.10 and 0.05", objective = "Se and Sp arms", N_total = fq1$N_total, ledger = "T5-17")
add5("D_comparison", "Buderer's formula (expected plug-in width), fixed values, two-sided 95% (same inputs)",
     criterion = "two-sided 95%, half-widths 0.10 and 0.05", objective = "Se and Sp arms", N_total = fq4$N_total, ledger = "T5-17")
add5("D_comparison", "Planned by Bulilete et al. (frequentist)", criterion = "prevalence 15%, margin 5%, Se 90%", N_total = pubN("Bulilete 2021 Methods"), ledger = "T5-03")
add5("D_comparison", "Enrolled by Bulilete et al.", N_total = pubN("Bulilete 2021 enrolled"), ledger = "T5-02")
add5("D_comparison", "Analysed by Bulilete et al.", N_total = pubN("Bulilete 2021 analysed"), ledger = "T5-02")
T5f <- do.call(rbind, T5)
T5d <- data.frame(Block = T5f$block, Label = T5f$label, Criterion = T5f$criterion, Objective = T5f$objective,
                  `Analysis prior` = T5f$analysis_prior, Value = T5f$value, `N total` = fi(T5f$N_total),
                  `A(N)` = fp(T5f$A_N), `A(N - 1)` = fp(T5f$A_N_minus_1), `Ledger ID` = T5f$ledger_id,
                  check.names = FALSE, stringsAsFactors = FALSE)
emit_main("R1_Table5_applied_example", T5f, T5d)

# =============================================================================
# Table 6: R1_Table6_operating_characteristics
# =============================================================================
cat("\n== Table 6: operating characteristics ==\n")
ocp <- rd("oc", "oc_summary_planning_point.csv")
ocd <- rd("oc", "oc_summary_by_design.csv")
mps <- rd("mp", "mp_summary.csv")
DES <- c("Buderer_N500_Wald", "Joint_N580_Wilson", "BAM_P1_N678", "BAM_P2_N672")
DLAB <- c("Buderer, N = 500, Wald interval", "Joint Wilson (joint_sample_size), N = 580",
          "BAM, vague Sp prior (P1), N = 678", "BAM, informative Sp prior (P2), N = 672")
mpv <- function(d, e, col) { r <- one(mps, mps$design == d & mps$scenario == e); num(r[[col]]) }
T6 <- data.frame(design = DES, design_label = DLAB, stringsAsFactors = FALSE)
T6$N <- vapply(DES, function(d) one(ocp, ocp$design == d)$N, numeric(1))
T6$p_success_planning <- vapply(DES, function(d) one(ocp, ocp$design == d)$p_success, numeric(1))
T6$cov_se_planning <- vapply(DES, function(d) one(ocp, ocp$design == d)$cov_se, numeric(1))
T6$cov_sp_planning <- vapply(DES, function(d) one(ocp, ocp$design == d)$cov_sp, numeric(1))
T6$p_success_se075 <- vapply(DES, mpv, numeric(1), e = "S01", col = "p_success")
T6$p_success_se065 <- vapply(DES, mpv, numeric(1), e = "S02", col = "p_success")
T6$p_success_prev008 <- vapply(DES, mpv, numeric(1), e = "S06", col = "p_success")
T6$detect_se_se065 <- vapply(DES, mpv, numeric(1), e = "S02", col = "det_se_two_sided_analysis")
T6$detect_any_se065 <- vapply(DES, mpv, numeric(1), e = "S02", col = "det_any_two_sided_analysis")
T6$detect_any_prior_predictive <- vapply(DES, function(d) {
  r <- mps[mps$design == d & mps$scenario == "PRIOR_PREDICTIVE", ]
  if (nrow(r) == 0) NA_real_ else num(r$det_any_two_sided_analysis)
}, numeric(1))
T6$ledger_id <- "T6-01;T6-02;T6-03;T6-05;T6-06;T6-07;T6-08;T6-12"
rownames(T6) <- NULL
anchor("T6-00", "N of the four designs", T6$N, c(500, 580, 678, 672))
anchor("T6-01", "P(success) at the planning point", T6$p_success_planning, c(0.5265, 0.8070, 0.9994, 0.9992), dec = 4)
anchor("T6-02", "coverage of the Se interval", T6$cov_se_planning, c(0.9364, 0.9511, 0.9621, 0.9622), dec = 4)
anchor("T6-03", "coverage of the Sp interval", T6$cov_sp_planning, c(0.9449, 0.9504, 0.9483, 0.9528), dec = 4)
anchor("T6-05", "P(success), true Se 0.75", T6$p_success_se075, c(0.0122, 0.0666, 0.8349, 0.8166), dec = 4)
anchor("T6-06a", "P(success), true Se 0.65 (BAM)", T6$p_success_se065[3:4], c(0.1649, 0.1464), dec = 4)
anchor("T6-06b", "P(success), true Se 0.65 (Buderer; joint) in 1 x 10^-n form", c(T6$p_success_se065[1] / 1e-6, T6$p_success_se065[2] / 1e-4), c(8.0, 1.5), dec = 1)
anchor("T6-07", "P(success), true prevalence 0.08", T6$p_success_prev008, c(0.0565, 0.0224, 0.1446, 0.1399), dec = 4)
anchor("T6-08", "Se quantity flagged (two-sided, analysis prior) with true Se 0.65 (P1; P2)", T6$detect_se_se065[3:4], c(0.521, 0.519), dec = 3)
anchor("T6-08a", "any of the three quantities flagged with true Se 0.65 (P1; P2)", T6$detect_any_se065[3:4], c(0.649, 0.519), dec = 3)
anchor("T6-08s", "Sp quantity flagged at the planning point and with true Se 0.65 (P1, S00; P1, S02)",
       c(mpv("BAM_P1_N678", "S00", "det_sp_two_sided_analysis"), mpv("BAM_P1_N678", "S02", "det_sp_two_sided_analysis")), c(0.268, 0.268), dec = 3)
anchor("T6-08n", "no detection value for the frequentist designs", T6$detect_se_se065[1:2], c(NA_real_, NA_real_))
anchor("T6-12", "prior predictive: any of three (P1; P2)", T6$detect_any_prior_predictive[3:4], c(0.133, 0.129), dec = 3)
anchor("T6-13", "Sp-width attainment at the planning point equals 1 (four designs)", vapply(DES, function(d) one(ocp, ocp$design == d)$p_width_sp, numeric(1)), rep(1, 4), tol = 1e-9)
anchor("T6-18", "share of the 504-point grid with P(success) >= 0.80 (Buderer Wald; joint; P1; P2)",
       vapply(DES, function(d) one(ocd, ocd$design == d)$frac_success_ge_080, numeric(1)), c(0.310, 0.452, 0.643, 0.643), dec = 3)
T6d <- data.frame(Design = T6$design_label, N = fi(T6$N), `P(success) at planning values` = fp(T6$p_success_planning),
                  `Coverage Se` = fp(T6$cov_se_planning), `Coverage Sp` = fp(T6$cov_sp_planning),
                  `P(success), true Se 0.75` = fp(T6$p_success_se075), `P(success), true Se 0.65` = fp(T6$p_success_se065),
                  `P(success), true prevalence 0.08` = fp(T6$p_success_prev008),
                  `Se flagged, true Se 0.65` = fp(T6$detect_se_se065),
                  `Any quantity flagged, prior predictive` = fp(T6$detect_any_prior_predictive),
                  check.names = FALSE, stringsAsFactors = FALSE)
emit_main("R1_Table6_operating_characteristics", T6, T6d)

# =============================================================================
# Additional file 3
# =============================================================================
cat("\n== Additional file 3 and Additional file 2 tables ==\n")
# --- S3.1 Wilson et al.: complete reproduction
s31 <- wt[wt$section != "d.points", c("key", "section", "case", "inputs", "quantity", "published", "published_source",
                                         "package", "emulation", "dataS2_original", "explanation")]
s31p <- wt[wt$section == "d.points", ]
anchor("S3.1a", "rows of the Wilson table: main; grid points", c(nrow(s31), nrow(s31p)), c(nrow(wt) - 126L, 126L))
emit_simple("R1_AF3_S3.1_wilson_reproduction", s31)
write_csv(s31p, "R1_AF3_S3.1_table2_grid_points.csv")
wl <- if (!is.null(named$witness)) named$witness else NA_character_
if (!is.na(wl) && file.exists(wl)) {
  cat(sprintf("INPUT --witness %s sha256=%s\n", wl, unname(tools::sha256sum(wl))))
  lines <- readLines(wl, encoding = "UTF-8", warn = FALSE)
  m <- regmatches(lines, regexec("k =\\s*(\\d+), x =\\s*(\\d+), W = ([0-9.]+), def = (True|False), code = (True|False), weight = ([-+0-9.e]+)", lines))
  m <- m[lengths(m) == 7L]
  wit <- data.frame(k = as.integer(vapply(m, `[`, "", 2)), x = as.integer(vapply(m, `[`, "", 3)),
                    width_W = as.numeric(vapply(m, `[`, "", 4)),
                    success_by_definition = vapply(m, `[`, "", 5) == "True", success_by_code_rule = vapply(m, `[`, "", 6) == "True",
                    weight = as.numeric(vapply(m, `[`, "", 7)), stringsAsFactors = FALSE)
  wit$N <- 104L; wit$source <- "wilson2022_rule_enumeration_enum.log, VAP, N = 104, ten witnesses of largest weight"
  anchor("S3.1b", "ten witness cases read from the rule-enumeration log; first witness k = 11, x = 11", c(nrow(wit), wit$k[1], wit$x[1]), c(10L, 11L, 11L))
  emit_simple("R1_AF3_S3.1_witness_cases", wit)
} else cat("NOTE: --witness not given or not found; witness table not written.\n")

# --- S3.2a, S3.2b Williamson et al.
drop_t <- function(df) df[, setdiff(names(df), "package_seconds"), drop = FALSE]
emit_simple("R1_AF3_S3.2a_williamson_illustrative", f32a)
emit_simple("R1_AF3_S3.2a_williamson_figure3", drop_t(f32b))
emit_simple("R1_AF3_S3.2a_williamson_figure3_variants", f32v)
b1 <- data.frame(config = "B", run = "B-1", max_n = 1000, app_N = NA, app_A = NA,
                 capture = "2026-11-14_2357_B_Both_maxn1000_output.png", config_typed = "B",
                 Se_prior = "Beta(17,3)", Sp_prior = "Beta(18,2)", prev_prior = "Beta(4,16)", halfwidth_Se = 0.07, halfwidth_Sp = 0.05,
                 assurance_target = 0.8, comment = "Run B-1, 2026-11-14 23:57: the application returned an error (An error has occurred); constant row transcribed from the authors' records of the application (screenshots kept by the authors, available on request)",
                 stringsAsFactors = FALSE)
f33x <- f33[, intersect(names(f33), c("config", "run", "max_n", "app_N", "app_A", "capture", "config_typed", "Se_prior", "Sp_prior",
                                      "prev_prior", "halfwidth_Se", "halfwidth_Sp", "assurance_target", "package_N", "package_A_at_N",
                                      "package_A_at_app_N", "diff_N_app_minus_package", "A_agrees_4dp", "N_agrees", "comment"))]
for (nm in setdiff(names(f33x), names(b1))) b1[[nm]] <- NA
b1 <- b1[, names(f33x)]
f33all <- rbind(f33x[f33x$run %in% c("B-2", "B-3", "B-4", "C-1", "C-2", "C-3", "C-4", "C-5"), ], b1)
f33all <- f33all[order(f33all$config, f33all$run), ]
f33all$outcome <- ifelse(f33all$run == "B-1", "error",
                   ifelse(f33all$N_agrees %in% TRUE & f33all$A_agrees_4dp %in% TRUE, "match",
                   ifelse(abs(f33all$diff_N_app_minus_package) == 1, "differs by one participant", "other")))
cnt <- table(factor(f33all$outcome, levels = c("match", "differs by one participant", "error", "other")))
cat("Runs in configurations B and C by outcome:", paste(names(cnt), cnt, sep = " = ", collapse = "; "), "\n")
anchor("T3-29", "runs B and C: match / one participant different / error (of 9)", c(cnt[["match"]], cnt[["differs by one participant"]], cnt[["error"]], nrow(f33all)), c(6, 2, 1, 9))
anchor("T3-24r", "configuration A run in the table", sum(f33$run == "A-1"), 1L)
emit_simple("R1_AF3_S3.2b_application_runs", f33all)
emit_simple("R1_AF3_S3.2b_application_all_runs_with_diagnostics", f33)
emit_simple("R1_AF3_S3.2b_application_restrictions", rd("williamson", "f33_restrictions.csv"))
write_csv(data.frame(outcome = names(cnt), runs = as.integer(cnt)), "R1_AF3_S3.2b_application_counts.csv")

# --- S3.3b matrix of capabilities of the CRAN packages
mx_path <- file.path(script_dir, "table_cran_feature_matrix_R1_source.csv")
mx <- utils::read.csv(mx_path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
cat(sprintf("INPUT (package tree) %s sha256=%s\n", mx_path, unname(tools::sha256sum(mx_path))))
prov <- rd("assets", "table3_provenance.csv")
anchor("S3.3b", "capability rows and package columns of the CRAN matrix (Table 4 of v0.6.9; Table S3.3e)", c(nrow(mx), ncol(mx)), c(7L, 8L))
pk <- sub(" \\(.*$", "", names(mx)[2:7]); ver <- sub("^.*\\((.*)\\)$", "\\1", names(mx)[2:7])
prov2 <- data.frame(package = pk, version_in_submitted_table = ver, stringsAsFactors = FALSE)
prov2 <- merge(prov2, prov, by = "package", all.x = TRUE, sort = FALSE)
# Versions whose help pages were read in November 2026 (tools::Rd_db(), R 4.6.2). The matrix source is the checked
# matrix, so the versions in its header must equal these; the August review stays in its own columns.
checked_versions <- c(presize = "0.3.11", pROC = "1.19.1", MKpower = "1.1", epiR = "3.0.1", MKmisc = "2.0", SampleSizeDiagnostics = "0.1.1")
anchor("S3.3b-v", "package versions in the matrix header equal the versions whose help pages were checked in November 2026",
       sum(prov2$version_in_submitted_table == unname(checked_versions[prov2$package])), nrow(prov2))
prov2$version_checked_2026_11 <- unname(checked_versions[prov2$package])
prov2$status <- ifelse(is.na(prov2$currently_installed_version), "no version recorded by the assets script",
                       ifelse(prov2$currently_installed_version == prov2$version_checked_2026_11,
                              "installed version equals the version checked", "installed version DIFFERS from the version checked"))
prov2$verification <- paste0("cells checked against the help pages of the installed version with tools::Rd_db() under R 4.6.2, November 2026 (version ",
                             prov2$version_checked_2026_11, ")")
emit_simple("R1_AF3_S3.3b_cran_feature_matrix", mx)
emit_simple("R1_AF3_S3.3b_cran_feature_matrix_provenance", prov2)

# --- S3.4 width attainment
ws <- rd("width", "width_attainment_summary.csv"); wc <- rd("width", "width_attainment_cells.csv"); wcv <- rd("width", "curve_Se085_d007_summary.csv")
anchor("S3.4a", "summary rows (8) and cell rows (224)", c(nrow(ws), nrow(wc)), c(8L, 224L))
anchor("S3.4b", "exact coverage at n = 100 (Wald truncated; Wilson)", c(wcv$coverage_n100_exact_comparison[wcv$interval == "wald_trunc"], wcv$coverage_n100_exact_comparison[wcv$interval == "wilson"]), c(0.9330, 0.9332), dec = 4)
emit_simple("R1_AF3_S3.4_width_attainment_summary", ws)
write_csv(wc, "R1_AF3_S3.4_width_attainment_cells.csv")
emit_simple("R1_AF3_S3.4_curve_Se085_d007_summary", wcv)

# --- S3.5 bridge
anchor("S3.5", "bridge table rows", nrow(br), 16L)
emit_simple("R1_AF3_S3.5_bridge", br)

# --- S3.6 success-set shortcut
ssc <- rd("shortcut", "scenario_summary.csv"); nca <- rd("shortcut", "noncontiguous_arm_sizes.csv"); ust <- rd("shortcut", "unimodality_stress.csv")
s36 <- ssc[ssc$rule %in% c("enumeration", "printed", "code"), ]
anchor("S3.6a", "noncontiguous arm sizes (rows); unimodality combinations; not quasi-concave", c(nrow(nca), nrow(ust), sum(ust$n_not_quasi_concave)), c(0L, 147L, 0L))
emit_simple("R1_AF3_S3.6_shortcut_scenarios", s36)
write_csv(data.frame(check = c("noncontiguous arm sizes (rows)", "unimodality combinations", "combinations not quasi-concave"),
                     value = c(nrow(nca), nrow(ust), sum(ust$n_not_quasi_concave))), "R1_AF3_S3.6_shortcut_counts.csv")

# --- S3.7 exact versus Monte Carlo
emit_simple("R1_AF3_S3.7_exact_distances", rd("exactmc", "distances_and_replicates.csv"))
emit_simple("R1_AF3_S3.7_compatible_range", rd("exactmc", "compatible_range.csv"))
emit_simple("R1_AF3_S3.7_mc_summary", rd("exactmc", "mc_summary.csv"))
emit_simple("R1_AF3_S3.7_timings", rd("exactmc", "timings.csv"), caption = "Timings depend on the machine and on the load at the time of the run.")

# --- S3.8 applied example, complete
# Table S3.8b: the blocks read "half-width", as in the article; applied_example.R is not re-run, and the Table 5 anchors read the scenario names
as_s38b <- as_
as_s38b$block <- gsub("semi-width", "half-width", as_s38b$block, fixed = TRUE)
emit_simple("R1_AF3_S3.8_sensitivity", as_s38b)
emit_simple("R1_AF3_S3.8_hellinger", rd("applied", "applied_example_hellinger_summary.csv"))
# Table S3.8d: the label of the main case reads "(main)", as in the article; applied_example.R is not re-run, and the Table 5 anchors read the label above
af_s38d <- af
af_s38d$case <- gsub("(principal)", "(main)", af_s38d$case, fixed = TRUE)
emit_simple("R1_AF3_S3.8_frequentist", af_s38d)
emit_simple("R1_AF3_S3.8_postconflict", rd("applied", "applied_example_postconflict.csv"))
emit_simple("R1_AF3_S3.8_initial_prior_check", rd("applied", "applied_example_initial_prior_check.csv"))

# --- S3.9 operating characteristics and misspecified priors
anchor("S3.9", "mp_summary rows", nrow(mps), 53L)
emit_simple("R1_AF3_S3.9_oc_by_design", ocd)
emit_simple("R1_AF3_S3.9_oc_planning_point", ocp)
emit_simple("R1_AF3_S3.9_mp_summary", mps)

# --- S3.9e and S3.9f: Monte Carlo cross-check (oc_montecarlo_crosscheck.R)
xc  <- rd("crosscheck", "oc_crosscheck_success_coverage.csv")
xd  <- rd("crosscheck", "oc_crosscheck_conflict_detection.csv")
ocw <- rd("oc", "oc_cells_wide.csv")
ocv <- function(d, se, sp, prev, col) {
  r <- ocw[ocw$design == d & abs(ocw$Se - se) < 1e-9 & abs(ocw$Sp - sp) < 1e-9 & abs(ocw$prev - prev) < 1e-9, , drop = FALSE]
  if (nrow(r) != 1L) stop("expected exactly one row of oc_cells_wide.csv, got ", nrow(r))
  num(r[[col]])
}
BU <- "Buderer_N500_Wald"; JO <- "Joint_N580_Wilson"; B1 <- "BAM_P1_N678"; B2 <- "BAM_P2_N672"
ex_success <- c(ocv(BU, 0.85, 0.90, 0.20, "p_success"), ocv(JO, 0.85, 0.90, 0.20, "p_success"),
                ocv(B1, 0.85, 0.90, 0.20, "p_success"), ocv(B2, 0.85, 0.90, 0.20, "p_success"),
                ocv(B1, 0.75, 0.90, 0.20, "p_success"), ocv(JO, 0.80, 0.90, 0.20, "p_success"),
                ocv(BU, 0.85, 0.90, 0.10, "p_success"), mpv(B1, "S02", "p_success"))
ex_cov_se <- c(ocv(BU, 0.85, 0.90, 0.20, "cov_se"), ocv(JO, 0.85, 0.90, 0.20, "cov_se"),
               ocv(B1, 0.85, 0.90, 0.20, "cov_se"), ocv(B2, 0.85, 0.90, 0.20, "cov_se"))
anchor("S3.9e-a", "exact success probabilities of the cross-check cells equal those of operating_characteristics.R / misspecified_priors.R at four decimals",
       xc$exact_success, round(ex_success, 4), tol = 1e-12)
anchor("S3.9e-b", "exact coverage of Se of the first four cells equal those of operating_characteristics.R at four decimals",
       xc$exact_cov_se[1:4], round(ex_cov_se, 4), tol = 1e-12)
mpd <- function(d, e, q) mpv(d, e, sprintf("det_%s_two_sided_analysis", q))
ex_det <- rbind(c(mpd(B1, "S00", "se"), mpd(B1, "S00", "sp"), mpd(B1, "S00", "any")),
                c(mpd(B1, "S02", "se"), mpd(B1, "S02", "sp"), mpd(B1, "S02", "any")),
                c(mpd(B2, "S04", "se"), mpd(B2, "S04", "sp"), mpd(B2, "S04", "any")))
anchor("S3.9f-a", "exact detection probabilities of the cross-check scenarios equal those of misspecified_priors.R at four decimals",
       c(xd$exact_se, xd$exact_sp, xd$exact_any), round(c(ex_det[, 1], ex_det[, 2], ex_det[, 3]), 4), tol = 1e-12)
anchor("S3.9e-z", "largest |z| over the success and coverage cells", max(abs(c(xc$z_success, xc$z_cov_se)), na.rm = TRUE), 1.66, dec = 2)
anchor("S3.9e-n", "success cells; coverage cells; detection scenarios", c(nrow(xc), sum(!is.na(xc$z_cov_se)), nrow(xd)), c(8L, 4L, 3L))
anchor("S3.9f-z", "all |z| <= 3 (success, coverage, detection)",
       max(abs(c(xc$z_success, xc$z_cov_se, xd$z_se, xd$z_sp, xd$z_any)), na.rm = TRUE) <= 3, TRUE)
fz <- function(v) ifelse(is.na(v), "", formatC(v, format = "f", digits = 2))
cell_lab <- c("Planning point (0.85, 0.90, 0.20)", "Planning point", "Planning point", "Planning point",
              "Se 0.75", "Se 0.80", "Prevalence 0.10", "Se 0.65 (S02)")
des_lab <- c("Buderer, Wald, N = 500", "Joint frequentist, Wilson, N = 580", "Bayesian P1, N = 678", "Bayesian P2, N = 672",
             "Bayesian P1", "Joint frequentist, Wilson", "Buderer, Wald", "Bayesian P1")
S39e <- data.frame(`Cell (Se, Sp, prevalence)` = cell_lab, Design = des_lab,
                   `Exact: success` = fp(xc$exact_success), `Monte Carlo: success` = fp(xc$success), z = fz(xc$z_success),
                   `Exact: coverage of Se` = fp(xc$exact_cov_se), `Monte Carlo: coverage of Se` = fp(xc$cov_se), `z (coverage)` = fz(xc$z_cov_se),
                   check.names = FALSE, stringsAsFactors = FALSE)
emit_main("R1_AF3_S3.9e_mc_crosscheck_success_coverage", xc, S39e)
sc_lab <- sub("^(P[12]) ", "\\1, ", xd$label)
rows <- list()
for (i in seq_len(nrow(xd))) {
  s1 <- sub(" \\(.*$", "", sc_lab[i])
  first <- TRUE
  add <- function(q, ex, mc, z) {
    rows[[length(rows) + 1L]] <<- data.frame(Scenario = if (first) sc_lab[i] else s1, Quantity = q,
      Exact = fp(ex), `Monte Carlo` = fp(mc), z = fz(z), check.names = FALSE, stringsAsFactors = FALSE)
    first <<- FALSE
  }
  if (xd$exact_se[i] > 0.001) add("Se", xd$exact_se[i], xd$det_se[i], xd$z_se[i])
  add("Sp", xd$exact_sp[i], xd$det_sp[i], xd$z_sp[i])
  # "any" is shown with a z only when it differs from the Sp row
  add("any", xd$exact_any[i], xd$det_any[i], if (xd$exact_any[i] > xd$exact_sp[i] + 1e-4) xd$z_any[i] else NA_real_)
}
S39f <- do.call(rbind, rows)
emit_main("R1_AF3_S3.9f_mc_crosscheck_detection", xd, S39f)

# --- Additional file 2: S2.3 negative controls, S2.4 empty arm
nc <- rd("negctl", "negative_control_halfwidth.csv")
nc$reading_label <- ifelse(startsWith(nc$reading, "(i) "), "(i) published full-width targets",
                    ifelse(startsWith(nc$reading, "(ii)"), "(ii) δ/2 used as the full-width target",
                    ifelse(startsWith(nc$reading, "(iii)"), "(iii) 2δ used as the full-width target", nc$reading)))
anchor("S2.3a", "joint assurance at the published N under reading (ii): P1; P2", c(nc$max_assurance_evaluated[nc$scenario == "P1" & startsWith(nc$reading, "(ii)")],
       nc$max_assurance_evaluated[nc$scenario == "P2" & startsWith(nc$reading, "(ii)")]), c(0.0032, 0.0381), dec = 4)
emit_simple("R1_AF2_S2.3_negative_controls", nc)
py05 <- file.path(script_dir, "python", "logs", "conf_PY_05_shared_assumptions.py.log")
pyl <- readLines(py05, warn = FALSE)
cat(sprintf("INPUT (package tree) %s sha256=%s\n", py05, unname(tools::sha256sum(py05))))
py_rows <- data.frame(case = rep(c("P1", "P2"), each = 3), variant = rep(c("V2 HPD (shortest) interval", "V5 n_d fixed at round(N*E[prev])", "V3 delta treated as HALF width"), 2),
                      N_star = c(666, 533, NA, 658, 533, NA), shift = c(-12, -145, NA, -14, -139, NA),
                      note = c("", "", "N* > 900; assurance at N = 678 is 0.003178", "", "", "N* > 900; assurance at N = 672 is 0.038115"),
                      source = "validation/python/logs/conf_PY_05_shared_assumptions.py.log, L2-L14", stringsAsFactors = FALSE)
okpy <- vapply(seq_len(nrow(py_rows)), function(i) {
  key <- sub(" .*$", "", py_rows$variant[i])
  ln <- pyl[startsWith(pyl, py_rows$case[i]) & grepl(paste0("\\b", key, "\\b"), pyl)]
  length(ln) == 1L && (is.na(py_rows$N_star[i]) || grepl(sprintf("\\b%d\\b", py_rows$N_star[i]), ln))
}, logical(1))
anchor("S2.3b", "constants of conf_PY_05 found in its log (HPD, fixed number of diseased, delta as half width)", sum(okpy), nrow(py_rows))
emit_simple("R1_AF2_S2.3_conf_PY_05_constants", py_rows)
ea <- rd("emptyarm", "edge_case_empty_arm.csv")
emit_simple("R1_AF2_S2.4_empty_arm", ea)

# =============================================================================
# Summary
# =============================================================================
cat("\n== INPUT FILES (path, SHA-256) ==\n")
for (i in seq_len(nrow(INPUT_LOG))) cat(sprintf("  --%s %s %s\n", INPUT_LOG$arg[i], INPUT_LOG$path[i], INPUT_LOG$sha256[i]))
cat(sprintf("\nelapsed: %.1f s\n", as.numeric(difftime(Sys.time(), T_START, units = "secs"))))
if (length(ANCH$fail) == 0L) {
  cat(sprintf("SUMMARY: %d anchors, all PASS\n", ANCH$n))
} else {
  cat(sprintf("SUMMARY: %d anchors, %d PASS, %d **CHECK**: %s\n", ANCH$n, ANCH$pass, length(ANCH$fail), paste(ANCH$fail, collapse = ", ")))
}
sink(); close(logcon)
if (length(ANCH$fail)) quit(save = "no", status = 1)
