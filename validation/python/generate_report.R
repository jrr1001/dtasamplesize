# Regenerate validation/python/report.html and the "Results" table in
# validation/python/README.md from the current run's logs.
#
# This replaces the undocumented `tools/generar_report_html.R` referenced by
# earlier revisions of README.md (under a `$WORK/tools/` path specific to the
# machine that ran the original 2026-08-27/2026-09-16 Python confirmatory
# validation; that script was never committed to this repository). This
# script is the public, versioned generator from here on: it reads the log
# files already committed under validation/logs/ and validation/python/logs/
# and writes report.html deterministically. It does not re-run any
# validation script itself -- run run_all.sh/run_all.ps1 and the top-level
# validation/*.R scripts first.
#
# Usage (from the package root):
#   Rscript validation/python/generate_report.R

repo_root <- normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."), "..", ".."),
                           mustWork = FALSE)
`%||%` <- function(a, b) if (is.null(a)) b else a

# Resolve paths relative to the repository root regardless of cwd.
find_root <- function() {
  d <- normalizePath(".", mustWork = TRUE)
  while (!file.exists(file.path(d, "DESCRIPTION")) && dirname(d) != d) {
    d <- dirname(d)
  }
  d
}
root <- find_root()
val <- file.path(root, "validation")
pyval <- file.path(val, "python")

read_log <- function(name, dir = file.path(val, "logs")) {
  f <- file.path(dir, name)
  if (!file.exists(f)) return(NA_character_)
  paste(readLines(f, warn = FALSE), collapse = "\n")
}

pkg_ver <- as.character(utils::packageVersion("dtasamplesize"))

check_log <- read_log("R_CMD_check_00check.log")
check_status <- {
  m <- regmatches <- regmatches(check_log, regexpr("Status:.*", check_log))
  if (length(m) && nzchar(m)) trimws(m) else "NOT FOUND"
}

repro <- read_log("reproduce_manuscript.log")
get_num <- function(text, pattern) {
  m <- regmatches(text, regexpr(pattern, text, perl = TRUE))
  if (length(m) && nzchar(m)) trimws(m) else "NOT FOUND"
}
n678_line <- get_num(repro, "N_total = 678[^\\n]*")
n672_line <- get_num(repro, "N_total = 672[^\\n]*")

exact_inv <- read_log("exact_invariance.log")
exact_inv_status <- if (grepl("PASS", exact_inv, fixed = TRUE) &&
                         !grepl("FAIL", exact_inv, fixed = TRUE)) "PASS" else "CHECK LOG"

exact_inv_hist <- read_log("exact_invariance_v065.log", dir = file.path(val, "historical_not_current"))
exact_inv_hist_status <- if (!is.na(exact_inv_hist) && grepl("FAIL", exact_inv_hist, fixed = TRUE)) {
  "FAIL (expected, historical 0.6.5 demonstration)"
} else "see historical_not_current/"

ts <- read_log("testthat.Rout")
testthat_summary <- get_num(ts, "\\[ FAIL [^\\n]*PASS [0-9]+ \\]")

out_html <- sprintf('<title>BAM Assurance Confirmation (%s)</title>
<style>
body{font-family:ui-monospace,Consolas,monospace;max-width:72ch;margin:2em auto;line-height:1.5}
table{border-collapse:collapse;width:100%%}
td,th{border:1px solid #888;padding:4px 8px;text-align:left}
.fail{color:#963528;font-weight:bold}
.pass{color:#2C6A43;font-weight:bold}
</style>
<h1>BAM joint-assurance confirmation -- dtasamplesize %s</h1>
<p>Regenerated automatically by validation/python/generate_report.R from the
logs under validation/logs/ produced by validation/python/run_all.sh (or
run_all.ps1) and the top-level validation/*.R scripts. Generated at
%s. This file replaces hand edits; to update it, re-run the
logged validation scripts and this generator, never edit the HTML
directly.</p>
<table>
<tr><th>Item</th><th>Value</th><th>Source log</th></tr>
<tr><td>Package version (sessionInfo)</td><td>%s</td><td>sessionInfo.txt</td></tr>
<tr><td>R CMD check status</td><td>%s</td><td>R_CMD_check_00check.log</td></tr>
<tr><td>testthat summary</td><td>%s</td><td>testthat.Rout</td></tr>
<tr><td>G04 exact_invariance (current, expected PASS)</td><td class="pass">%s</td><td>exact_invariance.log</td></tr>
<tr><td>G04 exact_invariance (0.6.5, historical, expected FAIL)</td><td class="fail">%s</td><td>historical_not_current/exact_invariance_v065.log</td></tr>
<tr><td>N=678 reproduced</td><td>%s</td><td>reproduce_manuscript.log</td></tr>
<tr><td>N=672 reproduced</td><td>%s</td><td>reproduce_manuscript.log</td></tr>
</table>
<p>See validation/python/README.md for the full historical provenance
(2026-08-27/2026-09-16/2026-09-28 Python confirmatory run, superseded
material under historical_not_current/) and validation/README.md for what
each script proves and does not prove.</p>
',
  pkg_ver, pkg_ver, format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
  pkg_ver, check_status, testthat_summary, exact_inv_status,
  exact_inv_hist_status, n678_line, n672_line)

writeLines(out_html, file.path(pyval, "report.html"))
cat("wrote", file.path(pyval, "report.html"), "\n")
cat("package version:", pkg_ver, "\n")
cat("check status:", check_status, "\n")
cat("testthat:", testthat_summary, "\n")
