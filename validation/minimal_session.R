# Minimal session of the R1 manuscript
# ---------------------------------------------------------------------
# The code between the two "MANUSCRIPT CODE" markers is the minimal
# session printed in the manuscript. Everything else here only sends the
# console output to a log file. Runs in well under one minute.
#
#   Rscript validation/minimal_session.R [output_dir]
#
# Output directory: first argument; default validation/output/ inside the
# package tree (created if needed). The log is minimal_session.log.

locate_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  if (basename(getwd()) == "validation") return(normalizePath(getwd()))
  if (dir.exists("validation")) return(normalizePath(file.path(getwd(), "validation")))
  stop("Run it via: Rscript validation/minimal_session.R [output_dir]",
       call. = FALSE)
}
cli_args <- commandArgs(trailingOnly = TRUE)
out_dir  <- if (length(cli_args) >= 1 && nzchar(cli_args[1])) cli_args[1] else
  file.path(locate_script_dir(), "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_dir  <- normalizePath(out_dir, winslash = "/")
log_file <- file.path(out_dir, "minimal_session.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)
options(width = 72)   # keeps printed tables within the PDF line width
t0 <- Sys.time()
cat("R version:        ", R.version.string, "\n")
cat("Output directory: ", out_dir, "\n\n")

# ===================== MANUSCRIPT CODE (start) =======================
library(dtasamplesize)
packageVersion("dtasamplesize")

# 1. The motivating result, exact (no simulation): how often does a
#    study sized with Buderer's formula reach the target CI width?
mc_validate_buderer(Se = 0.85, d = 0.07, method = "exact")

# 2. The worked example: Se ~ Beta(17, 3), Sp ~ Beta(2, 2),
#    prevalence ~ Beta(4, 16), total CI widths 0.14 (Se) and 0.10 (Sp),
#    95% intervals, joint assurance target 0.80, exact calculation.
bam <- bam_sample_size(
  prior_se = c(17, 3), prior_sp = c(2, 2), prior_prev = c(4, 16),
  delta_se = 0.14, delta_sp = 0.10, target_assurance = 0.80,
  method = "exact"
)
bam

# 3. A batch of prior scenarios (sensitivity to the prior):
#    vague vs informative prior for Sp, and a weaker prior for Se.
scenarios <- list(
  sp_vague       = list(prior_sp = c(2, 2)),
  sp_informative = list(prior_sp = c(18, 2)),
  se_weaker      = list(prior_se = c(8.5, 1.5))
)
bam_sensitivity(scenarios, prior_se = c(17, 3), prior_sp = c(2, 2),
                prior_prev = c(4, 16), delta_se = 0.14, delta_sp = 0.10)

sessionInfo()
# ====================== MANUSCRIPT CODE (end) ========================

cat(sprintf("\nElapsed: %.1f seconds\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
cat("Log written to", log_file, "\n")
