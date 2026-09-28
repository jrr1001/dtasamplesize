#!/usr/bin/env Rscript
## sweep_wald_wilson.R
##
## Purpose (H-02 of the 2026-09-27 blind audit): manuscript.md:75 asserts
## "A systematic sweep of 28 configurations (Se from 0.60 to 0.98, d from
## 0.03 to 0.10, each evaluated at its own Buderer sample size with 200000
## replicates)" for the Wald vs. Wilson comparison of P(CI width <= target).
## No script in this package (or in any earlier tagged version) generates
## that sweep -- it was typed by hand. This script is the generator.
##
## Grid declaration and justification
## -----------------------------------
## The manuscript names only the endpoints of each axis (Se: 0.60-0.98,
## d: 0.03-0.10) and the cell count (28 = 7 x 4), not the intermediate
## values. A grep of every Rd, R, NEWS.md and manuscript file in the
## package repository (see H-02 audit log) turns up no prior declaration
## of the intermediate grid points. We therefore choose a grid that:
##   (a) has exactly 7 Se values and 4 d values (7 x 4 = 28, as stated);
##   (b) includes both named endpoints on each axis (Se = 0.60 and 0.98;
##       d = 0.03 and 0.10);
##   (c) includes the single already-published control cell,
##       Se = 0.85, d = 0.07, so that cell reproduces P = 0.5730 (Wald)
##       inside the sweep as an internal check;
##   (d) spaces the remaining points at round, evenly-readable increments
##       typical of a sensitivity/precision sweep in this literature.
## Se: 0.60, 0.70, 0.75, 0.80, 0.85, 0.90, 0.98   (7 values)
## d :  0.03, 0.05, 0.07, 0.10                     (4 values)
## This is A reasonable grid consistent with every stated fact, not THE
## grid the original (unrecovered) script used -- there is no way to
## recover that with certainty, since the manuscript did not archive it.
## Per audit instructions, if the extremes below are grid-independent
## (driven only by Se and d at the boundary, not by the interior points),
## the published ranges should reproduce regardless of exactly which
## interior points were chosen.
##
## Method
## ------
## For each (Se, d) cell: n <- buderer_n(Se, d) (package function), then
## mc_validate_buderer(Se, d, n_diseased = n, B = 200000, ci_method = "both",
## seed = <fixed>) using the package's own Wald/Wilson width criterion:
## target full width = 2*d, P_width_target = mean(ci_width <= 2*d), where
## ci_width is computed exactly as in R/mc_validate_buderer.R (Wald normal
## approximation and Wilson score interval on Se_hat = x/n, x ~ Binomial(n,Se)).
## This is exactly the criterion used to obtain the published single-case
## value 0.5730 at Se=0.85, d=0.07 (validation/reproduce_manuscript.R,
## section 1), except that script used B=4000; this sweep uses B=200000 as
## stated in the manuscript sentence being audited. The same fixed seed
## (2026, the seed used throughout the manuscript; see reproduce_manuscript.R)
## is used for every cell, applied through mc_validate_buderer()'s own
## seed argument (which reseeds the Mersenne-Twister with explicit RNG kind
## before drawing, so results are reproducible independent of call order).
##
## A binomial Monte Carlo standard error is reported per cell/method:
## SE = sqrt(P*(1-P)/B).

suppressPackageStartupMessages({
  if (requireNamespace("pkgload", quietly = TRUE) &&
      file.exists(file.path("..", "DESCRIPTION"))) {
    pkgload::load_all("..", quiet = TRUE)
  } else if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(quiet = TRUE)
  } else {
    library(dtasamplesize)
  }
})

SEED <- 2026L
B    <- 200000L
ALPHA <- 0.05

Se_grid <- c(0.60, 0.70, 0.75, 0.80, 0.85, 0.90, 0.98)
d_grid  <- c(0.03, 0.05, 0.07, 0.10)

stopifnot(length(Se_grid) == 7, length(d_grid) == 4,
          length(Se_grid) * length(d_grid) == 28)
stopifnot(0.60 %in% Se_grid, 0.98 %in% Se_grid)
stopifnot(0.03 %in% d_grid, 0.10 %in% d_grid)
stopifnot(0.85 %in% Se_grid, 0.07 %in% d_grid)  # control cell present

cat("================================================================\n")
cat("sweep_wald_wilson.R -- H-02 generator for the 28-configuration sweep\n")
cat("================================================================\n\n")
cat("Grid (declared above the run, BEFORE inspecting published numbers):\n")
cat("  Se =", paste(Se_grid, collapse = ", "), "\n")
cat("  d  =", paste(d_grid, collapse = ", "), "\n")
cat("  cells =", length(Se_grid) * length(d_grid), "\n")
cat("  B =", B, ", seed =", SEED, ", alpha =", ALPHA, "\n\n")

grid <- expand.grid(Se = Se_grid, d = d_grid, KEEP.OUT.ATTRS = FALSE)
grid <- grid[order(grid$Se, grid$d), ]

rows <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  Se <- grid$Se[i]; d <- grid$d[i]
  n  <- buderer_n(Se, d, ALPHA)
  mv <- mc_validate_buderer(Se = Se, d = d, n_diseased = n, B = B,
                             ci_method = "both", alpha = ALPHA, seed = SEED)
  res <- mv$results
  p_wald   <- res$P_width_target[res$ci_method == "wald"]
  p_wilson <- res$P_width_target[res$ci_method == "wilson"]
  se_wald   <- sqrt(p_wald   * (1 - p_wald)   / B)
  se_wilson <- sqrt(p_wilson * (1 - p_wilson) / B)
  rows[[i]] <- data.frame(Se = Se, d = d, n = n,
                           P_wald = p_wald, MCse_wald = se_wald,
                           P_wilson = p_wilson, MCse_wilson = se_wilson)
}
tab <- do.call(rbind, rows)
rownames(tab) <- NULL

cat("---- Full 28-cell table ----\n")
print(tab, digits = 6)

cat("\n---- Observed ranges (recorded BEFORE looking at the published text) ----\n")
cat(sprintf("Wald:   min = %.4f  max = %.4f  (across %d cells)\n",
            min(tab$P_wald), max(tab$P_wald), nrow(tab)))
cat(sprintf("Wilson: min = %.4f  max = %.4f  (across %d cells)\n",
            min(tab$P_wilson), max(tab$P_wilson), nrow(tab)))
q <- quantile(tab$P_wald, probs = c(0.25, 0.5, 0.75))
cat(sprintf("Wald IQR: [%.4f, %.4f], median %.4f\n", q[1], q[3], q[2]))
cat("Max per-cell MC SE (Wald):   ", sprintf("%.5f", max(tab$MCse_wald)), "\n")
cat("Max per-cell MC SE (Wilson): ", sprintf("%.5f", max(tab$MCse_wilson)), "\n")

cat("\n---- Extreme cells ----\n")
cat("Wald min at: "); print(tab[which.min(tab$P_wald), c("Se","d","n","P_wald")])
cat("Wald max at: "); print(tab[which.max(tab$P_wald), c("Se","d","n","P_wald")])
cat("Wilson min at: "); print(tab[which.min(tab$P_wilson), c("Se","d","n","P_wilson")])
cat("Wilson max at: "); print(tab[which.max(tab$P_wilson), c("Se","d","n","P_wilson")])

cat("\n---- Se = 0.98 row (Wilson), all d ----\n")
print(tab[tab$Se == 0.98, c("Se","d","n","P_wilson")])

cat("\n---- Control cell: Se = 0.85, d = 0.07 (published: 0.5730) ----\n")
ctrl <- tab[tab$Se == 0.85 & tab$d == 0.07, ]
print(ctrl)
## The published 0.5730 was obtained with B = 4000 (reproduce_manuscript.R),
## not B = 200000 as used here for the sweep -- the two Monte Carlo runs
## consume the RNG stream differently even under the same seed, so the two
## estimates are independent draws of the SAME underlying probability, not
## expected to coincide exactly. The relevant tolerance is therefore the
## COMBINED MC standard error of the two runs (published B=4000 estimate
## and this B=200000 estimate), not the B=200000 SE alone.
se_b4000  <- sqrt(0.5730 * (1 - 0.5730) / 4000)
se_combined <- sqrt(ctrl$MCse_wald^2 + se_b4000^2)
ctrl_tol <- 3 * se_combined
cat(sprintf("MC SE at B=4000 (published run):  %.5f\n", se_b4000))
cat(sprintf("MC SE at B=200000 (this run):      %.5f\n", ctrl$MCse_wald))
cat(sprintf("Combined MC SE:                    %.5f  (tol = 3x = %.5f)\n",
            se_combined, ctrl_tol))
cat(sprintf("Control check: |P_wald - 0.5730| = %.5f vs tol = %.5f -> %s\n",
            abs(ctrl$P_wald - 0.5730), ctrl_tol,
            ifelse(abs(ctrl$P_wald - 0.5730) <= ctrl_tol, "PASS", "**CHECK**")))
stopifnot(abs(ctrl$P_wald - 0.5730) <= ctrl_tol)

cat("\n---- Comparison with the numbers currently in manuscript.md:75 ----\n")
cat("Published Wald range:   0.376 - 0.851 (concentrated 0.47-0.61)\n")
cat("Published Wilson range: 0.000 (Se=0.98, several d) to 1.000 (Se=0.60, d=0.10)\n")
cat(sprintf("This run   Wald range:   %.4f - %.4f\n", min(tab$P_wald), max(tab$P_wald)))
cat(sprintf("This run   Wilson range: %.4f - %.4f\n", min(tab$P_wilson), max(tab$P_wilson)))

cat("\n---- sessionInfo() ----\n")
print(sessionInfo())

cat("\nDone.\n")
