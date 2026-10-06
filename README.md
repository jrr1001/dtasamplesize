# dtasamplesize

Simulation-based sample size estimation for diagnostic test accuracy
(DTA) studies. Version 0.6.8.

Every simulation function takes a `seed`, is exactly reproducible, and
restores the caller's random-number state on exit, so calling it never
perturbs your own random stream.

## Installation

```r
# install.packages("devtools")
devtools::install_github("jrr1001/dtasamplesize")
```

## Quick start

```r
library(dtasamplesize)

# How often does the classical Buderer sample size deliver its target
# precision? (Only about 57% of the time.)
mc_validate_buderer(Se = 0.85, d = 0.07, B = 5000)

# Sample size to demonstrate positive net benefit. The default design is
# "cohort" (disease status random, as in a real prospective DTA study);
# use design = "fixed" only when the arm sizes are fixed by design.
ss_net_benefit(pt_range = c(0.10, 0.30), N_range = seq(50, 400, by = 10))
```

## Overview

Classical sample size formulas for DTA studies (Buderer 1996,
Hanley & McNeil 1982) assume ideal conditions that rarely hold in
practice. `dtasamplesize` provides simulation-based alternatives that
account for parameter uncertainty, imperfect reference standards,
time-to-event outcomes, and multiple accuracy metrics simultaneously.

## Functions

| Function | Purpose |
|----------|---------|
| `mc_validate_buderer()` | Validate classical formula via Monte Carlo |
| `bam_sample_size()` | Bayesian Assurance Method |
| `joint_sample_size()` | Joint Se + Sp + AUC optimization |
| `ss_imperfect_ref()` | Imperfect reference standard correction (**experimental**, see below) |
| `ss_time_dependent_roc()` | Time-dependent ROC sample size (**experimental**, see below) |
| `ss_adaptive_prevalence()` | Adaptive prevalence re-estimation (**experimental**, see below) |
| `ss_net_benefit()` | Net Benefit / Decision Curve Analysis (`design = "cohort"` or `"fixed"`) (**experimental**, see below) |
| `ss_unified()` | Unified framework combining all methods (**experimental**, see below) |
| `sensitivity_analysis()` | Recompute sample size across a grid of alternative prior scenarios |
| `plot_assurance_curve()` | Plot achieved assurance vs. sample size from a unified search (requires ggplot2) |
| `plot_inflation_heatmap()` | Plot the imperfect-reference variance inflation heatmap (requires ggplot2) |
| `plot_method_comparison()` | Plot required sample size by method from a unified search (requires ggplot2) |

## Validated scope and experimental modules

The validated scope of the package -- the methods checked against the
dtasamplesize software article's Monte Carlo/exact-mode validation
evidence (R CMD check, `tests/testthat`, and the confirmatory validation
pipeline under `validation/`) -- comprises `bam_sample_size(method =
"exact")` and `joint_sample_size()`, together with the helper functions
`mc_validate_buderer()`, `buderer_n()`, `wilson_ci()`, and `wald_ci()`.

The following five modules are **experimental**: they are exported and
documented, but are outside that validated scope and have not been
independently validated against an external reference implementation.
Do not rely on them for design decisions until validated separately.

- `ss_unified()`
- `ss_net_benefit()`
- `ss_imperfect_ref()`
- `ss_adaptive_prevalence()`
- `ss_time_dependent_roc()`

## References

- Buderer NMF (1996). Acad Emerg Med. doi:10.1111/j.1553-2712.1996.tb03538.x
- Wilson KJ et al. (2022). Stat Med. doi:10.1002/sim.9393
- Stark M, Zapf A (2020). Stat Methods Med Res. doi:10.1177/0962280220913588
- Staquet M et al. (1981). J Chronic Dis. doi:10.1016/0021-9681(81)90059-X
