# dtasamplesize

Sample size methods for diagnostic test accuracy (DTA) studies: exact
Bayesian assurance and related methods. Version 0.7.0.

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
# precision? Exact probability (binomial enumeration, no simulation) that
# the 95% CI is no wider than 2 * d: 0.5683 for Se = 0.85, d = 0.07.
mc_validate_buderer(Se = 0.85, d = 0.07, method = "exact")

# Bayesian assurance, exact (no simulation): smallest total N whose joint
# assurance for Se and Sp reaches 0.80; read N_total and joint_assurance.
res <- bam_sample_size(prior_se = c(17, 3), delta_se = 0.14, method = "exact")
res$N_total
res$joint_assurance

# The same calculation for several prior scenarios (N_total, A(N), A(N - 1)).
bam_sensitivity(list(vague = list(prior_sp = c(2, 2)),
                     informative = list(prior_sp = c(18, 2))))

# Sample size to demonstrate positive net benefit. The default design is
# "cohort" (disease status random, as in a real prospective DTA study);
# use design = "fixed" only when the arm sizes are fixed by design.
ss_net_benefit(pt_range = c(0.10, 0.30), N_range = seq(50, 400, by = 10))
```

## Overview

Classical sample size formulas for DTA studies (Buderer 1996,
Hanley & McNeil 1982) assume ideal conditions that rarely hold in
practice. `dtasamplesize` provides exact and simulation-based alternatives that
account for parameter uncertainty, imperfect reference standards,
time-to-event outcomes, and multiple accuracy metrics simultaneously.

## Functions

| Function | Purpose |
|----------|---------|
| `mc_validate_buderer()` | Probability that the classical Buderer sample size reaches its target width, by Monte Carlo (default) or exact enumeration (`method = "exact"`) |
| `bam_sample_size()` | Bayesian Assurance Method (`method = "exact"`: Beta-Binomial enumeration; separate design and analysis priors via `analysis_prior_se`/`analysis_prior_sp`; `width_type`; `assurance_for`; `priors`) |
| `bam_prior_from_pilot()` | Power prior with fixed weight `a0` from pilot counts, for `bam_sample_size(priors = )` |
| `bam_prior_data_conflict()` | Exact prior-data conflict tail probabilities for observed 2x2 counts |
| `bam_sensitivity()` | Exact sample size (`N_total`, `A_N`, `A_N_minus_1`) over a set of scenarios |
| `bam_prior_contour()`, `bam_sensitivity_local()` | Local prior sensitivity analysis (Wilson et al. 2022, section 4.2) |
| `bam_assurance_curve()` | Exact assurance A(N) over a range of total sample sizes |
| `joint_sample_size()` | Joint Se + Sp search by Monte Carlo, with a deterministic AUC gate |
| `ss_imperfect_ref()` | Imperfect reference standard correction (**experimental**, see below) |
| `ss_time_dependent_roc()` | Time-dependent ROC sample size (**experimental**, see below) |
| `ss_adaptive_prevalence()` | Adaptive prevalence re-estimation (**experimental**, see below) |
| `ss_net_benefit()` | Net Benefit / Decision Curve Analysis (`design = "cohort"` or `"fixed"`) (**experimental**, see below) |
| `ss_unified()` | Unified framework combining all methods (**experimental**, see below) |
| `sensitivity_analysis()` | Recompute `ss_unified()` across a grid of alternative prior scenarios (**experimental**; for the validated method use `bam_sensitivity()`) |
| `plot_assurance_curve()` | Plot assurance vs. sample size from `bam_assurance_curve()` (exact) or from a unified search (**experimental**) (requires ggplot2) |
| `plot_inflation_heatmap()` | Plot the imperfect-reference variance inflation heatmap (requires ggplot2) |
| `plot_method_comparison()` | Plot required sample size by method from a unified search (requires ggplot2) |

## Validated scope and experimental modules

The validated scope of the package, as declared in DESCRIPTION, comprises
`bam_sample_size(method = "exact")` and the other `bam_*` functions
(`bam_assurance_curve()`, `bam_sensitivity()`, `bam_sensitivity_local()`,
`bam_prior_contour()`, `bam_prior_from_pilot()` and
`bam_prior_data_conflict()`), `joint_sample_size()` and
`mc_validate_buderer(method = "exact")`. The helpers `buderer_n()`,
`wilson_ci()` and `wald_ci()` lie outside that scope and are checked
against external software and published values;
`mc_validate_buderer(method = "monte_carlo")`, the default, is auxiliary.

The Python implementation under `validation/python/` is a separate
re-implementation written by an LLM-based coding assistant under the authors'
supervision, not external software; agreement with it can reveal
implementation errors but cannot exclude a misreading of the specification
shared by both.

The following five modules are **experimental**: they are exported and
documented, but are outside that validated scope and have not been
validated against an external reference implementation. Do not rely on
them for design decisions until validated separately.

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
