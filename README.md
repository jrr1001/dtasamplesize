# dtasamplesize

Simulation-based sample size estimation for diagnostic test accuracy
(DTA) studies.

## Installation

```r
# install.packages("devtools")
devtools::install_github("jrr1001/dtasamplesize")
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
| `ss_imperfect_ref()` | Imperfect reference standard correction |
| `ss_time_dependent_roc()` | Time-dependent ROC sample size |
| `ss_adaptive_prevalence()` | Adaptive prevalence re-estimation |
| `ss_net_benefit()` | Net Benefit / Decision Curve Analysis |
| `ss_unified()` | Unified framework combining all methods |

## References

- Buderer NMF (1996). Acad Emerg Med. doi:10.1111/j.1553-2712.1996.tb03538.x
- Wilson KJ et al. (2022). Stat Med. doi:10.1002/sim.9393
- Stark M, Zapf A (2020). Stat Methods Med Res. doi:10.1177/0962280220913588
- Staquet M et al. (1981). J Chronic Dis. doi:10.1016/0021-9681(81)90059-X
