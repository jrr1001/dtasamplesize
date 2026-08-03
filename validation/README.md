# Validation scripts

Reproducible checks for `dtasamplesize`. They are not part of the installed
package (excluded from the build via `.Rbuildignore`); they document that
the package's numerical results are correct and reproducible.

Install the package first, then run from the repository root:

```sh
R -f validation/cross_validation.R    # core formulas vs independent references
R -f validation/integration_check.R   # end-to-end run + reproducibility
```

- **`cross_validation.R`** — checks the Wilson/Wald intervals against base R
  `prop.test` and `Hmisc::binconf`, `buderer_n` against published values,
  the Hanley–McNeil AUC variance against Monte Carlo, the net-benefit point
  estimate and standard error against the Vickers definition and Monte
  Carlo, the Beta–Binomial posterior against numerical integration, and the
  imperfect-reference variance inflation factor against the Rogan–Gladen
  formula. Optional reference packages: `Hmisc`, `pROC`.
- **`integration_check.R`** — runs every estimator, confirms the
  "more uncertainty → larger N" ordering, and verifies reproducibility under
  a fixed seed.

A narrative summary of the results is in `../validation_report.md`.
