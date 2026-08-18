# Validation scripts

Reproducible checks for `dtasamplesize`. They are not part of the installed
package (excluded from the build via `.Rbuildignore`); they document that
the package's numerical results are correct and reproducible.

Install the package first, then run from the repository root:

```sh
R -f validation/reproduce_manuscript.R  # every number reported in the article
R -f validation/cross_validation.R      # core formulas vs independent references
R -f validation/integration_check.R     # end-to-end run + reproducibility
```

- **`reproduce_manuscript.R`** — regenerates every headline number quoted in
  the accompanying article (the assurance gap and its coverage, the
  net-benefit sample size and its threshold profile, the fixed-margin versus
  cohort standard deviations, the closed-form cohort variance against Monte
  Carlo, the apparent sensitivity and its bias, and the five planning
  approaches at a common prevalence), printing each computed value next to
  the value the article reports. It requires version 0.3.0 or later and stops
  with an explanatory message on earlier versions, whose results differ.
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

## A note on data

The package ships no dataset, and the article analyses no empirical data. Every
quantity is computed from samples drawn inside the package's functions from the
stated parameters (sensitivity, specificity, prevalence, and the Monte Carlo
size `B`) under a fixed seed, so the parameters together with the seed determine
the data exactly. Re-running any of the scripts above regenerates the same
samples byte for byte.

If a materialised copy of the simulated samples is needed, set `DUMP_DATA <- TRUE`
at the top of `reproduce_manuscript.R`; it writes them as CSV to
`validation/simulated_data/`. Those files are a convenience only — the script,
not the CSV, is the authoritative record of how the data arise.
