# Validation scripts

Reproducible checks for `dtasamplesize`. They are not part of the installed
package (excluded from the build via `.Rbuildignore`); they document that
the package's numerical results are correct and reproducible.

Install the package first, then run from the repository root:

```sh
R -f validation/reproduce_manuscript.R    # every number reported in the article
R -f validation/make_manuscript_assets.R  # every figure and data-derived table
R -f validation/cross_validation.R        # core formulas vs independent references
R -f validation/integration_check.R       # end-to-end run + reproducibility
```

- **`reproduce_manuscript.R`** — regenerates every headline number quoted in
  the accompanying article (the assurance gap and its coverage, the
  net-benefit sample size and its threshold profile, the fixed-margin versus
  cohort standard deviations, the closed-form cohort variance against Monte
  Carlo, the apparent sensitivity and its bias, the five-step nested
  sequence behind Figure 4, and the harmonized comparison table of all five
  planning approaches under common assumptions), printing each computed
  value next to the value the article reports. It requires version 0.5.0 or
  later and stops with an explanatory message on earlier versions, whose
  results differ.
- **`make_manuscript_assets.R`** — the public generator behind every figure
  and data-derived table in the article: Figures 1–4 (written as a vector
  PDF and a 600 dpi TIFF to `submission-BMC-MRM/figures/`), the harmonized
  comparison table, Table 1 (derived directly from `NAMESPACE` and the
  installed help pages, not maintained by hand), and Table 3 (the feature
  matrix against `presize`, `pROC`, `MKpower`, `epiR` and `MKmisc`, read
  from the versioned `table3_feature_matrix_source.csv` in this directory,
  together with a provenance record of the comparator package versions
  reviewed). Table CSVs are written to `R-package/manuscript_assets/`.
- **`table3_feature_matrix_source.csv`** — the hand-curated data behind
  Table 3. Kept as a separate, versioned file (rather than embedded in the
  script) so the curation itself is visible and diffable on its own.
- **`cross_validation.R`** — checks the Wilson/Wald intervals against base R
  `prop.test` and `Hmisc::binconf`, `buderer_n` against published values,
  the Hanley–McNeil AUC variance against Monte Carlo, the net-benefit point
  estimate and standard error against the Vickers definition and Monte
  Carlo, the Beta–Binomial posterior against numerical integration, and the
  imperfect-reference variance inflation factor against the Rogan–Gladen
  formula. Alongside the printed report, it writes Table 2 ("Cross-validation
  of core formulas against independent references") to
  `R-package/manuscript_assets/table2_crossvalidation.csv`, with each row's
  result derived from the values the script itself computes -- including the
  observed value, the reference value, the tolerance used, and a pass flag,
  auditable against the printed report. `Hmisc` is a required reference
  package; `pROC` is listed as an optional one and is loaded only if
  available (no check currently calls it).
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
