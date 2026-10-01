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
R -f validation/rng_invariance.R          # RNG state invariance of the exact mode
R -f validation/exact_invariance.R        # G04: exact-mode N/A invariance across B/seed/N_range
```

- **`reproduce_manuscript.R`** — regenerates every headline number quoted in
  the accompanying article, now reduced to a single joint estimand: the
  assurance gap of the classical Buderer sample size and its coverage
  (Figure 1), the exact joint Beta-Binomial assurance for Se and Sp and its
  crossing at N = 678 (Figure 2), and Table 4's comparison of the three
  surviving methods (Buderer classical, BAM exact mode with harmonized
  priors, and `joint_sample_size`) under common assumptions — printing each
  computed value next to the value the article reports. `ss_unified()`, the
  AUC gate inside it, `ss_net_benefit()`, `ss_imperfect_ref()`,
  `ss_adaptive_prevalence()` and `ss_time_dependent_roc()` remain in the
  package but are no longer described in the article, so this script no
  longer checks their numbers. It requires version 0.6.6 exactly and stops
  with an explanatory message on any other installed version, whose results
  are not guaranteed to reproduce the article's bit-for-bit.
- **`make_manuscript_assets.R`** — the public generator behind every figure
  and data-derived table in the article's current, reduced scope: Figure 1
  (unchanged) and Figure 2 (the joint-assurance-vs-N curve), written as a
  vector PDF and a 600 dpi TIFF to `submission-BMC-MRM/figures/`; Table 1
  (derived directly from `NAMESPACE` and the installed help pages, not
  maintained by hand); Table 3 (the feature matrix against `presize`,
  `pROC`, `MKpower`, `epiR`, `MKmisc` and `SampleSizeDiagnostics`, read from
  the versioned `table3_feature_matrix_source.csv` in this directory,
  together with a provenance record of the comparator package versions
  reviewed); and Table 4 (the three surviving methods compared under common
  assumptions). Table CSVs are written to `R-package/manuscript_assets/`.
  Table 2 (cross-validation of core formulas) is written separately, by
  `cross_validation.R`.
- **`table3_feature_matrix_source.csv`** — the hand-curated data behind
  Table 3. Kept as a separate, versioned file (rather than embedded in the
  script) so the curation itself is visible and diffable on its own.
- **`cross_validation.R`** — checks the Wilson/Wald intervals against base R
  `prop.test` and `Hmisc::binconf`, `buderer_n` against published values,
  and the Beta posterior interval width against numerical integration (a
  separate check verifies the Beta-Binomial predictive mass itself, by
  finite sum, simulation and the `lbeta`-based closed form; it is not one
  of the four published rows). These four checks are the ones Table 2
  ("Cross-validation of core formulas against
  independent references") in the article publishes; alongside the printed
  report, they are written to
  `R-package/manuscript_assets/table2_crossvalidation.csv`, one row per
  check, with each row's result derived from the values the script itself
  computes -- including the observed value, the reference value, the
  tolerance used, and a pass flag, auditable against the printed report.
  The `package_version` column records the version of the package actually
  installed when the script runs (`packageVersion("dtasamplesize")`), not a
  fixed string, so it changes whenever the CSV is regenerated against a
  different release.
  The script also runs the same style of independent-reference check for
  three functions that remain in the package but are no longer described in
  the article -- the Hanley–McNeil AUC variance against Monte Carlo, the
  net-benefit point estimate and standard error against the Vickers
  definition and Monte Carlo, and the imperfect-reference variance
  inflation factor and apparent-sensitivity closed form against the
  Rogan–Gladen formula. Those checks still print PASS/`**CHECK**` to the
  console as an internal sanity check on package code that ships, but they
  are deliberately not written to the CSV, since they no longer correspond
  to a row the article publishes. In the Hanley–McNeil check the 0.20
  tolerance is applied only to the balanced configurations: that formula
  assumes a negative-exponential score model while the Monte Carlo
  reference is binormal, so the imbalanced configuration is reported as an
  expected deviation rather than as a failure. That is the criterion the
  script's own summary variable has always used. `Hmisc` is a required
  reference package; `pROC` is listed as an optional one and is loaded only
  if available (no check currently calls it).
- **`integration_check.R`** — a package-wide smoke test, not scoped to the
  article: it runs every exported estimator, including the ones the article
  no longer describes (`ss_net_benefit`, `ss_imperfect_ref`,
  `ss_adaptive_prevalence`, `ss_time_dependent_roc`, `ss_unified`), prints
  the sample size each one returns, and verifies reproducibility under a
  fixed seed. It also asserts the "more uncertainty → larger N" ordering,
  in a clearly delimited "Ordering check" section near the end of the
  script: it computes `bam_sample_size(prior_se = p, method = "exact")`
  (all other arguments default) for the Se prior triple Beta(34, 6),
  Beta(17, 3), Beta(8.5, 1.5), prints the three `N_total` values the script
  actually computes, and `stop()`s the script if they are not strictly
  increasing. (Earlier versions of this file and of the script's own header
  asserted the ordering only in prose, with no code behind the claim; that
  gap is what this section closes.)

- **`rng_invariance.R`** — confirms that `bam_sample_size(method = "exact")`
  and `joint_sample_size()` save and restore `.Random.seed`/`RNGkind()`
  exactly, and that the exact-mode closed-form path (no `B`/`seed`
  dependence) returns identical results across RNG kinds.
- **`exact_invariance.R`** (G04) — a dedicated check that the exact mode of
  `bam_sample_size()` is invariant to `B`, `seed` and `N_range` (since the
  exact Beta-Binomial calculation does not use Monte Carlo draws at all):
  for the published N = 678 case and the Table 4 harmonized N = 672 case,
  it sweeps `B` in `{5, 10, 50, 5000}` crossed with `seed` in
  `{1, 2, 3, 4, 2026}` and `N_range` in `{NULL, an explicit range
  containing the crossing}`, and requires identical `N_total` and
  `joint_assurance` and identical `.Random.seed`/`RNGkind()` before and
  after every call. It can also be pointed at the 0.6.5 library (see the
  script's own `--version-check=off` flag) to confirm that the same check
  FAILS there, demonstrating that it actually detects the defect this
  correction fixes rather than passing vacuously.

## A note on what these scripts prove

`reproduce_manuscript.R` and `integration_check.R` are reproducibility
checks: they confirm that the installed package, run today, prints the same
numbers the article reports (or that it runs end-to-end, self-consistently,
under a fixed seed). A PASS there is not a validity check -- a defect that
biases a result reproduces exactly as cleanly as a result that is correct,
and this history is not hypothetical: the numbers `reproduce_manuscript.R`
now checks were themselves revised twice across earlier package versions
after those versions passed their own reproducibility checks against
defective numbers (see `NEWS.md`). `cross_validation.R` is closer to a
genuine validity check for the specific formulas it covers, since it
compares the package's output against an independent implementation (base
R, `Hmisc`, or direct numerical integration) rather than against the
package's own past output -- but it still only covers the four quantities
in Table 2. Validity beyond what these scripts test is established
separately, release by release, by the package's own documentation and by
independent review of the methods, never by a script agreeing with itself.

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
