# Validation scripts

Reproducible checks for `dtasamplesize`. They are not part of the installed
package (excluded from the build via `.Rbuildignore`); they document that
the package's numerical results are reproducible and compare them with separate
computations by the same team and with external software (see "A note on what
these scripts prove").

The logs and output files in this directory were regenerated with
dtasamplesize 0.7.0 (built from the release tree and installed in a private
library, not the system library) on 2026-11-23. The current manifest of
SHA-256 hashes is `HASHES_2026-11-23_v0.7.0.txt`; the manifests of earlier
releases (`HASHES_2026-10-06_v0.6.9.txt`, `HASHES_2026-10-01_v2.txt` and the
files in `historical_not_current/`) are kept for provenance and are not
current evidence. The logs that print the package version show 0.7.0;
`rng_invariance.R`, `integration_check.R`,
`mc_validate_buderer_exact_vs_reference.R` and the registered scripts
`python/conf_R_03_properties.R` and `python/conf_R_04_mc_vs_exact.R` do not
print it and were run against the same installed library, which
`logs/sessionInfo.txt` records. The manifest hashes the files as stored in
the repository; `.gitattributes` in this directory turns off line-ending
conversion, so a clone or release archive on any platform reproduces them
(`cd validation && sha256sum -c HASHES_2026-11-23_v0.7.0.txt`).
The logs and outputs were produced on Windows, where R writes CRLF line endings,
and are stored with LF; the SHA-256 that a log prints for its input files refers
to those CRLF working copies, so compare after normalising line endings. The logs
record the commit of the working tree at the time of the run, before the release
commit was made.

Install the package first, then run from the repository root. Scripts that
write files take the output directory as their first argument (default
`validation/output/`, inside the package tree, ignored by git and by the
package build). No script writes outside the directory it is given.

```sh
# Scripts that were in 0.6.9 (console log in logs/<script>.log)
Rscript validation/reproduce_manuscript.R   [out]  # every headline number of the article
Rscript validation/rng_invariance.R                # RNG state invariance of the exact mode
Rscript validation/exact_invariance.R              # G04: exact-mode invariance across B/seed/N_range
Rscript validation/integration_check.R             # end-to-end run + reproducibility
Rscript validation/sweep_wald_wilson.R             # sweep of 28 configurations, Wald and Wilson
Rscript validation/cross_validation.R       <out>  # table2_crossvalidation.csv (Table S2.5 of Additional file 2)
Rscript validation/make_manuscript_assets.R <out>  # Figures 1-2 and the CSV of Table 1, the method comparison and the feature matrix

# Scripts added in 0.7.0 (each writes its log and tables to <out>)
Rscript validation/minimal_session.R        <out>  # the minimal session printed in the article
Rscript validation/comparator_calls.R                       <out>
Rscript validation/width_attainment_exact.R                 <out>
Rscript validation/mc_validate_buderer_exact_vs_reference.R <out>/curve_Se085_d007.csv
Rscript validation/negative_control_halfwidth.R             <out>
Rscript validation/edge_case_empty_arm.R                    <out>
Rscript validation/compare_wilson2022.R      <out> [dataS2_results.csv]
Rscript validation/compare_williamson2023.R  <out>
Rscript validation/exact_vs_mc_resolution.R  <out>
Rscript validation/bridge_frequentist_bayesian.R <out>
Rscript validation/success_set_shortcut.R    <out>
Rscript validation/applied_example.R         <out>      # about 25 minutes
Rscript validation/applied_example_snippet.R <out> <out of applied_example.R>/applied_example_main_table.csv
Rscript validation/operating_characteristics.R <out>
Rscript validation/misspecified_priors.R     <out>
Rscript validation/oc_montecarlo_crosscheck.R <out>  # Monte Carlo cross-check of operating_characteristics.R and misspecified_priors.R
Rscript validation/fixed_prevalence_control.R <out>  # sample size with the prevalence fixed at 0.20 (about 2 minutes)
Rscript validation/compat_v069_check.R <library> <out>  # twice (0.6.9 and 0.7.0), then: --compare <dir_v069> <dir_new> <log>
python -I validation/python/wilson2022_rule_enumeration.py <out> enum   # and, in a second run: mc

# Generators of the tables and of Figure 3 of the revised article, last:
Rscript validation/make_manuscript_tables_R1.R <out> --assets=<dir> --comparators=<dir> --width=<dir> \
  --negctl=<dir> --emptyarm=<dir> --wilson=<dir> --williamson=<dir> --bridge=<dir> --exactmc=<dir> \
  --shortcut=<dir> --applied=<dir> --oc=<dir> --mp=<dir> --crosscheck=<dir> [--witness=<wilson2022_rule_enumeration_enum.log>]
Rscript validation/make_figure3_applied.R <out> --applied=<dir>
```

The scripts added in 0.7.0 prepend the libraries named in the environment
variables `DTA_PKG_LIB` and `DTA_EXTRA_LIB` to `.libPaths()` when they are set,
except `comparator_calls.R`, `mc_validate_buderer_exact_vs_reference.R` and
`minimal_session.R`; those and the scripts of 0.6.9 use `R_LIBS`, and
`compat_v069_check.R` takes the library as an argument. The Python channel and the confirmatory scripts are
described in `python/README.md`.

## Where each output is

The console log of each script that was already in 0.6.9 is
`logs/<script>.log`; the files it writes, when it takes an output directory,
are in `logs/<script>/`. For the scripts added in 0.7.0 the output directory
is `logs/<script>/`: it contains the script's own log, its CSV tables and
`_console.log`, the console output as captured by the caller.
`logs/compat_v069_check/` holds the comparison with 0.6.9;
`logs/make_manuscript_tables_R1/` the tables of the revised article and
`logs/make_figure3_applied/` its Figure 3; `logs/wilson2022_rule_enumeration/`
the two logs of the rule enumeration. The Python channel writes
`python/logs/`, `python/run_all_top.log` and the versioned JSON files.
`logs/R_CMD_build.log`, `logs/R_CMD_check_00check.log`, `logs/testthat.Rout`,
`logs/testthat_NOT_CRAN.log` and `logs/sessionInfo.txt` are the build, check and
test logs of 0.7.0.

## Table and figure labels

Labels in scripts, logs and file names follow the numbering of the revised
article. Table 1 lists the functions; Table 2 gives the sample size according
to the planning question; Table 3 reproduces the results of Wilson et al. and
Williamson et al.; Table 4 is the functional comparison; Table 5 is the applied
example; Table 6 gives the exact operating characteristics. Figure 1 is the
exact probability of reaching the width with the Buderer formula; Figure 2 is
the joint assurance curve of the worked example; Figure 3 shows the curves of
the applied example. The tables moved to the Additional files are numbered S2.x
(Additional file 2) and S3.x (Additional file 3). A "Table n" followed by "of
Wilson et al." or "of Bulilete et al." refers to that article, not to this one.

The output files of `make_manuscript_tables_R1.R` are labelled by subsection
(S3.2a = S3.2.1, S3.2b = S3.2.2, S3.3b = S3.3.2), not by table:
`R1_AF3_S3.2a_williamson_figure3.csv` feeds Table S3.2b,
`R1_AF3_S3.2b_application_runs.csv` Table S3.2d and
`R1_AF3_S3.3b_cran_feature_matrix.csv` Table S3.3e. The two Monte Carlo
cross-check tables are named after their tables
(`R1_AF3_S3.9e_mc_crosscheck_success_coverage` and
`R1_AF3_S3.9f_mc_crosscheck_detection`).

Some file names keep the numbering that the scripts used when they were
written, which does not coincide with the numbering of the 0.6.9 article
either:

| File | Content | Table in the 0.6.9 article | Table in the revised article |
|---|---|---|---|
| `table1_functions.csv` | functions described in the 0.6.9 article | Table 1 | not used (Table 1 of the revised article is descriptive and is not generated by a script) |
| `table4_method_comparison.csv` | Buderer, BAM and joint Se/Sp under common assumptions | Table 2 | Table 2 (built by `make_manuscript_tables_R1.R`) |
| `table2_crossvalidation.csv` | cross-validation of core formulas | Table 3 | Table S2.5 (Additional file 2) |
| `table3_feature_matrix_source.csv`, `table3_feature_matrix.csv`, `table3_provenance.csv` | feature matrix against CRAN packages, 0.6.9 version | Table 4 | Table S3.3e (Additional file 3), built from `table_cran_feature_matrix_R1_source.csv`; Table 4 condenses the functional comparison |
| `table2_grid_points.csv` | the 63 grid points of Table 2 of Wilson et al. | not in the article | Additional file 3, Section S3.1 |

The registered files (`python/ref_bam.py`, `python/LOCKED_confirmatory_grid.json`,
`python/conf_*`, `python/dev_*`) and their logs contain no table or figure
numbers of the article.

The scripts inherited from 0.6.9 (`make_manuscript_assets.R`, `reproduce_manuscript.R`,
`rng_invariance.R` and `exact_invariance.R`) call the prior P2 (Sp Beta(18, 2)) "harmonized
priors"; the revised article calls it the "informative Sp prior". It is the same prior.

## Column "Data S2" of the comparison with Wilson et al.

`compare_wilson2022.R` takes as optional second argument the file
`dataS2_results.csv`, the results of running the supplementary code of Wilson
et al. (Data S2) as distributed. That file is not in this repository, because
the licence of the third-party code that produces it has not been verified. The
column `dataS2_original` of `compare_wilson2022_table.csv` and of
`suppl_C_cells.csv`, which holds the Data S2 results, therefore requires
obtaining Data S2 and running it as described in the header of the script;
without the file it is `NA`. The column that emulates the
programmed rule of Data S2 (`emulation`) is computed by the script without
Data S2 and reproduces the same values. The logs in this directory were
produced with a `dataS2_results.csv` obtained outside the repository. No script
of this repository copies, adapts or loads the Data S2 code.

## Scripts

- **Generators of the article's numbers, figures and tables**
  (`reproduce_manuscript.R`, `make_manuscript_assets.R`, `cross_validation.R`,
  `minimal_session.R`). All are run with `Rscript` against the installed package
  (dtasamplesize >= 0.7.0 for `reproduce_manuscript.R`). `make_manuscript_assets.R`
  writes Figures 1 and 2 (TIFF 600 dpi and vector PDF, 170 mm wide), the data of
  both curves, and the CSV files of Table 1, of the method comparison and of the
  feature matrix. Figure 1 is the exact probability that the full width of the
  95% interval for sensitivity is at most 0.14 against the number of diseased
  participants (Se = 0.85, d = 0.07, n = 40 to 200; truncated Wald and Wilson),
  computed with `mc_validate_buderer(method = "exact")`; Figure 2 is the exact
  joint assurance A(N) of the worked example, computed with
  `bam_assurance_curve()`. The script stops, without drawing, if an anchor fails
  (P(100) = 0.5683; first n with P >= 0.80 equal to 114 for both intervals;
  smallest n staying at or above 0.80 up to 200 equal to 117 for Wald and 114
  for Wilson; N = 678, A(678) = 0.8003489948 and A(677) = 0.7996848824).
  `table1_functions.csv` (Table 1 of the 0.6.9 article, not used in the revised
  article) is derived from `NAMESPACE` and the installed help pages; the feature
  matrix is read from the versioned `table3_feature_matrix_source.csv`, with a
  provenance record of the comparator package versions reviewed.
  `reproduce_manuscript.R` checks every headline number of the article and
  writes `reproduce_manuscript_checks.csv`. `minimal_session.R` is the minimal
  session printed in the article (about 45 seconds) and saves its console output
  in `minimal_session.log`. Each script prints the R and package versions and
  the output directory at the start, and a SUMMARY line with the result of every
  check at the end.
- **`cross_validation.R`** checks the Wilson/Wald intervals against base R
  `prop.test` and `Hmisc::binconf`, `buderer_n` against published values, and the
  Beta posterior interval width against numerical integration (a separate check
  verifies the Beta-Binomial predictive mass itself, by finite sum, simulation
  and the `lbeta`-based closed form; it is not one of the four published rows).
  The references are base R functions, `Hmisc` and direct numerical integration,
  run by the package authors; none is an external implementation of the Bayesian
  assurance method. The rows are written to `table2_crossvalidation.csv`, one
  row per check, with each result derived from the values the script computes
  (observed value, reference value, tolerance and a pass flag). The
  `package_version` column records the version of the package actually
  installed when the script runs. The script also runs the same style of check
  for three functions that remain in the package but are no longer described in
  the article (the Hanley-McNeil AUC variance against Monte Carlo, the
  net-benefit point estimate and standard error, and the imperfect-reference
  variance inflation factor and apparent-sensitivity closed form). Those checks
  print PASS or `**CHECK**` to the console as an internal sanity check on
  package code that ships, and are deliberately not written to the CSV. In the
  Hanley-McNeil check the 0.20 tolerance is applied only to the balanced
  configurations: that formula assumes a negative-exponential score model while
  the Monte Carlo reference is binormal, so the imbalanced configuration is
  reported as an expected deviation rather than as a failure. The `**CHECK**`
  mark in `logs/cross_validation.log` (line 6 in the 0.6.9 log, line 9 in this one) is documented and expected: at
  x = 19, n = 20 the package matches base `prop.test` while `Hmisc` differs at
  the upper limit. `Hmisc` is a required reference package; `pROC` is optional
  and loaded only if available.
- **`integration_check.R`** is a package-wide smoke test, not scoped to the
  article: it runs every exported estimator, prints the sample size each one
  returns, and verifies reproducibility under a fixed seed. It also asserts the
  "more uncertainty, larger N" ordering in a delimited section near the end: it
  computes `bam_sample_size(prior_se = p, method = "exact")` for the Se prior
  triple Beta(34, 6), Beta(17, 3), Beta(8.5, 1.5) and `stop()`s if the three
  `N_total` values are not strictly increasing.
- **`rng_invariance.R`** confirms that `bam_sample_size(method = "exact")` and
  `joint_sample_size()` save and restore `.Random.seed` and `RNGkind()` exactly,
  and that the exact-mode closed-form path returns identical results across RNG
  kinds.
- **`exact_invariance.R`** (G04) checks that the exact mode of `bam_sample_size()`
  is invariant to `B`, `seed` and `N_range`: for the published N = 678 case and
  the method-comparison N = 672 case it sweeps `B` in `{5, 10, 50, 5000}`
  crossed with `seed` in `{1, 2, 3, 4, 2026}` and `N_range` in `{NULL, an
  explicit range containing the crossing}`, and requires identical `N_total` and
  `joint_assurance` and identical `.Random.seed` and `RNGkind()` before and
  after every call. Run with `--version-check=off` against the 0.6.5 library it
  fails, which shows that it detects the defect it was written for (log in
  `historical_not_current/`).
- **`sweep_wald_wilson.R`** is the generator of the sweep of 28 configurations
  that compares the probability that the Wald and Wilson interval widths meet
  their target.
- **Comparison and evidence scripts of 0.7.0**: `comparator_calls.R` (calls to
  the comparator packages and to two integer-grid rules),
  `width_attainment_exact.R` and `mc_validate_buderer_exact_vs_reference.R`
  (exact probability that the realised width meets its target, and the Monte
  Carlo mode against it), `negative_control_halfwidth.R` (the width target read
  as a half-width or as a double width), `edge_case_empty_arm.R` (the
  convention for an arm with no participants), `compare_wilson2022.R` and
  `compare_williamson2023.R` (reproduction of the published results, with the
  difference and its explanation stated for each value),
  `exact_vs_mc_resolution.R`, `bridge_frequentist_bayesian.R`,
  `success_set_shortcut.R`, `applied_example.R` (with `applied_example_snippet.R`),
  `operating_characteristics.R` and `misspecified_priors.R`,
  `oc_montecarlo_crosscheck.R` (Monte Carlo cross-check of the two previous
  scripts) and `fixed_prevalence_control.R` (sample size with the prevalence
  fixed at 0.20). Each states in its
  header the conventions it uses (total or half width, one- or two-sided, level,
  joint or marginal objective, design or analysis prior, success rule) and the
  sources of the published figures it compares with. `compat_v069_check.R` shows
  that the default arguments of 0.7.0 give results identical to 0.6.9; its log
  is `logs/compat_v069_check/compat_v069_check.log`.
- **`python/wilson2022_rule_enumeration.py`** enumerates, for the examples of
  Wilson et al., the success set under three rules (the definition in the
  article, the rule as printed, and the rule as programmed in Data S2) and, in a
  second part (`mc`), simulates the generative model directly. It is a check by
  the same team, written with SciPy; it does not run or load the Data S2 code.
  Its `enum` log is the input of the table of witness cases of Additional file 3
  (`make_manuscript_tables_R1.R --witness=`).
- **`make_manuscript_tables_R1.R` and `make_figure3_applied.R`** generate the
  tables and Figure 3 of the revised article from the outputs of the scripts
  above. Every number that the text quotes is checked against the expected value
  (an "anchor"); the script prints PASS or `**CHECK**` for each and exits with a
  non-zero code if one fails.
- **`table3_feature_matrix_source.csv`** is the hand-curated data behind the
  0.6.9 feature matrix; `table_cran_feature_matrix_R1_source.csv` is the one of
  the revised article: the matrix whose cells were checked in November 2026
  against the help pages of the installed versions (`tools::Rd_db()`, R 4.6.2;
  presize 0.3.11, pROC 1.19.1, MKpower 1.1, epiR 3.0.1, MKmisc 2.0,
  SampleSizeDiagnostics 0.1.1), identical to Table S3.3e. The generator records
  those versions in `R1_AF3_S3.3b_cran_feature_matrix_provenance.csv`, next to
  the versions reviewed in August 2026, and an anchor (S3.3b-v) checks them
  against the header of the matrix. They are kept as separate, versioned files so the
  curation is visible and diffable on its own.

## Package check logs (0.7.0, current)

- `logs/R_CMD_check_00check.log`: `R CMD check --as-cran --no-manual` of the
  0.7.0 source tarball, run in a folder on the local system drive with
  `NOT_CRAN=true` and the CRAN incoming checks that need the network switched
  off; result `Status: OK` (no errors, warnings or notes). The build is in
  `logs/R_CMD_build.log`.
- `logs/testthat.Rout`: the test suite run inside that check
  (`[ FAIL 0 | WARN 0 | SKIP 3 | PASS 2141 ]`).
- `logs/testthat_NOT_CRAN.log`: the complete suite run from the source tree with
  `NOT_CRAN=true` (`[ FAIL 0 | WARN 0 | SKIP 2 | PASS 2147 ]`). It has one skip
  fewer and six more expectations than the run inside the check because the test
  that needs the `man/` folder runs here and is skipped inside the check.
- `historical_not_current/R_CMD_check_with_pdf_manual_00check_v0.6.9.log`: an
  `R CMD check --as-cran` of 0.6.9 that also attempted the PDF manual. It ended
  with 1 ERROR and 1 WARNING because the local TeX installation lacks the
  Courier font metrics (`pcrr8t`). It cannot be regenerated on this machine
  (`makeindex` is missing) and is not evidence for 0.7.0.

## A note on what these scripts prove

`reproduce_manuscript.R` and `integration_check.R` are reproducibility checks:
they confirm that the installed package, run today, prints the same numbers the
article reports (or that it runs end to end, self-consistently, under a fixed
seed). A PASS there is not a validity check: a defect that biases a result
reproduces exactly as cleanly as a result that is correct, and the numbers
`reproduce_manuscript.R` checks were themselves revised twice across earlier
package versions after those versions passed their own reproducibility checks
against defective numbers (see `NEWS.md`). `cross_validation.R` is closer to a
genuine validity check for the specific formulas it covers, since it compares the
package's output against a separate computation (base R, `Hmisc`, or direct
numerical integration, all run by the package authors) rather than against the
package's own past output, but it only covers the four quantities in its table.
The Python channel and `wilson2022_rule_enumeration.py` are checks by the same
team. Validity beyond what these scripts test is established separately, release
by release, by the package's own documentation and by external review of the
methods, never by a script agreeing with itself.

## A note on data

The package ships no dataset, and the article analyses no empirical data. Every
quantity is computed from samples drawn inside the package's functions from the
stated parameters (sensitivity, specificity, prevalence, and the Monte Carlo
size `B`) under a fixed seed, so the parameters together with the seed determine
the data exactly. The applied example uses the counts published by the sources
named in the header of `applied_example.R`.

If a materialised copy of the simulated samples is needed, set `DUMP_DATA <- TRUE`
at the top of `reproduce_manuscript.R`; it writes them as CSV to
`simulated_data/` inside its output directory. Those files are a convenience
only: the script, not the CSV, is the authoritative record of how the data arise.
