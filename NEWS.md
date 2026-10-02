# dtasamplesize (development version)

**Test robustness across R versions (R-devel/R >= 4.6.0 RNG-stream
change).** R-devel (tested: R 4.6.0 r90605, 2026-09-30) changed the RNG
stream consumed internally by some sampling primitives relative to
R 4.5.x, so several Monte Carlo-derived test figures that were pinned to
an exact decimal value captured under R 4.5.x (e.g.
`test-joint_contract.R`'s `joint_prob_se_sp == 0.8041`,
`joint_prob_mcse == 0.002806`) failed on R-devel for the *same*,
unaffected `n_total`/grid outcome -- only the continuous Monte Carlo
estimate at that `N` shifted, by less than its own reported MCSE. Every
such assertion across the suite (`test-joint_contract.R`,
`test-rng_state.R`, `test-ss_imperfect_ref.R`, `test-ss_unified.R`) now
checks the value against a tolerance derived from its own reported Monte
Carlo standard error (typically 3x MCSE), with the exact historical
figure kept as an additional check guarded by
`if (getRversion() < "4.6.0")` so the published-manuscript numbers remain
pinned on the R series they were produced on. Deterministic exact-mode
tests (`bam_sample_size(method = "exact")`'s finite sum,
`joint_sample_size()`'s grid/AUC-gate structure, `nb_assurance_ceiling()`'s
closed-form quadrature) are unaffected and were left untouched -- they do
not depend on R's RNG at all. No published figure changes on the R series
it was computed on.

**Non-crossing contract extended to `ss_time_dependent_roc()`.** Applies
the same non-crossing rule introduced for `bam_sample_size()`/
`joint_sample_size()` in 0.6.6 (and already respected by
`ss_unified()`/`ss_net_benefit()`, reviewed and confirmed compliant in
this release): if no candidate `N` in `N_range` reaches `target_prob` for
a given `censoring_rates` entry, `ss_time_dependent_roc()` versions
`<= 0.6.6` set `N_required` to `max(N_range)` and reported the (failing)
probability observed there in `results`, with no field distinguishing
that outcome from a genuine solution short of reading the warning text --
and if ANY entry failed to converge, the top-level `n_total` (the max
over `results$N_required`) was itself silently pinned to that
`max(N_range)` value rather than reflecting that the true worst case was
unknown. Confirmed against the installed 0.6.6
(`CORRECCION_INTEGRAL_2026-09-30/tools/L067_tests_contra_066` log): a
deliberately impossible target (`delta_auc = 0.001`, `target_prob =
0.999`, `N_range = seq(100, 140, by = 20)`) returned `n_total = 140` with
`prob_achieved = 0` for every censoring rate, as if `N = 140` were a
validated design. `ss_time_dependent_roc()` now sets, per censoring rate,
`N_required = NA_integer_` and `prob_achieved = NA_real_` (never
`max(N_range)`) whenever that rate's search does not cross `target_prob`,
and reports the new `target_reached`, `max_prob_evaluated` and
`N_at_max_prob` diagnostic columns in `results`; the top-level `n_total`/
`n_diseased` are `NA` and a new top-level `target_reached` field is
`FALSE` whenever ANY censoring rate failed to converge, mirroring
`bam_sample_size()`/`joint_sample_size()`/`ss_net_benefit()`'s
`N_conservative`. See `?ss_time_dependent_roc`, `@return`.
`ss_unified()`, `ss_net_benefit()`, `ss_imperfect_ref()` and
`ss_adaptive_prevalence()` were reviewed against the same non-crossing
contract in this release: `ss_unified()` already reports `status =
"unreachable"`/`"grid_exhausted"` with `N_effective = NA_integer_`;
`ss_net_benefit()` already reports `N_conservative = NA_integer_` on a
non-crossing or structurally-infeasible threshold; `ss_imperfect_ref()`
and `ss_adaptive_prevalence()` compute sample size in closed form (no
`N_range` grid search), so the non-crossing failure mode does not apply
to them. No published figure changes.

**Sole package authorship.** `Authors@R` now lists only Jesús D.
Rojas (`aut`, `cre`, ORCID 0000-0003-0912-490X); Rafael
Pichardo-Rodriguez is no longer listed as an author of this package (he
was not assigned any other role). The maintainer email
(`jrojasrivero@gmail.com`) is unchanged -- that is a separate decision,
not made here. `LICENSE`/`LICENSE.md` already named the generic
"dtasamplesize authors" as copyright holder and needed no change; a
search of `README.md`, the vignettes and `inst/` found no other mention
of Rafael Pichardo-Rodriguez as a package author. Added `.zenodo.json`
(title "dtasamplesize", creator Jesús D. Rojas with the same ORCID,
`upload_type: software`, `license: MIT`) for future Zenodo deposits;
listed in `.Rbuildignore`. No published figure, API, or behavior changes.

# dtasamplesize 0.6.6

**M-01 correction: `bam_sample_size()`'s exact-mode joint search.**
Responds to the integral audit of the BMC MRM submission package
(`AUDITORIA_INTEGRAL_2026-09-29`), finding M-01. Two defects in the
"headline" joint search (`N_total` / `joint_assurance`) under
`method = "exact"`, the package's documented deterministic, B/seed-free
calculation:

* **Spurious B/seed dependence when `N_range = NULL`.** Versions `<= 0.6.5`
  built the automatic search grid (`N_range`) from `n_se`, `n_sp` and
  `N_total_P90` -- all three Monte Carlo quantities depending on `B` and
  `seed` -- so the exact calculation could silently return a DIFFERENT
  `N_total` for the SAME deterministic priors and targets, purely because
  `B` or `seed` changed the grid it happened to search. Confirmed by
  `CORRECCION_INTEGRAL_2026-09-30/tools/L01_reproducir_defecto_065.R`
  against 0.6.5: with the article's worked-example priors
  (`prior_se = c(17, 3)`, `prior_sp = c(2, 2)`, `prior_prev = c(4, 16)`,
  `delta_se = 0.14`, `delta_sp = 0.10`, `target_assurance = 0.80`),
  `B = 5, seed = 1` returned `N_total = 631` (joint assurance 0.764813,
  below target) while `B = 5000, seed = 2026` (and most other `B`/`seed`
  combinations) returned the correct `N_total = 678`. `N_range = NULL`
  under `method = "exact"` now scans every integer `N` starting at 2,
  ascending, in doubling blocks capped at a new `N_max` argument (default
  3000), independent of `n_se`, `n_sp`, `N_total_P90`, `B` and `seed`; see
  `?bam_sample_size`. The exact per-arm cache is extended block by block
  (`.bam_exact_extend_cache()`), not rebuilt from scratch, so growing the
  search stays cheap.
* **Non-crossing searches returned `max(N_range)` as a solution.** Under
  both `method = "exact"` and `method = "monte_carlo"`, if no candidate `N`
  reached `target_assurance`, versions `<= 0.6.5` set `N_total` to
  `max(N_range)` and returned a real (if below-target) `joint_assurance`,
  with no field distinguishing that outcome from a genuine solution short
  of inspecting the warning text. Confirmed by the same script: an
  insufficient `N_range = 100:200` returned `N_total = 200` with
  `joint_assurance = 0.0112`, as if `N = 200` were a validated design.
  `bam_sample_size()` now sets `N_total = NA_integer_`, `n_total =
  NA_integer_`, `joint_assurance = NA_real_` and the new field
  `target_reached = FALSE` whenever no candidate reaches the target, under
  either method, together with the new diagnostic fields
  `max_assurance_evaluated` and `N_at_max_assurance` (the highest
  assurance actually seen, and at which `N`) and a warning naming the
  ceiling reached and how to search further. `print.dtasamplesize()` now
  reports this case explicitly ("Target assurance NOT reached -- no sample
  size returned") instead of printing `N_total: NA` under the normal,
  converged-looking label.
* **New fields, always present:** `target_reached` (logical),
  `N_range_used` (the candidate `N` actually evaluated, in order),
  `search_type` (`"integer_scan_auto"`, `"user_N_range"`, or
  `"monte_carlo_auto"`), and `N_max` (the ceiling in effect for the
  automatic integer scan; `NA_integer_` otherwise). New argument: `N_max`
  (default `3000L`), the ceiling for the automatic exact-mode integer
  scan.
* **Explicit `N_range` under `method = "exact"`** is now always searched
  as `sort(unique(as.integer(N_range)))`, ascending, regardless of the
  order or duplication the caller supplied (this was already the effective
  behavior for a strictly ascending, duplicate-free `N_range`, so no
  published result changes; see below).
* **No published headline number changes.** All three published exact-mode
  results were re-verified under the corrected code, including against the
  full `B` in `{5, 10, 50, 5000}` x `seed` in `{1, 2, 3, 4, 2026}` battery
  that exposed the defect under 0.6.5:
  `CORRECCION_INTEGRAL_2026-09-30/tools/L01_verificar_cifras_066.R`
  confirms, byte-identically across every `B`/`seed` combination,
  N = 678 (joint assurance 0.8003489948, prior_prev = c(4, 16)), N = 672
  (joint assurance 0.8002692084, harmonized priors), and the package
  default N = 583 (`prior_prev = c(6, 14)`, the formal default). Minimality
  was independently re-verified via the internal exact Beta-Binomial
  helpers: A(677) = 0.7996848824 < 0.80 <= A(678) = 0.8003489948, and no
  integer in 2..677 reaches 0.80 either.
* **`joint_sample_size()` (Monte Carlo, unrelated search code): M-03
  correction, Option A (lote 01b).** `CORRECCION_INTEGRAL_2026-09-30/
  DECISION_JOINT_SAMPLE_SIZE.md` analyzed this function's non-crossing
  behavior and its "smallest N" wording (P-DISEÑO) and presented two
  options; the author approved **Option A** (2026-10-01): keep the
  step-10 default grid (`N_range = seq(100, 800, by = 10)`) and its
  published figure unchanged, reinterpret `N = 580` as "the first
  candidate of the step-10 grid reaching the target" rather than "the
  smallest sample size" (it is not: a finer, integer-step search crosses
  the same target as early as N = 578-579 at the manuscript's operating
  point), plus the output changes common to both options:
    - **Non-crossing searches no longer return `max(N_range)` as a
      solution** (the same class of defect M-01 corrected in
      `bam_sample_size()` above). When no candidate in `N_range` reaches
      `target_prob` for Se/Sp, `joint_sample_size()` now sets `n_total =
      NA_integer_`, `target_reached = FALSE`, `joint_prob_se_sp =
      NA_real_` (and its alias `joint_prob`), `n_diseased = NA_integer_`,
      `n_non_diseased = NA_integer_`, and reports `max_joint_prob_evaluated`
      / `N_at_max_joint_prob` as diagnostics, together with a warning
      distinguishing "the AUC gate blocked every candidate" from "the AUC
      gate passed somewhere but Se/Sp never reached target_prob".
    - **New fields, always present:** `target_reached`, `N_range_used`,
      `search_type` (always `"grid_first_candidate"`), `seed` (previously
      not echoed back), `joint_prob_mcse` (the Monte Carlo standard error
      of `joint_prob_se_sp` at `n_total`), `auc_gate` (a list with `AUC`,
      `delta_auc`, a `table` of the deterministic Hanley-McNeil gate's
      outcome at every candidate N -- `NA`/`FALSE`/`TRUE` for
      skipped/blocked/passed, never a probability of 0 standing in for
      "blocked" -- and `first_N_auc_pass`), `max_joint_prob_evaluated`,
      and `N_at_max_joint_prob`. `auc_gate_passed` is kept for backward
      compatibility.
    - **The search itself is unchanged**: same candidate order, same
      per-candidate `set.seed()`, same `joint_prob >= target_prob` rule, so
      **no published number changes**. The published N = 580 (joint
      probability 0.8041, now also reported with `joint_prob_mcse` =
      0.002806 under `B = 20000`, `seed = 2026`, the values used for that
      figure -- the function's own `B` default remains 5000) is
      byte-identical before and after this change
      (`tools/L01b_verificar_cifra.R`).
    - **`print.dtasamplesize()`** now reports a non-crossing
      `joint_sample_size()` result explicitly ("Target NOT reached -- no
      sample size returned") instead of printing `N_total: NA` under the
      normal, converged-looking label, and a successful result now shows
      `joint_prob_se_sp`, an explicit "first candidate in N_range"
      statement (never "smallest"), `joint_prob_mcse`, and the first N
      passing the AUC gate.
    - **RNG preservation was already correct** (`save_rng_state()` /
      `restore_rng_state()`, unchanged by this release) and is now also
      pinned by a dedicated test on the non-crossing branch.
  See `CORRECCION_INTEGRAL_2026-09-30/analisis_joint/EVIDENCIA_P-DISENO.md`
  for the original analysis.
* **Tests.** `tests/testthat/test-joint_contract.R` (new) locks in the
  published-case figure (580 / 0.8041 / MCSE 0.002806) together with
  `search_type`, `N_range_used`, `seed`, `auc_gate`, both non-crossing
  branches (AUC-gate-blocked-everywhere vs. AUC-gate-passed-but-below-
  target), the AUC-gate-blocked-is-not-probability-0 distinction, RNG
  preservation on the non-crossing branch, and that `print()` never shows
  an N as a solution when `target_reached = FALSE` nor the word
  "smallest" when it is `TRUE`. `tests/testthat/test-joint_sample_size.R`'s
  tests that used a single-value `N_range` purely as a Monte Carlo
  probability calculator (not a search) now pass an explicit, very low
  `target_prob` so the single candidate is still accepted and the real
  computed probability is reported (`joint_prob_se_sp`) rather than `NA`;
  the one case where the true probability is exactly 0 (and so can never
  cross any positive `target_prob`) now reads `max_joint_prob_evaluated`
  instead, which is exactly the value that field used to be reported as
  before this release. The "expected margin of exactly 0" test is updated
  to probe `target_reached` / `auc_gate` instead of the removed
  `max(N_range)` fallback it used to rely on.
* **Tests (unrelated to the above).** `tests/testthat/test-bam_exact_contract.R` (new) locks in
  B/seed invariance (reduced battery inline, full battery under
  `skip_on_cran()`), the explicit-`N_range` sort/unique/ascending
  contract, the shared no-crossing contract for both methods, RNG
  preservation under the integer scan (including on non-crossing),
  minimality of the returned N* via the internal exact helpers, the N = 1
  (always degenerate) vs. N = 2 (arm of size 1, not degenerate) boundary,
  `.bam_exact_extend_cache()` block-extension correctness, and
  `print.dtasamplesize()` on a non-crossing result.
  `tests/testthat/test-bam_sample_size.R`'s tests that exercised a
  single-value or otherwise non-crossing `N_range` and read the resulting
  `joint_assurance` directly (`res_591`, `res_677`, the N = 1 degenerate
  test, the H-06 arm-of-size-1 test, and the N in {10, 20, 30} degenerate
  exact-mode test) are updated to read `max_assurance_evaluated` /
  `N_at_max_assurance` / `target_reached` instead: those tests codified the
  defective "return max(N_range) as a solution" behavior this release
  corrects, and their actual diagnostic purpose (verifying the assurance
  value computed AT a given N) is preserved under the new field names.

# dtasamplesize 0.6.5

**Blind-audit correction release.** Responds to a blind audit of the
BMC Medical Research Methodology software-article submission package
(`AUDITORIA_CIEGA_2026-09-27/INFORME.md`, findings H-02, H-03, H-05, H-06,
H-08, H-10). Only two published numbers change, and both are outputs of
validation scripts, not of the package itself: the textual Wald-sweep range
in the manuscript, and the reversal count reported for
`conf_PY_06_monotone_in_N.py`. Every package-computed headline number is
unchanged: N = 678 (0.8003489948), N = 672 (0.8002692084), N = 580 (0.8041),
and the single-case Wald/Wilson coverage 0.5730 (Se = 0.85, d = 0.07).

## Validation scripts

* **H-02.** Added `validation/sweep_wald_wilson.R`, the previously missing
  generator for the manuscript's "systematic sweep of 28 configurations"
  claim (Se in {0.60, 0.70, 0.75, 0.80, 0.85, 0.90, 0.98} x d in
  {0.03, 0.05, 0.07, 0.10}, B = 200000, seed 2026), with the single
  previously published case (Se = 0.85, d = 0.07, P(width <= target) =
  0.5730) reproduced inside the grid as an internal control. Output
  archived to `validation/logs/sweep_wald_wilson.log`. The reproduced
  ranges (Wald 0.469-0.872, concentrated 0.51-0.56; Wilson 0.000 at
  Se = 0.98 for d = 0.05/0.07/0.10, up to 1.000 at Se = 0.60, d = 0.10)
  differ from the previously published, undeclared-grid ranges (Wald
  0.376-0.851, concentrated 0.47-0.61); the manuscript text is corrected to
  the reproduced grid and cites the new script. The qualitative conclusion
  is unchanged: the ~50% coverage is characteristic of the Wald interval,
  and the Wilson interval spans the full [0, 1] range depending on
  configuration.
* **H-03.** `validation/python/conf_PY_06_monotone_in_N.py` now uses a
  declared tolerance of 1e-9 (was an undeclared `d < -1e-15`, which counted
  floating-point rounding noise near P approx 1 as reversals) and clips P to
  [0, 1] before comparison. For m <= 800, this finds **73** true reversals
  (`|Delta| > 1e-9`), not the previously published 395 - the same 73 recur
  for any tolerance between 1e-9 and 1e-4, and for m <= 900. The worked
  example `p(26) = 0.068288 -> p(27) = 0.063834` is unchanged and still
  correct. `validation/python/HASHES_2026-09-28.txt` documents the change to
  the hash-locked script with a dated note; `validation/python/README.md`
  and `report.html` are updated to the new count.
* **H-08.** `validation/make_manuscript_assets.R` now exports Figure 1 and
  Figure 2 at a width of 170 mm (was a narrower default), keeping aspect
  ratio, TIFF resolution (>= 600 dpi), and font legibility unchanged.
* **H-05.** Re-verified `validation/python/conf_R_07_reverify_version.R`
  (and its accompanying reverification scripts) against the tarball built
  from this release and installed to a private library; `requirements.txt`
  updated from `dtasamplesize 0.6.3` to `dtasamplesize 0.6.5`. See
  `validation/python/logs/` for the archived log.
* **H-10.** `R CMD build` + `R CMD check --as-cran --no-manual` re-run from
  a neutral directory outside this repository (no `audit`/`_trabajo`/internal
  path names). See `validation/logs/R_CMD_check_00check.log` and the
  testthat log for the exact `Status:` line and expectation/skip counts.

## Bug fixes

* Unified the degenerate-replicate rule between `joint_sample_size()` and
  `bam_sample_size()`: both now treat an arm as degenerate (forced failure)
  only when it received zero subjects (n = 0), not n < 2.
  `joint_sample_size()`'s Wilson interval is well defined at n = 1 (width
  approx 0.79). `bam_sample_size()` already used n = 0 and is unchanged. No
  published number changes (Table 4: N = 580, 0.8041; N = 672; N = 678).
  Rd pages for both functions updated to document the shared rule; new
  tests in `tests/testthat/test-joint_sample_size.R` and
  `tests/testthat/test-bam_sample_size.R` pin n = 0, 1 and 2 in each arm.

# dtasamplesize 0.6.4

**Documentation and validation-tooling release.** No computed result
changes: every quantity reported in the manuscript and reproduced by
`validation/reproduce_manuscript.R` is identical to 0.6.3.

## Documentation

* `bam_sample_size()`'s `@details` (and the corresponding `?bam_sample_size`
  help page) no longer describe the two worked examples' `prior_prev =
  c(4, 16)` as "the defaults": that prior is the accompanying article's
  worked-example prevalence prior, not `bam_sample_size()`'s actual default
  of `prior_prev = c(6, 14)`. The formal `@param prior_prev` default
  documentation was already correct; only the `@details` prose was
  overstated.
* `ss_imperfect_ref()`'s worked-example section heading now reads "Worked
  example (package defaults, except `prev = 0.20` instead of the default
  0.30)", rather than the previous "(package defaults with `prev = 0.20`)",
  which could be misread as `prev = 0.20` also being a package default.
* `validation/reproduce_manuscript.R` and `validation/make_manuscript_assets.R`
  no longer label the article's worked-example priors (`prior_prev =
  c(4, 16)`, etc.) as package "defaults" in comments; the object previously
  named `bam_default` is now `bam_example`. Figures 1 and 2's y-axis labels
  are now drawn horizontal (`las = 1`) rather than rotated.
* `validation/table3_feature_matrix_source.csv` now marks the package's
  less mature modules (joint Se+Sp+AUC search, imperfect-reference
  correction, adaptive prevalence re-estimation, time-dependent ROC, net
  benefit, and the unified framework) as "Experimental" rather than "Yes",
  distinguishing them from the mature CI-width sizing capability.

## Tests

* New test pins the documented default `prior_prev = c(6, 14)` and the
  `N_total = 583` it produces under `bam_sample_size(method = "exact")`
  with otherwise-default arguments, guarding against a repeat of the
  defaults/worked-example mismatch described above.

## Validation scripts

* `validation/integration_check.R` now asserts the "more uncertainty ->
  larger N" ordering as executable code (an "Ordering check" section
  computing `bam_sample_size(prior_se = p, method = "exact")` for the Se
  priors Beta(34,6), Beta(17,3), Beta(8.5,1.5) and `stop()`ping if
  `N_total` is not strictly increasing), rather than only asserting it in
  prose as before; `validation/README.md` is updated to match.
* `validation/cross_validation.R` now applies the 0.20 Hanley-McNeil
  tolerance only to the balanced AUC configurations; the imbalanced
  configuration is reported as an expected deviation (the Hanley-McNeil
  formula assumes a negative-exponential score model, the Monte Carlo
  reference here is binormal, and the gap between the two widens as the
  groups become imbalanced). `validation/README.md` documents the same
  criterion.
* New `validation/rng_invariance.R` regenerates the article's published
  numbers under five selected RNG kinds (Mersenne-Twister, L'Ecuyer-CMRG,
  Wichmann-Hill, Marsaglia-Multicarry, and Knuth-TAOCP-2002), as evidence
  for the manuscript's RNG-invariance claim; its output is logged under
  `validation/logs/`.
* New `validation/python/` directory: an independent Python reference
  implementation, a locked confirmatory grid, and recorded hashes and logs
  used to cross-check the R package's results against a second,
  independently written implementation.

## Compatibility

* `ss_time_dependent_roc()` no longer lets an upstream `timeROC` warning
  introduced in R-devel ("object length is not a multiple of subscript
  length") escape; only that exact message is muffled.

## Infrastructure

* New `.github/workflows/R-CMD-check.yaml` runs `R CMD check` on Windows,
  macOS, and Ubuntu on push and pull request.

# dtasamplesize 0.6.3

Corrects four defects reported by an independent review of 0.6.2.
Published sample sizes are unchanged.

* The net-benefit feasibility ceiling was estimated by Monte Carlo with a
  default of 200,000,000 draws, which exhausted memory on a machine with
  ordinary RAM. It is now computed by deterministic integration: it agrees
  with an independent quadrature to six decimal places, is identical across
  seeds, and needs no large allocation. A memory budget now rejects an
  oversized request instead of attempting it.
* `ss_unified()` no longer reports a sample size when the search does not
  converge. `N_effective` and `n_total` are `NA` in that case, and the
  returned object carries an explicit `status` field (`converged`,
  `unreachable`, `grid_exhausted`).
* README, DESCRIPTION, the reproduction script and the accompanying article
  agree on the version; `reproduce_manuscript.R` requires the exact version
  whose numbers it checks.
* Four tests asserted the absence of warnings at a `B` the package now warns
  about, or matched the small-B warning where the ceiling message was
  intended.

# dtasamplesize 0.6.2

**Validity release.** Continued adversarial verification of `ss_unified()`
found nine further defects: how its three `decision` rules handle a
candidate-N grid that is not already sorted ascending or that contains
repeated values, an unvalidated `N_range`, a small-`B` gap in the isotonic
margin's block-pooling, an unsupported precision claim for the `check_nb`
ceiling, and a gap in 0.6.1's own RNG-generator fix that left the *normal*
generator (as opposed to the uniform one) still inherited from the caller.
None of this release's fixes move the manuscript's headline sample sizes
computed under the package's own default RNG configuration, realistic `B`,
and `N_range` given already sorted ascending -- the calling pattern used
throughout `validation/make_manuscript_assets.R` and
`validation/reproduce_manuscript.R` -- for every affected function.

## Validity fixes (results change for affected callers)

* **`ss_unified()`'s `decision \%in\% c("point", "lower_bound")` accepted
  whichever N reached `target_assurance` FIRST IN `N_range`'s OWN ORDER,
  not the smallest N that did.** The search loop evaluated `N_range` in
  the order given and stopped (or, under `full_grid = TRUE`, still
  recorded only the first acceptance) at the first N to pass -- correct
  only when `N_range` happens to already be ascending. Reproduced directly
  (identical priors, `B`, `seed`, and *set* of 17 candidate N throughout):
  the same search returned `N_effective = 850` (the true minimum) with
  `N_range` given ascending, `1200` (+41\%) given descending, and a third,
  different value given shuffled -- with no warning of any kind. Both
  rules now traverse the DISTINCT candidate N in ascending order
  internally (via the same value-keyed RNG-stream lookup already used to
  make the simulated data itself order-invariant), regardless of how
  `N_range` was arranged, so the smallest N reaching `target_assurance` is
  always the one returned. `decision = "isotonic"` (the default since
  0.5.1) already evaluated the complete range and sorted before fitting,
  so it was already invariant to `N_range`'s order and is unaffected. See
  `?ss_unified`, `@param full_grid`.

* **`ss_unified()`'s "did not converge" fallback reported the assurance of
  whichever N the search loop happened to evaluate LAST, not of the N it
  actually returns (`max(N_range)`).** When no N in `N_range` reached
  `target_assurance`, `optimal_N` is set to `max(N_range)`, but
  `joint_assurance`/`assurance_mcse`/`assurance_lower` were carried over
  from local variables left by the search loop's final iteration -- the
  last element of `N_range` IN THE ORDER GIVEN, not necessarily its
  largest value. Reproduced directly (identical non-converging search,
  same seed/`B`/priors/*set* of candidate N, `decision = "isotonic"`): the
  same returned `N_effective = 400` was reported with `joint_assurance =
  0.37850` given `N_range` ascending, `0.06400` given descending (the true
  assurance at N = 200, the smallest candidate, not at N = 400), and
  `0.30000` given shuffled -- three different numbers attached to the same
  N, none of them wrong in isolation but only one of them actually
  describing the N returned. The reported quantities are now always looked
  up from `grid_results` AT `N = max(N_range)` specifically, for every
  `decision` rule, regardless of the order `N_range` was searched in.

* **`ss_unified()`'s `decision = "isotonic"` let a duplicated N in
  `N_range` fabricate precision it did not have.** A repeated N always
  replays the identical L'Ecuyer-CMRG stream (by design, so that asking
  about the same design twice gives the same answer), so its assurance
  estimate is byte-identical to the first copy's -- zero new information.
  The margin's block-pooling step, however, counted each repeated N as an
  additional independent observation when computing a pooled block's
  effective sample size (`B * block_len`), so replicating grid points
  fabricated statistical precision from nothing. Reproduced directly (a
  16-point base grid, otherwise identical `B = 4000`, seed, and priors,
  each point replicated 1x/2x/3x/5x/10x): `N_effective` fell from 857 to
  835 as the replication count rose, with `assurance_lower` correspondingly
  -- and wrongly -- INCREASING (0.800163 to 0.800275), purely from
  re-counting the same evidence. The grid is now deduplicated by N value
  before the isotonic fit is computed, so a repeated N changes nothing
  about the result (as it should not, carrying no new evidence), and only
  genuinely distinct N contribute to a pooled block's effective sample
  size.

* **`ss_unified()`'s `decision = "isotonic"` margin could report a
  deceptively strong lower bound from a minimal `B`.** The block-pooling
  step (previous item) cannot distinguish a block that `isoreg()`
  genuinely pooled to enforce monotonicity from a run of already-distinct,
  non-duplicated N whose RAW assurance simply happened to coincide -- both
  look like an identical run of fitted values. At a realistic `B` this
  distinction rarely matters (raw assurance is close to continuous, so
  coincidental ties are negligible), but at a very small `B` raw assurance
  can only take `B + 1` possible values -- at `B = 1`, just `{0, 1}` -- so
  long flat runs of DISTINCT N arise routinely by chance, not because the
  underlying curve is genuinely flat there. Reproduced directly: a
  generous scenario at `B = 1`, one pass/fail replicate per N and no
  duplicated N anywhere in a wide `N_range`, pooled dozens of consecutive,
  genuinely distinct N into single blocks and reported `assurance_lower`
  near 0.81 -- built from exactly one bit of information per contributing
  N, with no non-convergence warning. The pooling multiplier is now capped
  at `min(block_len, B)`: at `B = 1` no block's effective sample size can
  exceed `1 * 1 = 1`, so `assurance_lower` cannot exceed roughly 0.27
  regardless of how wide a spurious flat run happens to be. The cap only
  bites when `block_len` exceeds `B`, which requires more grid points
  pooled together than `B` itself -- every scenario in this package's own
  tests and the manuscript uses a `B` in the hundreds to tens of thousands
  against grids of well under 150 points, so the cap is a complete no-op
  there.

* **A gap in 0.6.1's own RNG-generator fix: `set.seed(..., kind =
  "Mersenne-Twister")` fixes the UNIFORM generator but not the NORMAL or
  discrete-sampling ones.** `RNGkind()` selects three independent
  generators (uniform, normal, discrete-sampling); naming only `kind`
  leaves `normal.kind`/`sample.kind` inherited from the caller, exactly
  the class of defect 0.6.1 set out to close. `ss_time_dependent_roc()` is
  the only function in the package that draws from the normal generator
  (`stats::rnorm()`, for the simulated biomarker) and was therefore the
  only one practically affected: reproduced directly, identical arguments
  and seed gave `n_total = 420` under the caller's `normal.kind =
  "Inversion"` (R's default) but `400` under `"Box-Muller"`, with no
  warning. `ss_imperfect_ref()` has the same gap independently of this
  release's other fixes (it predates 0.6.1's RNG work and was not part of
  it): reproduced directly with the exact call used in
  `validation/reproduce_manuscript.R` (`B = 6000`, `seed = 2026`), the
  manuscript's published apparent sensitivity of 0.764 (bias -0.086)
  reproduces exactly under Mersenne-Twister but rounds to 0.765 under
  Wichmann-Hill. Every `set.seed()` call that seeds a generator in this
  package now names `kind`, `normal.kind` and `sample.kind` explicitly
  (`"Mersenne-Twister"`/`"Inversion"`/`"Rejection"`, or, for
  `ss_unified()`'s own per-N streams, `"L'Ecuyer-CMRG"`/`"Inversion"`/
  `"Rejection"`), including in the five functions that do not themselves
  draw from the normal or discrete-sampling generators (`bam_sample_size()`,
  `joint_sample_size()`, `mc_validate_buderer()`,
  `ss_adaptive_prevalence()`, `ss_net_benefit()`, and the internal
  `nb_assurance_ceiling()` and `ss_unified()` search itself) -- a no-op for
  their own results, but it makes the reproducibility guarantee general
  rather than dependent on which generator a given function happens to
  consume. `ss_imperfect_ref()`'s RNG-state restoration is also upgraded
  from restoring only `.Random.seed`'s value (its previous, pre-0.6.1-style
  mechanism) to the full `save_rng_state()`/`restore_rng_state()` used
  elsewhere, so the caller's complete RNG configuration -- not just the
  seed -- is restored on exit, including on error. Verified: `0.76430`
  (rounds to 0.764) and bias `-0.08570` reproduce identically across
  Mersenne-Twister, Wichmann-Hill, Marsaglia-Multicarry, Super-Duper,
  Knuth-TAOCP-2002 and L'Ecuyer-CMRG; `ss_time_dependent_roc()`'s
  `n_total = 420` reproduces identically across Inversion, Box-Muller,
  Kinderman-Ramage and Ahrens-Dieter.

* **`nb_assurance_ceiling()`'s precision claim rested on an observation,
  not a bound.** The previous default, `B_ceiling = 30000000`, was
  justified by the largest error OBSERVED in a six-scenario cross-check
  against an independently derived (Gauss-Legendre quadrature) reference --
  a sample statistic from six data points, not a guarantee that covers the
  next call. Measured directly for the package's own default `pt_range =
  c(0.15, 0.40)`: the true ceiling there (Gauss-Legendre reference,
  0.879606) sits only `1.06e-4` above `0.8795`, the boundary at which the
  value published to three decimals flips from 0.879 to 0.880, while the
  Monte Carlo standard error of the estimator at the previous default is
  close to `5.6e-05` there -- under two standard errors from that
  boundary, putting an estimated 3\% of seeds on the wrong side of the
  published rounding. `B_ceiling` is now 200000000 (2e8), which keeps the
  standard error (a real bound on the estimator, `sqrt(p(1-p)/B_ceiling)`,
  not an observation) below `2.6e-05` across the typical 0.7-0.95 ceiling
  range, pushing the same boundary more than four standard errors away.
  This does not move the published `nb_ceiling` table (0.880, 0.969, 0.956,
  0.695, 0.753, 0.000 to three decimals): the new default's value for the
  package's own `pt_range` is 0.8795889, 1.71e-05 from the Gauss-Legendre
  reference and still comfortably on the 0.880 side of the rounding
  boundary. The added precision costs roughly two minutes per call on
  ordinary hardware (draws are fully vectorized, so this does not scale
  with `B` or `N_range`), against a few seconds at the previous default;
  a caller who only needs an approximate ceiling can still pass a smaller
  `nb_B_ceiling`/`B_ceiling` for speed.

## Bug fixes (published numbers unaffected)

* **`ss_unified()` did not validate `N_range` at all.** A negative, `NA`,
  infinite, or non-integer element was not rejected up front; it instead
  reached `rbinom()`/`rmultinom()` deep inside the search loop and died
  with a generic, uninformative `"missing value where TRUE/FALSE needed"` --
  a crash identifying neither the offending argument nor which of its
  elements was at fault. `N_range` is now validated on entry, with a
  message naming the specific problem (and, where useful, the offending
  position or value) for each of: not numeric / zero-length, containing
  `NA`, containing a non-finite value, containing a value below 1, and
  containing a non-integer value. Separately, an `N_range` extended with a
  degenerate value such as `0` (previously accepted, contributing to the
  same block-pooling mechanism as the duplicated-N defect above) is now
  rejected outright rather than silently shifting `N_effective`.

## Documentation

* `man/ss_unified.Rd` and `R/ss_unified.R` claimed that under
  `decision = "isotonic"` both `assurance_mcse` AND `assurance_lower` are
  reported as `NA`. Only `assurance_mcse` is; `assurance_lower` is the
  margin-adjusted (block-pooled Wilson) curve actually inverted to select
  `N_effective`, and is always a real number when `N_effective` is found.
  Only the documentation changes; the returned values were already
  correct.

# dtasamplesize 0.6.1

**Validity release.** Five further defects were found: four downstream of how
0.6.0 seeded its Monte Carlo draws or bounded a small-sample estimate, and one
in how `decision = "isotonic"` (new in 0.6.0) paired `stats::isoreg()`'s
fitted curve back up with the candidate N it belonged to. None of this
release's fixes move the manuscript's headline sample sizes (the
Buderer/BAM/joint/imperfect-reference/unified comparisons, Figure 4's nested
sequence, or `ss_net_benefit()`'s N = 240) computed under the package's own
default RNG configuration and realistic `B`; the `nb_ceiling` table's values
do change, becoming more accurate (see below).

## Validity fixes (results change for affected callers)

* **`ss_unified()`'s `decision = "isotonic"` silently paired each N with the
  wrong fitted assurance whenever `N_range` was not already sorted
  ascending.** `stats::isoreg()` returns `$x` in the ORIGINAL order of its
  input but `$yf` in the SORTED order -- an internal detail of its
  implementation, not part of its documented contract. The fit was built
  directly from `grid_results$N`/`grid_results$assurance`, which follow
  `N_range`'s own order, and `$x`/`$yf` were then read off positionally, so
  every N was silently paired with a fitted value belonging to a
  *different* N whenever `N_range` was given in any order other than
  already-increasing (descending, shuffled, a hand-typed `c()` not given in
  order, ...). On one reproducible scenario (`prior_prev = c(4, 16)`,
  `delta_auc = 0`, `check_nb = FALSE`, `target_assurance = 0.80`, identical
  `B`, `seed` and *set* of candidate N throughout) this selected
  `N_effective = 729` (`joint_assurance = 0.8104`) given `N_range` ascending,
  `677` given descending, and, given a shuffled order, `415` with
  `joint_assurance = 0.4156` and `assurance_lower = 0.4041` -- **below**
  `target_assurance = 0.80` -- with no warning of any kind. The fit is now
  built from an explicitly pre-sorted grid (`order(grid_results$N)` applied
  before calling `isoreg()`), which removes the ambiguity entirely: once the
  input is already increasing, `isoreg()`'s "original order" and "sorted
  order" coincide by construction, so `$x` and `$yf` come back aligned no
  matter what order `N_range` itself was in. New defensive checks abort
  outright (rather than warn) if this alignment is ever violated, since the
  original failure mode was itself a warning-free silent mismatch. This is
  the most serious defect of this release: unlike the four below, it could
  select an N whose true assurance falls *under* the target with no
  diagnostic at all. The manuscript's own `N_range` arguments are always
  given already sorted ascending, so none of its published numbers were
  affected.

* **`ss_unified()`'s internal `nb_assurance_ceiling()` helper inherited the
  caller's active RNG generator.** It saved and restored `.Random.seed` but
  called `set.seed(seed)` with no `kind`, so the draws behind the
  `check_nb` ceiling -- and therefore whether the grid search was skipped
  as unreachable -- silently depended on whatever generator the calling
  session had active, the same class of defect 0.6.0 fixed in
  `ss_unified()` itself but left uncorrected in this helper. With
  identical arguments and seed, the reported ceiling ranged from 0.28253
  (Marsaglia-Multicarry) to 0.28539 (Knuth-TAOCP-2002) across generators;
  near a ceiling boundary, some generators let the search run while others
  correctly skipped it as unreachable. `set.seed()` now fixes the kind
  explicitly to `"Mersenne-Twister"` (R's own default) and both the kind
  and the seed are restored on exit, including on error.

* **Six other exported functions had the same generator-inheritance
  defect: `bam_sample_size()`, `joint_sample_size()`,
  `mc_validate_buderer()`, `ss_adaptive_prevalence()`, `ss_net_benefit()`,
  and `ss_time_dependent_roc()`.** Each saved and restored
  `.Random.seed`'s value but called `set.seed(seed)` with no `kind`, so
  their results depended on the caller's active generator rather than on
  `seed` alone -- for a package whose stated purpose is reproducibility,
  the published numbers were not actually reproducible by a reader whose
  session had a different default generator active. Measured with the
  manuscript's own arguments and seed, `ss_net_benefit()`'s conservative N
  ranged from 230 (Wichmann-Hill, Knuth-TAOCP-2002) to 240
  (Mersenne-Twister, L'Ecuyer-CMRG), and `mc_validate_buderer()`'s
  `P_width_target` (the 57% figure in Figure 1) ranged from 0.56625 to
  0.57575. Every `set.seed()` call in these six functions now fixes
  `kind = "Mersenne-Twister"` explicitly, and each function restores both
  the caller's RNG kind and `.Random.seed` on exit (not just the seed
  value, as before), via the same mechanism `ss_unified()` already used
  (new internal `save_rng_state()` / `restore_rng_state()`). Because
  Mersenne-Twister is R's own default generator, **every number in the
  manuscript and in this package's examples is unchanged** by this fix.
  Naming `kind` fixes the *uniform* generator only; five of these six
  functions draw only from it (`stats::rbeta()`/`rbinom()`/`rmultinom()`),
  so their results are now reproducible regardless of the reader's own RNG
  configuration -- but `ss_time_dependent_roc()` also draws from the
  *normal* generator (`stats::rnorm()`, for the simulated biomarker), which
  `normal.kind` controls independently of `kind`, and 0.6.1 left that
  unfixed; see 0.6.2 for the correction. `ss_unified()` itself remains
  unaffected here, as it already fixed its generator (L'Ecuyer-CMRG,
  needed for its independent per-N streams) explicitly.

* **`ss_unified()` could report a vacuous `assurance_lower = 1.000` (and
  `joint_assurance = 1.000`) at a very small `B`.** All three decision
  rules computed the lower confidence bound on the Monte Carlo assurance
  with the ordinary normal approximation (Wald), whose standard error is
  *exactly* 0 whenever the estimated proportion is 0 or 1, regardless of
  `B`: at `B = 1`, a single replicate that happened to pass every active
  target gave `joint_assurance = 1` with a Wald standard error of
  `sqrt(1 * 0 / 1) = 0`, so `assurance_lower` came out as `1.000` --
  manufactured certainty from one replicate -- and every decision rule
  (`"point"`, `"lower_bound"`, `"isotonic"`) accepted the first N tried on
  that basis. `assurance_lower` (and, under `decision = "isotonic"`, the
  margin subtracted from the fitted curve before inversion) is now a
  one-sided Wilson score bound instead, which stays strictly inside
  `(0, 1)` for any finite `B` -- at `B = 1` and `joint_assurance = 1` it
  gives approximately 0.27, not 1.000. The Wilson and Wald bounds agree to
  `O(z^2 / B)`, on the order of `1e-4` or smaller already at the `B =
  20000` used for the manuscript's own results, so **no result computed
  at a realistic `B` is affected**; only the small-`B` behaviour is.
  `assurance_mcse` is unchanged (still the descriptive Wald standard
  error, documented as such).

* **`nb_ceiling`'s six published values (Table 5) were systematically low
  by about 0.001, and their errors were perfectly correlated with each
  other.** `nb_assurance_ceiling()` used a single shared sample of 200,000
  prior draws, generated once from `seed` alone; a caller evaluating
  several `pt_range` values under the same `seed` (as the manuscript's
  Table 5 does, one row per range) therefore drew the *exact same*
  `(prev, Se, Sp)` triplet for every row, since `pt_range` only entered
  the calculation after the draws. Whichever direction that one sample's
  noise happened to point, every row inherited it in the same direction
  and to a similar relative size -- a table of nominally independent
  estimates whose Monte Carlo errors did not average out across rows. Two
  changes address this: (1) `seed` is now combined with a deterministic,
  order-independent fingerprint of `pt_range` before seeding, so that
  different `pt_range` values draw statistically independent samples
  under the same nominal `seed`, while a given `(seed, pt_range)` pair
  stays exactly reproducible; and (2) the default number of prior draws
  rises from 200,000 to 30,000,000, which keeps the Monte Carlo standard
  error under about `1e-4` for ceilings in the typical 0.7-0.95 range
  (verified against an independently derived reference computed by
  tensor-product Gauss-Legendre quadrature: the largest observed absolute
  error across six `pt_range` scenarios was `6.1e-5`). The manuscript's
  Table 5 values move accordingly, e.g. `pt_range = c(0.15, 0.40)` from
  0.87864 to approximately 0.8796 (independently verified value:
  0.879606); `validation/make_manuscript_assets.R`'s anchor check for
  this table now compares against that independently derived value
  instead of a number read off the package's own prior output, which was
  circular (it confirmed only that the code reproduces itself, not that
  it is correct).

# dtasamplesize 0.6.0

**Validity release, breaking.** Five further defects were found, three of
which change published numbers. This release does not touch the numbers
`bam_sample_size()` and `joint_sample_size()` returned in the manuscript
(678, 672, 580), but it does change what `ss_adaptive_prevalence()`,
`ss_imperfect_ref()`, and `ss_unified()` return, and it changes them in
different directions -- some numbers go up, `ss_unified()`'s search also
gets tighter around the target it claims to reach. We state this plainly
because users may have planned studies with the old numbers: if you used
`ss_adaptive_prevalence()`, `ss_imperfect_ref()`, or `ss_unified()` under
0.5.0 or earlier, recompute your sample size.

## Validity fixes (results change)

* **`ss_adaptive_prevalence()` discarded the stage-1 pilot and re-sampled
  the full stage-2 sample from scratch.** Stage 1 was simulated only to
  produce the interim prevalence estimate and then thrown away; stage 2
  then drew the *entire* re-estimated `N_final` again, so the stage-1
  subjects were recruited but never analysed and never counted. This is
  not what an internal-pilot design does (Stark & Zapf 2020): the whole
  point is that the pilot's subjects are folded into the final analysed
  sample, and stage 2 tops up rather than starting over. Stage 1 is now a
  genuine internal pilot -- its subjects, and the diseased/non-diseased
  counts observed in them, are retained and analysed alongside stage 2.
  The recruitment the old code reported **understated the true N needed
  between 21% and 36%**; the reported N now matches actual recruitment.

* **`ss_imperfect_ref()` multiplied the sample size for the index test's
  sensitivity and specificity by the Rogan-Gladen variance-inflation
  factor, `1 / (Se_ref + Sp_ref - 1)^2`.** That factor is a correct
  identity for the variance of a *prevalence* estimated from an imperfect
  screening test; it is not the right multiplier for the variance of an
  index test's Se/Sp evaluated against an imperfect reference, which the
  package now derives from the exact, identified 2x2-table model instead
  of borrowing a formula from a different problem. A single inflation
  factor also cannot represent this correctly, because an imperfect
  reference forces a choice of *what* is being sized: the study can be
  powered for the *apparent* Se/Sp (what a naive analysis against the
  imperfect reference actually estimates, `P(T+|R+)`/`P(T-|R-)`) or for
  the misclassification-*corrected* Se/Sp (the true index-test accuracy,
  recovered from the full 2x2 table). New argument
  `estimand = c("apparent", "corrected")`, defaulting to `"apparent"`;
  **both sample sizes are always computed and returned**
  (`N_apparent`/`N_apparent_loss` and `N_corrected`/`N_corrected_loss`), so
  a caller who reads only the generic `N_adjusted`/`n_total` cannot
  silently under-power a corrected-estimand analysis by relying on the
  default. In the manuscript's configuration, `N_adjusted` moves from 695
  to **732** (apparent estimand, the new default) or to **1235**
  (corrected estimand); with 10% expected loss to follow-up, to **814**
  and **1373** respectively. The old `VIF` field is kept, unchanged, as an
  informative diagnostic of the reference standard's information about
  *prevalence* -- it no longer multiplies any sample size.

* **`ss_unified()` accepted the first N in the search grid whose lower
  Monte Carlo confidence bound reached the target.** Because the grid is
  evaluated at successive N with the generator re-seeded from the same
  `seed` at every N, the per-N replications were not independent draws --
  they agreed on the first replication and diverged unpredictably from
  the second, worse than either genuine common random numbers or genuine
  independence. Repeating the search across seeds, the selected N had
  standard deviation **96.8**, and only **84%** of the selected N's
  actually reached `target_assurance` when checked against an independent
  high-precision reference -- an 18-point undershoot of the nominal 95%
  guarantee the old `decision = "lower_bound"` rule was supposed to
  provide.

  `ss_unified()` now seeds an L'Ecuyer-CMRG generator once and advances to
  a fresh, independent stream per candidate N (via
  `parallel::nextRNGStream()`), and gains `decision = "isotonic"`
  (**new default**): it evaluates the whole `N_range`, fits a monotone
  non-decreasing curve to the assurance-by-N points with
  `stats::isoreg()`, and inverts that fitted curve at `target_assurance`
  after subtracting a conservative Monte Carlo margin, rather than reading
  off the first noisy crossing. Measured the same way across 17 seeds,
  the selected N now has standard deviation **33.5**, and **100%** of the
  selected N's reach the target. `decision = "lower_bound"` (the previous
  default) and `decision = "point"` (the 0.4.0 behaviour) remain available
  unchanged.

* **`ss_unified()`'s `check_nb = TRUE` search had no way to detect an
  unreachable target and would simply exhaust `N_range`, suggesting the
  caller expand the grid** -- advice that cannot help when the criterion
  cannot be met at *any* N. The function now computes, before the grid
  search, the criterion's achievable **ceiling** as N -> infinity (a fast
  vectorized calculation over the priors, `Se_ref`, `Sp_ref` and
  `pt_range` alone, independent of `B` and `N_range`) and returns it as
  `nb_ceiling`. If `target_assurance` exceeds this ceiling, the grid
  search is skipped and the function warns that the target is unreachable
  under these priors, rather than recommending a wider grid. With
  `ss_unified()`'s default priors and `pt_range = c(0.15, 0.40)` the
  ceiling is **0.880**; narrowing the top of the range to
  `pt_range = c(0.10, 0.40)` lowers it to **0.695**, and widening it to
  `pt_range = c(0.05, 0.50)` lowers it to **0.000** -- no N reaches the
  target under that range, because a wide `pt_range` squeezes the
  criterion's two inequalities from both ends at once.

## Bug fixes (published numbers unaffected)

* **`bam_sample_size(method = "exact")` did not apply the
  degenerate-replicate rule its own documentation declared.** A
  replicate with no diseased or no non-diseased subjects is supposed to
  count as a failure regardless of how narrow the prior-only credible
  interval happens to be -- the same convention `method = "monte_carlo"`
  already enforced. The exact calculation instead evaluated the
  credible-interval width at the *prior alone* for a zero-size arm, which
  is not 0 in general and, for a sufficiently informative prior, can be 1.
  With an informative prior this let a degenerate N = 10 replicate count
  as a success purely because the prior interval already met the target:
  the exact joint guarantee at N = 10 came out as **1.000**, against a
  true value of **0.344**. Degenerate replicates now always score as
  failures, matching `method = "monte_carlo"` exactly. **The manuscript's
  published numbers do not change** (678 and 672): in that scenario the
  degenerate-replicate contribution was negligible.

* **`ss_net_benefit()` returned `NA` with no explanation when the
  criterion could not be met at any sample size.** The joint criterion
  requires, at every threshold, that the test beat both the treat-none
  and the treat-all strategy. Whether it can do so at all is a structural
  property of `Se`, `Sp`, `prev` and the threshold, not of the sample
  size: if the limiting net benefit is not positive, no `N` helps. The
  function now computes those limiting values in closed form, names the
  thresholds that cannot be reached and which of the two comparisons
  fails at each, and says plainly that expanding `N_range` will not help.
  The per-threshold verdicts are exposed in `N_by_pt` so they can be
  inspected rather than only read from a message. A search that fails
  because the grid was too small still gets the original advice to widen
  it. Published numbers are unchanged (the conservative `N` = 240 and the
  per-threshold sizes).

* **`ss_unified()` now restores the random-number generator *kind*, not
  only `.Random.seed`.** Drawing an independent stream per candidate `N`
  requires the L'Ecuyer-CMRG generator. Restoring the seed vector alone
  would leave that generator selected for the rest of the session, so any
  function called afterwards would return different results even with an
  explicit seed of its own. Both the kind and the seed are now restored
  on exit, including when the function exits with an error.

## Documentation

* `ss_time_dependent_roc()`'s `censoring_rates` was documented as
  "censoring proportions." It is not: each entry is the marginal
  probability that the *censoring time* precedes `t_horizon`, calibrating
  an exponential censoring model. Because censoring competes with the
  event, the censoring actually observed in the simulated data is lower
  than this nominal value -- roughly 80% of it, for the package's default
  `lambda_event` and `t_horizon` (the ratio is not universal and shifts
  with those parameters). The documentation now states this explicitly
  and the function's warning message says "nominal censoring rate"
  instead of "censoring rate." **Only the documentation changes; the
  simulation and every numeric result are identical.**

## Dependencies

* `DESCRIPTION` now lists `parallel` (an R base package) under `Imports`,
  used by `ss_unified()` to generate independent per-N random streams via
  `parallel::nextRNGStream()`.

# dtasamplesize 0.5.0

**Validity release.** `bam_sample_size()` searched sensitivity and specificity
requirements independently and never checked that both credible intervals
reached their target in the same study; `ss_unified()` accepted the first N at
which a Monte Carlo search happened to cross the target, without accounting
for the sampling error of that crossing; and `joint_sample_size()` fixed the
number of diseased subjects at its expected value instead of letting it vary
as it does in a real cohort. These are not cosmetic: the sample sizes
`bam_sample_size()` and `ss_unified()` returned in 0.4.0 do not reach the
assurance they declared, and a study planned with those values is
under-sized. We state this plainly because users may have planned studies
with the old numbers: if you used `bam_sample_size()` or `ss_unified()` under
0.4.0 or earlier, recompute your sample size.

## Validity fixes (results change)

* **`bam_sample_size()` reported a marginal per-arm guarantee, not a joint
  one.** It searched `n_se` and `n_sp` independently, each against its own
  marginal credible-interval target, and returned the median of
  `max(n_se / prev, n_sp / (1 - prev))`. It never checked that both credible
  intervals met their target in the *same* study. For the manuscript
  scenario, the true joint assurance of the N it returned (591) was
  **0.7245**, not 0.80; the first N whose joint assurance reaches 0.80 is
  **678** -- an understatement of 87 participants (14.7%).

  It now searches directly for the smallest total N whose **joint**
  assurance reaches the target, under a cohort design (the number of
  diseased subjects is random, `Binomial(N, prev)`), and degenerate
  replications count as failures with denominator `B`. It also gains
  `method = c("exact", "monte_carlo")`, **defaulting to `"exact"`**: the
  joint assurance has closed form (a Beta-Binomial sum), so it is computed
  deterministically, with no Monte Carlo error and no dependence on `B` or
  the seed. The legacy fields (`n_diseased`, `n_non_diseased`,
  `N_total_median`, ...) are kept and documented as per-arm diagnostics.

* **`joint_sample_size()` fixed the number of diseased subjects at
  `floor(N * prev)`.** Conditioning on that expected value overstates the
  assurance in a prospective cohort, where the count is random -- the same
  defect fixed in `ss_net_benefit()` in 0.3.0. For the manuscript scenario
  (N = 580) the reported assurance was 0.856 and the true cohort assurance
  is **0.808**: the N remains adequate, but the real margin is far smaller
  than declared. In an adversarial case (Se = 0.70, Sp = 0.80, prev = 0.20)
  the difference is 0.807 against 0.716, enough to overturn the conclusion.

  New argument `design = c("cohort", "fixed")`, **defaulting to
  `"cohort"`**. `"fixed"` reproduces the previous behaviour and is
  documented as valid only when the two groups are recruited separately
  with pre-specified sizes.

* **`ss_unified()` accepted the first crossing of a noisy grid.** The
  algorithm and the success criterion were correct, but at `B = 1200` the
  Monte Carlo standard error of the assurance near 0.80 is **≈0.0115**.
  Repeating the same search with ten seeds, the selected N ranged from
  **900** to **940**. Evaluated at `B = 300000`, the N that had been
  published reached **0.793**, not 0.80.

  An N is now accepted only when the **lower bound** of the one-sided 95%
  interval of the assurance reaches the target (`decision = "lower_bound"`,
  the default; `"point"` recovers the previous behaviour). The result now
  also carries `assurance_mcse` and `assurance_lower`, and the function
  warns when the Monte Carlo error is large relative to the distance to the
  target.

## Reproducibility

* New `validation/make_manuscript_assets.R`, the public generator of every
  figure and of the manuscript's derived tables. This script was not
  previously part of the repository.
* Table 1 is now derived from the `NAMESPACE` instead of being maintained by
  hand.
* The capability matrix comparing this package against alternatives now
  lives in a version-controlled CSV, accompanied by a provenance log
  recording the exact versions of `presize`, `pROC`, `MKpower`, `epiR`, and
  `MKmisc` examined, and the date they were examined.
* `validation/reproduce_manuscript.R` now covers 20 checks, including the
  five steps of the nested comparison.

# dtasamplesize 0.4.0

**Usability release.** No published number changes: every quantity reported
in the manuscript and reproduced by `validation/reproduce_manuscript.R` is
identical to 0.3.0. This release adds a sensitivity-analysis helper, plotting
functions, several fields that were previously computable only by re-deriving
them by hand, backward-compatible naming aliases, and an opt-out advisory
warning for under-powered Monte Carlo runs.

## New features

* **`sensitivity_analysis()`** runs `ss_unified()` once for each of a set of
  prior scenarios (e.g. optimistic / neutral / pessimistic assumptions about
  Se, Sp, and prevalence) and collects the resulting sample sizes into a
  single data frame (`scenario`, `N_total`, `joint_assurance`, `N_buderer`,
  `n_diseased`). Scenarios can be supplied as a named list or as a data
  frame with one row per scenario. Replaces the hand-written `lapply()` loop
  this comparison otherwise requires.
* **Three plotting functions**, all requiring `ggplot2` (in `Suggests`,
  never `Imports`; each fails with an informative message if it is not
  installed):
  * `plot_assurance_curve()` plots achieved joint assurance against total
    sample size from an `ss_unified()` result's `grid_results`, with the
    target assurance and the optimal N marked.
  * `plot_inflation_heatmap()` plots the variance inflation factor from an
    `ss_imperfect_ref()` result's `sensitivity_table` as a heatmap over
    reference-standard sensitivity and specificity.
  * `plot_method_comparison()` plots the required total N by method from an
    `ss_unified()` result's `comparison` table.
* **`ss_unified()` gains `full_grid`** (default `FALSE`, preserving current
  behaviour and cost). The search over `N_range` normally stops at the first
  N that reaches `target_assurance`; with `full_grid = TRUE` it evaluates
  every N in `N_range` instead, so that the new `grid_results` element (a
  data frame of every `(N, assurance)` pair evaluated) holds the complete
  assurance curve rather than one truncated at the optimum.
* **`ss_unified()` also returns `N_buderer`, `N_imperfect`, `seed`, and
  `target_assurance`**, so the classical and imperfect-reference comparison
  values, the search seed, and the assurance target no longer have to be
  re-read out of the `comparison` table or the call.
* **A small-`B` advisory warning.** Every exported function that runs a
  Monte Carlo search (`mc_validate_buderer()`, `bam_sample_size()`,
  `joint_sample_size()`, `ss_net_benefit()`, `ss_adaptive_prevalence()`,
  `ss_imperfect_ref()` when `B > 0`, `ss_time_dependent_roc()`, and
  `ss_unified()`) now warns when `0 < B < 1000`, since the Monte Carlo error
  of the reported assurance can be substantial at low `B`. The warning is
  controlled by the `dtasamplesize.warn_small_B` option (default `TRUE`);
  set `options(dtasamplesize.warn_small_B = FALSE)` to silence it, e.g. for
  scripts that deliberately use a small `B` for speed.

## Backward-compatible aliases

* `ss_imperfect_ref()` accepts `delta_se` / `delta_sp` as aliases of
  `d_se` / `d_sp`, matching the `delta_*` naming used elsewhere in the
  package.
* `ss_net_benefit()` accepts `target_assurance` as an alias of
  `target_prob`, matching the `target_assurance` naming used in
  `ss_unified()` and `bam_sample_size()`.
* All aliases are non-breaking: the original argument names keep working
  unchanged, and when both an old and a new name are supplied the new name
  takes precedence.

## Enhancements

* `print.dtasamplesize()` now also shows `joint_assurance`, `N_buderer`,
  `N_imperfect`, `seed`, and `B` when the object carries them.
* `ss_imperfect_ref()`'s `sensitivity_table` gains an `inflation_factor`
  column, an identically-valued but more descriptive alias of `VIF`; both
  are kept so existing code reading `VIF` still works.
* `ss_unified()` and `ss_time_dependent_roc()` gain `@references` (O'Hagan
  et al. 2005, Wilson et al. 2022, Rogan & Gladen 1978, and Hanley & McNeil
  1982 for the former; Heagerty et al. 2000 and Blanche et al. 2013 for the
  latter), bringing every exported sample-size function up to the same
  documentation standard.

## Documentation

* New vignette section, "Sensitivity analysis across prior scenarios",
  demonstrating `sensitivity_analysis()` over an optimistic / neutral /
  pessimistic set of priors and showing how to read the resulting table.
* `@param B` across all Monte Carlo functions now documents the small-`B`
  warning and how to silence it.

## Tests

* New tests for `grid_results` and `full_grid` (row counts, and that
  `N_effective` / `n_total` / `joint_assurance` are unaffected by
  `full_grid`), for the new `ss_unified()` fields and their coherence with
  `comparison`, for the `delta_se` / `delta_sp` and `target_assurance`
  aliases, for `sensitivity_analysis()` (both scenario formats, column
  names and row order, and its error messages), for the small-`B` warning
  (on and off), and for the three plotting functions.
* `tests/testthat/setup.R` sets `dtasamplesize.warn_small_B = FALSE` for the
  suite, since many tests deliberately use a small `B` for speed; the
  warning itself is tested explicitly with the option restored to `TRUE`.

# dtasamplesize 0.3.0

**Statistical-validity release.** An independent review of version 0.2.0 found a
**validity error** in `ss_net_benefit()` and `ss_adaptive_prevalence()`, plus a
geometrically impossible default in `joint_sample_size()`. These were not
cosmetic: they made the package **understate the required sample size** and, in
one case, report that *more* dropout buys *more* precision. Results produced with
0.2.0 using these functions should be recomputed. We state this plainly because
users may have planned studies with the old numbers.

## Validity fixes (results change)

* **`ss_net_benefit()` applied a fixed-margin variance to prospective cohorts.**
  Versions <= 0.2.0 fixed the number of diseased subjects at `floor(N * prev)`
  and used the variance *conditional* on that, which is correct only when the two
  disease groups are recruited separately with pre-specified sizes. In a
  prospective cohort disease status is **random**, and conditioning it away
  discards a large variance component: at
  `Se = 0.85, Sp = 0.90, prev = 0.20, pt = 0.20, N = 500` the standard error used
  was **0.0077** against a true cohort sampling SD of **0.0174** -- 2.25x too
  small. The declared assurance was therefore well above the real one, and the
  recommended `N` too small.

  New argument `design = c("cohort", "fixed")`, **default `"cohort"`**. Under the
  cohort design the multinomial variance is used, and because the treat-all net
  benefit is then **also estimated**, the treat-all comparison is made on the
  confidence interval of the *difference* `NB - NB_all`, which simplifies exactly
  to `(w * TN - FN) / N`. `design = "fixed"` reproduces the old behaviour and is
  documented as **not valid for a prospective cohort**.

  Required `N` roughly doubles at the extreme thresholds; the conservative `N`
  over the default threshold range goes from **110** to **240**. The `N_range`
  default was widened to `seq(50, 2500, by = 10)`.

* **`ss_adaptive_prevalence()` inflated for losses but never applied them.** The
  function computed `N_final_adj = N_final / (1 - loss_rate)` and then simulated
  and analysed **all** `N_final_adj` subjects, so the inflation became free extra
  data. The consequence was absurd: `precision_achieved` **rose** with the dropout
  rate (0.70 at 0% losses, 0.83 at 10%, 0.98 at 30%). Losses are now actually
  simulated and only the surviving sample is analysed; precision is flat in
  `loss_rate`, as it must be. New column `N_analysed_median` reports the post-loss
  sample size; the `N_final_*` columns remain the **recruited** sizes.

* **`joint_sample_size()` shipped a geometrically impossible default.**
  `AUC = 0.80` with `Se = 0.85` and `Sp = 0.90` describes a test that cannot
  exist: any concave ROC through the operating point `(1 - Sp, Se)` has area at
  least `0.5*(1-Sp)*Se + 0.5*Sp*(Se+1)` = **0.875**. The package computed a sample
  size for it regardless. Such combinations now raise an error naming the minimum
  attainable AUC, and the **default is now `AUC = 0.90`** (above the geometric
  minimum, below the binormal value of 0.949 implied by the same operating point,
  hence conservative).

* **`ss_unified()` used the net-benefit criterion it documents as abandoned.**
  `check_nb = TRUE` required only that the *point estimate* of NB beat the default
  strategies -- a condition met with near-certainty whenever the test is useful, so
  it barely moved `N`. It now uses the same CI-based criterion as
  `ss_net_benefit(design = "cohort")`. With `prior_sp = c(18, 2)` the required `N`
  moves from ~700 to **~1580**.

* **`ss_unified()` reported a conditional assurance.** Degenerate replications
  (fewer than 5 diseased or non-diseased) were dropped from the denominator,
  yielding an assurance *conditional* on the study not being degenerate. They now
  count as failures; the denominator is `B`.

## Breaking changes

* `ss_net_benefit()` gains `design`, defaulting to `"cohort"`. To reproduce 0.2.0
  exactly, pass `design = "fixed"`.
* `joint_sample_size()`: default `AUC` changed from 0.80 to 0.90; incompatible
  `Se`/`Sp`/`AUC` combinations now **error** instead of returning a number.
* `joint_sample_size()`: `joint_prob` is **renamed `joint_prob_se_sp`**, because it
  never was a three-way joint probability -- it is the Monte Carlo probability over
  Se and Sp only, while the AUC enters as a *deterministic* plug-in gate (a
  median-style criterion at roughly 50% assurance, not an assurance-style one).
  The documentation now says so explicitly. `joint_prob` is kept as a
  **deprecated alias**.
* `joint_sample_size()`: when the AUC gate is never passed, `joint_prob_se_sp` is
  now `NA` (the Se/Sp Monte Carlo was never run). Previously the field leaked its
  initialisation value of `0`, indistinguishable from a computed zero. New fields
  `AUC_min` and `auc_gate_passed`.
* `ss_imperfect_ref()`: the `method` string now says **"Rogan-Gladen variance
  inflation"** instead of "(Staquet Correction)", because that is what the code
  implements. New field `youden_ref`.
* `ss_adaptive_prevalence()`: `results` gains an `N_analysed_median` column.

## Bug fixes

* `ss_adaptive_prevalence()`: the interim prevalence estimate was bounded below
  (0.05) but **not above**, so an all-diseased stage-1 sample gave `prev_hat = 1`,
  hence `n_sp / (1 - prev_hat) = Inf` and a results table of `Inf`/`NA`. It is now
  truncated to `[0.05, 0.95]`; the truncation is documented, and when it binds
  the function now emits a warning so the caller knows the diseased arm may be
  under-sized if the true prevalence lies outside that range.
* `ss_imperfect_ref()`: `Se_ref = Sp_ref = 0.51` is a **legal** input (it passes
  `Se_ref + Sp_ref > 1`) but gives a variance inflation factor of 2500 and made the
  Monte Carlo validation attempt a ~9 Gb allocation. Two documented, overridable
  guards were added: `min_youden` (default 0.5, capping the inflation factor at 4)
  and `max_mc_cells` (default 2e7, bounding `B * N_adjusted`; this also catches an
  extreme `prev`, which inflates `N` without touching the Youden index). Both fail
  with a message naming the offending quantities.

## RNG hygiene (CRAN policy)

* All eight exported functions that call `set.seed()` internally
  (`mc_validate_buderer()`, `bam_sample_size()`, `joint_sample_size()`,
  `ss_net_benefit()`, `ss_adaptive_prevalence()`, `ss_imperfect_ref()`,
  `ss_time_dependent_roc()`, `ss_unified()`) **overwrote the caller's global
  `.Random.seed` and never restored it**, silently changing the random numbers a
  user obtained after calling them. They now save and restore the RNG state on
  exit (and leave no `.Random.seed` behind in a session that never had one).

## Documentation

* `ss_unified()` gains a `@note` on the **apparent-sensitivity bias**: it estimates
  `P(T+ | R+)`, not `P(T+ | D+)`. Sizing by this framework buys *precision* about
  the apparent sensitivity; it does not remove the *bias*.
* `ss_net_benefit()` documents both designs, their variances, and why the choice
  matters.
* `joint_sample_size()` documents the geometric constraint linking Se, Sp and AUC.

## Tests

* Seven tests in 0.2.0 were green while exercising the **non-converged fallback
  branch** (search grids that never reached the real answer, e.g.
  `bam_sample_size(n_range = 20:200)` when the real `n_sp` is 373). Grids were
  widened past the true answers and `expect_no_warning()` added to pin convergence;
  the fallback branches are still tested, but deliberately, with
  `expect_warning()`.
* New regression tests for every fix above, including `test-rng_state.R` and
  independent re-simulations that verify the declared cohort assurance is the one
  actually attained.
* Suite: **0 failures, 0 warnings, 192 passing, 1 skipped** (0.2.0: 0 failures, 7 warnings,
  96 passing).

## Vignette and Shiny app

* **The introduction vignette had the same two defects.** Its BAM chunk used
  `n_range = 20:200` while the real `n_sp` is 370, so the vignette ran on the
  non-converged fallback and **published an understated total N (388 instead of
  563)**. And its comparison table put Buderer, BAM and the imperfect-reference
  correction at prevalence 0.30 but the joint calculation at 0.20 -- the same
  method-by-prevalence confound as Figure 4. Both are fixed: the grid is widened,
  the joint call takes `prev = 0.30`, and the text now declares the common
  prevalence. The vignette's comparison table moves from
  (334, 388, 600, 516) to **(334, 563, 400, 516)**.
* **The Shiny app's own default was the impossible combination.** Its AUC input
  defaulted to 0.80 while Se/Sp defaulted to 0.85/0.90, so after the geometric-constraint
  validation the Joint tab would have errored on its default settings. The
  default is now 0.90, with help text stating the geometric constraint. The app
  also had no error handling whatsoever, so any input the package now rightly
  rejects would have crashed it; errors and warnings are now surfaced as
  notifications.

## Manuscript assets

* **Figure 4 was a confounded comparison.** Its five bars were computed at
  **three different prevalences** (Buderer 0.20, Joint 0.20, Unified 0.25,
  BAM 0.30, imperfect-reference 0.30), so part of the difference between methods
  was really the prevalence -- and the caption declared none. All five are now
  computed at a **common prevalence of 0.20** (the Bayesian methods using a
  `Beta(4, 16)` prevalence prior, mean 0.20), and the caption states it. The
  values move from (500, 556, 600, 516, 800) to
  **(500, 591, 580, 773, 1000)** for Buderer / BAM / Joint / imperfect-reference /
  unified. The unified framework still demands the most; the imperfect-reference
  correction moves from 4th to 2nd once the confound is removed.

# dtasamplesize 0.2.0

Statistical-correctness release. Several functions were corrected for
statistical coherence and robustness. The most consequential change is to
`ss_net_benefit()`, whose results differ materially from 0.1.0.

## Breaking changes

* `ss_net_benefit()` now uses an **inference-based** sample-size criterion.
  A study is counted as successful only when the lower limit of the
  `(1 - alpha)` confidence interval for net benefit exceeds both 0
  (treat-none) and the treat-all net benefit. The previous version used the
  net-benefit *point estimate*, a criterion satisfied with near-certainty
  whenever the test was useful, so it returned the smallest `N` in the
  search range regardless of the parameters. The new criterion yields
  meaningful, threshold-dependent (U-shaped) sample sizes. A new `alpha`
  argument controls the confidence level, and the default `N_range` now
  starts at 50.
* `ss_net_benefit()` flags thresholds at which the test is not useful under
  the assumed parameters (net benefit not above 0 or not above treat-all)
  as `feasible = FALSE` with `N_required = NA`, instead of silently
  returning `max(N_range)`. The `N_by_pt` data frame gains a `feasible`
  column.

## Bug fixes

* `ss_net_benefit()`: fixed an error (`object 'prob' not found`) that
  occurred when every candidate `N` was skipped (e.g. very low prevalence).
* `bam_sample_size()`: the total sample size (`N_total_median`, `_P75`,
  `_P90`) now accounts for **both** the sensitivity (diseased) and
  specificity (non-diseased) requirements, i.e.
  `max(n_se / prev, n_sp / (1 - prev))`. Previously only the sensitivity
  arm was used, underestimating the total when the specificity arm was
  binding (with the vague default `prior_sp`, the underestimate exceeded
  30%).
* `ss_time_dependent_roc()`: the "best probability achieved" is now tracked
  explicitly (rather than relying on a leaked loop variable) and a warning
  is emitted when the target precision is not reached within `N_range`. The
  reported `(N_required, prob_achieved)` pair is consistent — the assurance
  at the largest `N` tried — with the grid-wide best conveyed in the warning.
* `ss_net_benefit()`: `N_conservative` (and the derived `n_total` /
  `n_diseased`) is now reported as `NA` with a warning when a *feasible*
  threshold could not reach the target assurance within `N_range`. Earlier
  it silently took the maximum over the thresholds that did converge,
  understating the worst-case requirement.
* `bam_sample_size()`: non-finite total-N draws (possible with very diffuse
  prevalence priors whose draws reach the 0/1 boundary) are dropped before
  computing the median and percentiles, preventing contaminated upper
  quantiles.

## Improvements

* `ss_imperfect_ref()`: the Monte Carlo validation now reports the mean
  apparent sensitivity `P(index+ | ref+)`, the true sensitivity, and the
  verification `bias` between them. The documentation clarifies that the
  variance inflation factor `1 / (Se_ref + Sp_ref - 1)^2` is the
  Rogan-Gladen misclassification-correction variance multiplier, used here
  as an approximation; it restores the *precision* of the apparent
  estimator but does not remove its *bias*. Reference to
  Rogan & Gladen (1978) added.
* `ss_unified()`: documentation default for `N_range` corrected to
  `seq(200, 1000, by = 20)` (matching the implementation); defensive
  initialisation added.

## Documentation / tests

* Added `NEWS.md`.
* Tests for `ss_net_benefit()` rewritten to check feasibility flagging, the
  non-trivial (U-shaped) sample-size pattern, and monotonicity in the
  assurance target. Tests for `ss_imperfect_ref()` now check that the
  adjusted sample size improves precision over the unadjusted one and that
  the apparent-sensitivity bias is exposed.

# dtasamplesize 0.1.0

* Initial release.
