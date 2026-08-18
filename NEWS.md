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
