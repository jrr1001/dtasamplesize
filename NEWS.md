# dtasamplesize 0.2.0

Statistical-audit release. Several functions were corrected for
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
