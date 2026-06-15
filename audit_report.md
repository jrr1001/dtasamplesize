# Audit Report — dtasamplesize

**Package:** dtasamplesize
**Audited version:** 0.1.0 → **0.2.0 (corrected)**
**Date:** 2026-06-15
**Scope:** Full statistical and software audit (correctness, coherence,
consistency, robustness) of every exported function, its documentation,
and its tests. All findings were reproduced empirically by running the
package under R 4.5.2; the validation script is provided in `validation/`.

> This report supersedes the previous `audit_report.md`, which concluded
> "0 issues found." That conclusion was incorrect: it
> reflected only that `R CMD check` and the test suite passed. Neither
> tool verifies statistical validity, and several of the issues below
> were latent precisely because the tests asserted the *buggy* behaviour.

---

## 1. Method of audit

`R CMD check --as-cran` and `testthat` confirm that the package *builds*,
*installs*, and *runs*. They do **not** confirm that the numbers returned
are statistically correct. The audit therefore consisted of:

1. Reading every function against its design specification and against the
   cited literature.
2. Re-deriving each statistical formula and checking the code matches.
3. **Executing** each function across normal, boundary, and degenerate
   inputs and inspecting the returned objects.
4. Prototyping each proposed correction and confirming the corrected
   behaviour numerically.
5. Re-running the full test suite and `R CMD check` after the fixes.
6. Cross-validating the core formulas against independent references
   (`validation/cross_validation.R`; see `validation_report.md`).

---

## 2. Findings and resolutions

Severity key — **CRITICAL**: wrong statistical result a user would act on;
**HIGH**: crash / silently misleading output; **MEDIUM**: coherence /
honesty of reported quantities; **LOW**: documentation / consistency.

### C1 (CRITICAL) — `bam_sample_size()` total N ignored the specificity arm

- **Symptom.** `N_total_median`, `N_total_P75`, `N_total_P90` were computed
  as `n_se / prev`, using only the sensitivity (diseased) requirement and
  ignoring the specificity (non-diseased) requirement `n_sp / (1 - prev)`.
- **Evidence.** With package defaults the specificity arm is binding:
  `n_se = 107`, `n_sp = 373`. The reported total was **368**; the value
  required to obtain 373 non-diseased at prevalence 0.30 is
  `373 / 0.70 ≈ 533`. Reported median after fix: **556**. The original
  understated the total sample size by more than 30%.
- **Why it matters.** A user planning a study from this output would
  enrol far too few participants to meet the specificity precision target.
- **Fix.** `N_total_draws <- pmax(n_se / prev_draws, n_sp / (1 - prev_draws))`,
  consistent with every other function in the package (all of which take
  the per-arm maximum). Documentation updated.
- **Status: FIXED.**

### H1 (HIGH) — `ss_net_benefit()` crashed on low prevalence

- **Symptom.** `Error: object 'prob' not found` whenever every candidate
  `N` was skipped (e.g. `prev = 0.002`), because the not-found branch read
  a loop variable that had never been assigned. This is the *same* class
  of bug the previous audit claimed to have fixed in `joint_sample_size()`
  ("H1"), but it persisted here and in `ss_time_dependent_roc()`.
- **Evidence.** `ss_net_benefit(prev = 0.002, ...)` → reproducible error.
- **Fix.** The not-found probability is now tracked in an explicitly
  initialised `best_prob`; infeasible/skip cases return gracefully (see
  H2/C2). The analogous fragile pattern in `ss_time_dependent_roc()` was
  hardened the same way and now emits a warning.
- **Status: FIXED.**

### C2 (CRITICAL) — `ss_net_benefit()` sample-size criterion was non-informative

- **Symptom.** "Success" was defined on the **point estimate** of net
  benefit: `(NB_hat > 0) & (NB_hat > NB_all_hat)`. Whenever the test was
  genuinely useful, the point estimate cleared the thresholds with
  probability ≈ 1 even at the smallest `N`, so the function returned the
  **floor of `N_range`** for every threshold.
- **Evidence.** Under defaults (Se 0.85, Sp 0.90, prev 0.20) every
  threshold returned `N_required = 100` (the old floor) with
  `prob_achieved ≈ 1.0`. The design spec expected ~180–370; that target
  was never attainable with the implemented criterion. The existing tests
  ("N increases with threshold", "N_conservative in [100,600]") passed
  *vacuously* (100 is in range; all-equal satisfies `diff ≥ 0`).
- **Fix.** Replaced with an **inference-based** criterion: a study succeeds
  when the lower limit of the `(1 − alpha)` confidence interval for net
  benefit exceeds both 0 (treat-none) and the treat-all net benefit, using
  the closed-form variance
  `Var(NB) = (n_d/N)² Se(1−Se)/n_d + (n_nd/N)² (pt/(1−pt))² Sp(1−Sp)/n_nd`.
  The required `N` is the smallest value at which this happens with
  probability ≥ `target_prob`. Added an `alpha` argument.
- **Resulting behaviour.** Sample size is now meaningful and
  **U-shaped in the threshold** (largest where the test is hardest to
  separate from a default strategy: low `pt`, where treat-all competes;
  high `pt`, where the net benefit itself is small) — the statistically
  correct shape, not the monotone increase the spec assumed.
- **Status: FIXED.** *(Behaviour change — see `NEWS.md`.)*

### H2 (HIGH) — `ss_net_benefit()` silently sized infeasible thresholds

- **Symptom.** When the test is **not** useful at a threshold under the
  assumed parameters (`NB_true ≤ 0` or `NB_true ≤ NB_treat_all`), no sample
  size can demonstrate superiority, yet the function returned
  `max(N_range)` as if it were a valid answer.
- **Evidence.** `Se = Sp = 0.55, prev = 0.10, pt = 0.05`:
  `NB_true = 0.034 < NB_treat_all = 0.053` (treat-all dominates), but the
  function returned `N_required = 300`.
- **Fix.** Feasibility is now tested first; infeasible thresholds are
  flagged `feasible = FALSE` with `N_required = NA`. A new `feasible`
  column is added to `N_by_pt`, and `N_conservative` is the max over
  *feasible* thresholds (`NA` if none).
- **Status: FIXED.**

### M1 (MEDIUM) — `ss_imperfect_ref()` Monte Carlo "validation" measured a biased quantity

- **Symptom.** The MC block reported the probability that the CI **width**
  of the *apparent* sensitivity `P(index+ | ref+)` met the target. The
  whole purpose of the imperfect-reference correction is that the apparent
  sensitivity is **biased**: against an imperfect reference it converges to
  a value well below the true Se. The MC therefore certified the precision
  of an estimator that is centred on the wrong value, with no disclosure.
- **Evidence.** True Se = 0.85; mean apparent Se = **0.764** (theoretical
  0.7639); empirical bias **−0.086** at both the unadjusted and adjusted
  sample sizes. The variance inflation makes the apparent-Se interval tight
  — but tight around 0.764, not 0.85. The existing test only checked
  `nrow == 2`, so it never detected the issue.
- **Fix.** (i) The `mc_validation` table now reports `se_apparent`,
  `se_true`, and `bias` for both scenarios, so the bias is visible.
  (ii) Documentation now states explicitly that the factor
  `1/(Se_ref + Sp_ref − 1)²` is the **Rogan–Gladen** misclassification
  variance multiplier used as an approximation; it restores *precision*
  but not *unbiasedness*, and recovering an unbiased estimate requires a
  bias-correction or latent-class analysis at the analysis stage. Added
  Rogan & Gladen (1978) to the references. (iii) Tests now assert that the
  adjusted sample size improves precision over the unadjusted one and that
  the bias is exposed and negative.
- **Status: FIXED (made honest).** The underlying VIF heuristic is retained
  but its scope and assumptions are now correctly stated — see the
  software report for the methodological discussion a reviewer will expect.

### L1 (LOW) — `ss_unified()` documented default did not match the code

- **Symptom.** `@param N_range` documented the default as
  `seq(200, 800, by = 20)`; the implementation used
  `seq(200, 1000, by = 20)`.
- **Fix.** Documentation corrected to match the code (`1000`). Defensive
  initialisation of `joint_assurance` added.
- **Status: FIXED.**

---

## 3. Items reviewed and judged correct (no change)

These were checked against the literature and/or numerically and are
sound; they are listed so the audit is transparent about what was *not*
changed.

- **`buderer_n()`** `= ceil(z² p(1−p)/d²)`: matches Buderer (1996);
  reproduces 100 (Se 0.85, d 0.07) and 139 (Se 0.90, d 0.05).
- **`wilson_ci()` / `wald_ci()`**: match Wilson (1927) / Newcombe (1998)
  and the textbook Wald interval; clamped to [0, 1].
- **`hanley_mcneil_var()`**: `Q1 = AUC/(2−AUC)`, `Q2 = 2AUC²/(1+AUC)`,
  variance formula matches Hanley & McNeil (1982).
- **`mc_validate_buderer()` headline (~57%)**: empirically `P_width_target
  = 0.5714`, coverage 0.931 (Wald) at B = 20000. The result is correct and
  is a general property of sample sizes targeting the *expected* CI width
  (~50% assurance), not a coding artefact. See the software report for the
  recommended framing in the manuscript.
- **`joint_sample_size()`**: the AUC component is a deterministic
  Hanley–McNeil gate (documented); the reported `joint_prob` is the Se∧Sp
  Monte Carlo probability conditional on that gate. The H1-pattern crash
  here was already fixed in 0.1.0 (`joint_prob <- 0` init); confirmed it
  returns gracefully at `prev = 0.01`.
- **`ss_adaptive_prevalence()`**: blinded re-estimation, "never reduce N"
  rule, and the per-arm `max(n_se/prev, n_sp/(1−prev))` aggregation are
  correct. The `precision_achieved` varying with the true prevalence
  (~0.70–0.99) is honest expected-width behaviour, not a bug.
- **`ss_unified()`**: the T×R 2×2 construction under conditional
  independence, the binormal AUC approximation, and the comparison table
  are coherent and produce the intended monotone "N grows with
  uncertainty" pattern.

---

## 4. Known design choices documented for the manuscript (not defects)

- **Confidence level.** `joint_sample_size()`, `ss_adaptive_prevalence()`,
  and `ss_unified()` hard-code the 95% level (`alpha = 0.05`); only
  `mc_validate_buderer()`, `ss_imperfect_ref()`, and now `ss_net_benefit()`
  expose `alpha`. 95% is the universal DTA convention, so this is a
  deliberate API simplification, not an error — but it is an inconsistency
  worth either documenting or unifying in a future minor release.
- **Default `prior_sp = c(2, 2)`** in `bam_sample_size()` / `ss_unified()`
  is *vague* (mean 0.5). It is the worst case for specificity precision and
  makes the specificity arm dominate the total N (now correctly reflected
  after C1). Users with a realistic prior should supply, e.g.,
  `prior_sp = c(18, 2)`. Recommend documenting this in the vignette.

---

## 5. Verification after fixes

- **Test suite:** `FAIL 0 | PASS 95` (was 85; coverage of the changed
  behaviour added). The 7 warnings are intentional/expected (small `B` or
  `N_range` in fast tests, plus the new informative "target not reached"
  warnings).
- **`R CMD check --as-cran`:** `Status: 2 NOTEs`, both environment-only
  (no `pandoc`/`qpdf` in the audit machine; "new submission"). No ERRORs,
  no WARNINGs from the package itself.
- **Reproducibility:** identical results on repeated calls with the same
  seed, for every changed function.
- **Vignette:** R chunks execute end-to-end with the corrected numbers
  (HTML rendering needs `pandoc`, available in the original build
  environment).

---

## 6. Independent code-review pass (second pass)

After the fixes above, the diff was put through a second, independent code
review (line-by-line, removed-behaviour, cross-file, and cleanup angles).
It surfaced four additional
*real* issues — three of them in code the first pass had just touched — and
refuted two false positives. All four were fixed and verified.

- **R1 (HIGH).** `ss_net_benefit()`: `N_conservative` silently dropped any
  *feasible* threshold that failed to reach the assurance target within
  `N_range`, taking the maximum over the easier thresholds and so reporting
  a worst-case `N` that a per-threshold warning had just flagged as
  insufficient. **Fix:** report `NA` (with a warning) when a feasible
  threshold is unmet; otherwise the max over met feasible thresholds.
- **R2 (MEDIUM, regression introduced by the first pass).**
  `ss_time_dependent_roc()`: when the target was not reached it paired
  `N_required = max(N_range)` with `prob_achieved = best_prob`, where the
  best could have occurred at a *different* `N` — an inconsistent
  `(N, prob)` pair. **Fix:** report the assurance at the largest `N` tried;
  the grid-wide best is kept in the warning.
- **R3 (LOW).** `bam_sample_size()`: `pmax(n_se/prev, n_sp/(1-prev))` can be
  non-finite if a prevalence draw reaches the 0/1 boundary (diffuse
  priors), contaminating the upper percentiles. **Fix:** drop non-finite
  draws before `median`/`quantile`.
- **R4 (test robustness).** The U-shaped net-benefit test compared fixed
  thresholds without guarding against `NA`; an `expect_false(anyNA(...))`
  guard was added so it fails cleanly rather than erroring.

**Refuted (not changed):** (i) the claim that testthat edition 3 escalates
warnings to failures — disproved by the actual run (`WARN 7 | FAIL 0`);
(ii) the claim that the net-benefit CI criterion is invalid because
treat-all net benefit is "estimated and correlated" — in this design the
prevalence split `n_d/N` is fixed, so treat-all net benefit is a known
constant and the criterion is a proper confidence-interval test of the
difference. A noted-but-not-fixed item: the per-`N` `set.seed()` (common
random numbers) is a deliberate, package-wide variance-reduction choice,
documented rather than altered.

## 7. Change inventory

| File | Change |
|------|--------|
| `R/ss_net_benefit.R` | Rewritten: CI-based criterion, `alpha` arg, feasibility flagging, crash fix, `best_prob` tracking, wider default `N_range`. |
| `R/bam_sample_size.R` | `N_total` now `pmax(n_se/prev, n_sp/(1−prev))`; docs updated. |
| `R/ss_imperfect_ref.R` | `mc_validation` reports apparent Se / true Se / bias; `@note` + Rogan–Gladen reference. |
| `R/ss_time_dependent_roc.R` | Explicit `best_prob`; warning when target not reached. |
| `R/ss_unified.R` | `@param N_range` default corrected to 1000; defensive init. |
| `tests/testthat/test-ss_net_benefit.R` | Rewritten for new behaviour (feasibility, U-shape, assurance monotonicity). |
| `tests/testthat/test-ss_imperfect_ref.R` | Now checks precision gain and bias exposure. |
| `DESCRIPTION` | Version 0.1.0 → 0.2.0. |
| `NEWS.md` | Added (full changelog). |
