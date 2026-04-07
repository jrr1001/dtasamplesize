# Audit Report: dtasamplesize v0.1.0

**Date:** 2026-04-07
**Auditor:** Claude Code (automated)
**Package:** dtasamplesize 0.1.0

---

## CRITICAL (statistical bugs)

*None found.*

All formulas match published references:
- `buderer_n()`: z^2 * p * (1-p) / d^2 matches Buderer 1996 exactly (verified: 0.85/0.07 -> 100, 0.90/0.05 -> 139)
- `wilson_ci()`: Wilson score formula matches Wilson 1927 / Newcombe 1998
- `hanley_mcneil_var()`: Q1, Q2, variance formula matches Hanley & McNeil 1982 (verified numerically)
- `ss_imperfect_ref()` VIF: 1/(Se_ref + Sp_ref - 1)^2 matches Staquet 1981 (verified: 0.90+0.95 -> 1.384)
- `bam_sample_size()`: Beta-Binomial conjugate update matches Wilson et al. 2022
- Numerical stability: P_width_target varies < 2.5% across seeds (0.558-0.571 with B=10000)
- Reproducibility: identical results with same seed confirmed

---

## HIGH (crashes / runtime errors)

### H1. `joint_sample_size()` crashes with `prev=0.01`

- **File:** `R/joint_sample_size.R`
- **Description:** When prevalence is extremely low, all N values are skipped (n_d < 2), and `joint_prob` is never assigned. The fallback at line 125 references the undefined variable.
- **Error:** `object 'joint_prob' not found`
- **Fix:** Initialize `joint_prob <- 0` before the loop.
- **Status:** FIXED

### H2. `ss_net_benefit()` no `alpha` parameter

- **File:** `R/ss_net_benefit.R`
- **Description:** Unlike other functions, `ss_net_benefit()` hardcodes alpha=0.05 implicitly (via the NB formula which doesn't use alpha). This is correct by design since NB is not a CI-based metric, but should be documented. Not a crash, but inconsistency.
- **Severity downgraded to:** MEDIUM (consistency)

---

## MEDIUM (documentation / consistency)

### M1. Parameter naming inconsistency: `d` vs `d_se`/`d_sp` vs `delta_se`/`delta_sp`

- **Files:** Multiple
- **Description:** Precision parameter names differ across functions:
  - `mc_validate_buderer`: `d` (half-width)
  - `ss_imperfect_ref`, `ss_adaptive_prevalence`: `d_se`, `d_sp`
  - `bam_sample_size`, `joint_sample_size`, `ss_unified`: `delta_se`, `delta_sp`
  - Note: `delta_se` in BAM is full width (0.14), while `d_se` in imperfect_ref is half-width (0.07)
- **Impact:** Users may confuse half-width vs full-width. Could lead to incorrect sample sizes if wrong value supplied.
- **Suggested fix:** Add clarifying `@details` note in each function stating whether delta/d is half-width or full-width.

### M2. No test files for `launch_app.R` and `print_methods.R`

- **Files:** `tests/testthat/`
- **Description:** `launch_app.R` and `print_methods.R` have no dedicated test files. `print.dtasamplesize` is indirectly tested via other tests but has no explicit test.
- **Suggested fix:** Add basic tests for print method; launch_app test would require shiny so can be skipped.

### M3. Vignette title mismatch with VignetteIndexEntry

- **File:** `vignettes/introduction.Rmd`
- **Description:** The YAML title is longer than the VignetteIndexEntry. Causes a warning during vignette building.
- **Suggested fix:** Make VignetteIndexEntry match the full title.

---

## LOW (style / conventions)

### L1. Empty `inst/shiny/www/` directory removed during build

- **Description:** The `www/` subdirectory is empty and gets removed by R CMD build. Harmless.
- **Suggested fix:** Add a `.gitkeep` or remove the empty dir.

### L2. `alpha` parameter inconsistency

- **Description:** `mc_validate_buderer` and `ss_imperfect_ref` expose `alpha` parameter. `joint_sample_size`, `ss_net_benefit`, `ss_adaptive_prevalence` hardcode alpha=0.05 internally.
- **Impact:** Minor; most DTA studies use alpha=0.05.

### L3. `doc/` and `Meta/` in .gitignore

- **Description:** Build artifacts `doc/` and `Meta/` are in .gitignore, which is correct practice.

---

## Summary

| Severity | Count | Status |
|----------|-------|--------|
| CRITICAL | 0 | - |
| HIGH | 1 | H1 FIXED |
| MEDIUM | 3 | M1 noted, M2-M3 to fix |
| LOW | 3 | Noted |

## Test results

- **75 tests pass, 0 failures**
- **9 test files covering 9 of 11 R source files** (launch_app.R and print_methods.R lack dedicated tests)
- Edge cases: 18/19 passed (1 fixed: H1)
- Numerical stability: variation < 2.5% across seeds with B=10000
- Reproducibility: exact match confirmed
