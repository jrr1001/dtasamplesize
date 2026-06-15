# External Cross-Validation Report — dtasamplesize 0.2.0

**Date:** 2026-06-15
**Goal:** verify the core numerical machinery of `dtasamplesize` against
*independent* references — not by re-reading the code, but by comparing its
output to established R implementations, numerical integration, and Monte
Carlo ground truth. Reproduction script: `validation/cross_validation.R`
(run under R 4.5.2 with `Hmisc`, `pROC`).

This pass is qualitatively different from the code audit: it asks "are the
formulas numerically correct against an outside source?", which internal
review cannot answer.

---

## Summary

| # | Quantity | Reference | Result |
|---|----------|-----------|--------|
| V1 | Wilson score CI | base `prop.test(correct=FALSE)` | **Exact match** (all cases) |
| V2 | Wald CI | `Hmisc::binconf("asymptotic")` | **Exact match** |
| V3 | `buderer_n` | published values (100, 139) | **Exact** |
| V4 | Hanley–McNeil AUC variance | Monte Carlo (binormal) | Match for balanced; conservative for imbalanced (expected) |
| V5 | Net benefit point estimate | Vickers (2006) definition | **Exact identity** |
| V6 | **Net benefit variance (new in 0.2.0)** | Monte Carlo SD | **Match, ratio 0.995–0.998** |
| V7 | Beta–Binomial posterior CI (BAM) | numerical integration | **Match to 1e-3** |
| V8 | Imperfect-ref VIF & apparent Se | Rogan–Gladen / closed form | **Exact** |

No defect was found. Two "non-matches" were investigated and are
**reference-side or documented-approximation effects, not package errors**
(V1 boundary, V4 imbalance) — explained below.

---

## Detail and interpretation

**V1 — Wilson CI.** `wilson_ci()` reproduces base R's `prop.test(x, n,
correct = FALSE)$conf.int` (the canonical Wilson score interval) to 1e-6 for
x/n ∈ {85/100, 45/50, 2/20, 19/20, 0/30}. At the boundary case 19/20,
`Hmisc::binconf("wilson")` returns a different upper limit (0.99744 vs
0.99112); base R and `dtasamplesize` agree (0.99112). The divergence is in
Hmisc's boundary handling, not in the package — validated against the
canonical implementation.

**V2 — Wald CI.** Exact agreement with `Hmisc::binconf("asymptotic")` on
interior cases. (`wald_ci()` additionally clamps to [0, 1], which only
matters near the boundary.)

**V4 — Hanley–McNeil AUC variance.** The formula matches the Monte Carlo
variance of the empirical AUC within ~4% for balanced designs (100/100,
80/80) but is ~30% larger for an imbalanced, high-AUC design (60/140,
AUC 0.85). This is expected: the Hanley–McNeil (1982) variance uses
`Q1 = AUC/(2−AUC)`, `Q2 = 2·AUC²/(1+AUC)`, derived under a
negative-exponential score model; for binormal data the true Q-terms differ,
and the gap grows with class imbalance and AUC. For **sample-size planning
this is conservative** (it predicts more variance ⇒ a larger, safe N), so it
errs in the protective direction. The implementation faithfully reproduces
the cited 1982 formula; the approximation is documented as a limitation
rather than changed.

**V6 — Net benefit variance.** The closed-form standard error introduced in
v0.2.0,
`SE(NB) = sqrt((n_d/N)²·Se(1−Se)/n_d + (n_nd/N)²·(pt/(1−pt))²·Sp(1−Sp)/n_nd)`,
matches the Monte Carlo standard deviation of the simulated net benefit to
within 0.5% (ratios 0.995 and 0.998 across two scenarios). This confirms the
new criterion's variance is correct — the key validation of the rewritten
`ss_net_benefit()`.

**V7 — Beta–Binomial conjugacy.** The posterior credible-interval width used
by `bam_sample_size()` (`qbeta` on `Beta(a+x, b+n−x)`) equals the width from
a 200,001-point numerical posterior (likelihood × prior, normalised) to
1e-3, confirming the conjugate update.

**V8 — Imperfect reference.** The VIF equals the Rogan–Gladen value
`1/(Se_ref+Sp_ref−1)²` exactly; the simulated apparent sensitivity
(0.7642) matches the closed-form `P(T+|R+)` (0.7639).

---

## Methodological defense of the two novel components

A peer reviewer will scrutinise the components without a single canonical
precedent. The defense:

### `ss_net_benefit()` — sample size for net benefit

- **Net benefit is standard** (Vickers & Elkin 2006); only its use *as a
  sample-size target* is novel. The estimand is exactly the established one
  (V5), and its sampling variance is derived analytically and confirmed
  against Monte Carlo (V6).
- **The criterion is inference-based**, not a point-estimate heuristic: a
  study "succeeds" when the lower (1−α) confidence limit of net benefit
  exceeds both default strategies. Because the prevalence split `n_d/N` is
  fixed by design, the treat-all net benefit is a known constant, so the
  comparison is a proper one-sided CI test of the difference (this is why
  the criterion is valid — a point a reviewer might otherwise question).
- **Feasibility is explicit:** thresholds where the test cannot beat a
  default strategy are flagged, not silently sized.
- *Honest limitation to state:* the normal approximation to the net-benefit
  distribution degrades at very small arm sizes; report the achieved
  assurance and avoid extrapolating below a few events per arm.

### `ss_imperfect_ref()` — imperfect reference standard

- The inflation factor is the **Rogan–Gladen (1978)** misclassification
  variance multiplier (V8), a published, citable quantity, applied under the
  conditional-independence assumption shared with the Staquet (1981)
  framework.
- *Honest limitation, now surfaced by the package itself:* the factor
  restores **precision** but not **unbiasedness** — the apparent sensitivity
  estimated against an imperfect reference is biased (≈0.76 vs true 0.85 in
  the default scenario; the `mc_validation` table reports this bias).
  Unbiased estimation requires a bias-correction or latent-class analysis at
  the analysis stage. The package should be presented as a *planning* tool
  that sizes for adequate precision, paired with an appropriate analysis
  method.

---

## Bottom line

The numerical core is correct against independent references; the new
net-benefit variance is validated to <1%; the only deviations are a
reference-package boundary quirk and the documented conservativeness of the
Hanley–McNeil approximation. The residual work for publication is
*methodological framing* (the two points above), not further bug-hunting.
