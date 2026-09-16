"""
NEGATIVE CONTROLS / limits of the validation.

My reference and the package share, by construction, the SPECIFICATION of the
estimand.  This script quantifies how far the answer moves under plausible
alternative readings of that specification -- i.e. how large an error my
agreement with the package would NOT detect, because I would have made the
same modelling choice.

Variants probed on the two published scenarios (target 0.80):
  V0  equal-tailed CI, width <= delta            (the assumed convention)
  V1  equal-tailed CI, width <  delta            (strict inequality)
  V2  HPD (shortest) credible interval           (different interval type)
  V3  delta read as HALF width                   (parameterisation slip)
  V4  degenerate arm allowed to succeed on the prior-only width
  V5  n_d fixed at round(N*E[prev]) instead of Binomial(N, prev)
      (i.e. ignoring the randomness of the cohort split)
"""
import numpy as np
from scipy.special import betaincinv
import ref_bam as R

LEVEL, TARGET = 0.95, 0.80
NMAX = 900
CASES = [("P1", (17, 3), (2, 2), (4, 16), 0.14, 0.10, 678),
         ("P2", (17, 3), (18, 2), (4, 16), 0.14, 0.10, 672)]


def widths_et(nmax, a, b):
    m = np.arange(nmax + 1)[:, None].astype(float)
    x = np.arange(nmax + 1)[None, :].astype(float)
    valid = x <= m
    xs = np.where(valid, x, 0.0)
    return R.ci_width(xs, m, a, b, LEVEL), valid


def widths_hpd(nmax, a, b, iters=70):
    """Shortest-interval width for every (m, x), by vectorised bisection on
    the lower tail probability lp in [0, 1-LEVEL]."""
    m = np.arange(nmax + 1)[:, None].astype(float)
    x = np.arange(nmax + 1)[None, :].astype(float)
    valid = x <= m
    xs = np.where(valid, x, 0.0)
    ap = a + xs
    bp = b + m - xs
    lo_p = np.zeros_like(ap)
    hi_p = np.full_like(ap, 1.0 - LEVEL)
    logd = lambda t, A, Bp: (A - 1) * np.log(np.clip(t, 1e-300, 1)) \
                            + (Bp - 1) * np.log1p(-np.clip(t, 0, 1 - 1e-16))
    for _ in range(iters):
        mid = 0.5 * (lo_p + hi_p)
        lo = betaincinv(ap, bp, np.clip(mid, 1e-15, 1))
        hi = betaincinv(ap, bp, np.clip(mid + LEVEL, 0, 1 - 1e-15))
        # density at lo minus density at hi; HPD when equal
        f = logd(lo, ap, bp) - logd(hi, ap, bp)
        go_up = f < 0                      # density still lower at lo -> move right
        lo_p = np.where(go_up, mid, lo_p)
        hi_p = np.where(go_up, hi_p, mid)
    mid = 0.5 * (lo_p + hi_p)
    w = betaincinv(ap, bp, np.clip(mid + LEVEL, 0, 1 - 1e-15)) \
        - betaincinv(ap, bp, np.clip(mid, 1e-15, 1))
    # monotone posteriors: the shortest interval is one-sided
    w = np.minimum(w, betaincinv(ap, bp, LEVEL))          # lower-anchored
    w = np.minimum(w, 1.0 - betaincinv(ap, bp, 1 - LEVEL))  # upper-anchored
    return w, valid


def arm_prob_from_ok(ok, a, b, nmax):
    out = np.zeros(nmax + 1)
    for m in range(0, nmax + 1):
        idx = np.nonzero(ok[m, : m + 1])[0]
        if idx.size:
            out[m] = np.exp(R.log_betabinom_pmf(idx, m, a, b)).sum()
    return out


def assur(N, Pse, Psp, ap, bp):
    k = np.arange(N + 1)
    pk = np.exp(R.log_betabinom_pmf(k, N, ap, bp))
    return float(np.sum(pk * Pse[k] * Psp[N - k]))


def search(Pse, Psp, ap, bp, nmax):
    for N in range(2, nmax + 1):
        if assur(N, Pse, Psp, ap, bp) >= TARGET:
            return N
    return None


print(f"{'case':5} {'variant':44} {'N*':>7} {'shift':>8}  note")
for tag, pse, psp, pprev, dse, dsp, Nclaim in CASES:
    w_se, v_se = widths_et(NMAX, *pse)
    w_sp, v_sp = widths_et(NMAX, *psp)
    h_se, _ = widths_hpd(NMAX, *pse)
    h_sp, _ = widths_hpd(NMAX, *psp)

    def run(ok_se, ok_sp, zero_ok_se=False, zero_ok_sp=False, fixed=False):
        ok_se = ok_se.copy(); ok_sp = ok_sp.copy()
        if not zero_ok_se:
            ok_se[0, :] = False
        if not zero_ok_sp:
            ok_sp[0, :] = False
        Pse = arm_prob_from_ok(ok_se, pse[0], pse[1], NMAX)
        Psp = arm_prob_from_ok(ok_sp, psp[0], psp[1], NMAX)
        if fixed:
            Ep = pprev[0] / (pprev[0] + pprev[1])
            for n in range(2, NMAX + 1):
                k = int(round(n * Ep))
                if Pse[k] * Psp[n - k] >= TARGET:
                    return n, Pse, Psp
            return None, Pse, Psp
        return search(Pse, Psp, pprev[0], pprev[1], NMAX), Pse, Psp

    base = None
    rows = [
        ("V0 equal-tailed, w <= delta  [ASSUMED]",
         v_se & (w_se <= dse), v_sp & (w_sp <= dsp), dict()),
        ("V1 equal-tailed, w <  delta",
         v_se & (w_se < dse), v_sp & (w_sp < dsp), dict()),
        ("V2 HPD (shortest) interval",
         v_se & (h_se <= dse), v_sp & (h_sp <= dsp), dict()),
        ("V4 degenerate arm may succeed",
         v_se & (w_se <= dse), v_sp & (w_sp <= dsp),
         dict(zero_ok_se=True, zero_ok_sp=True)),
        ("V5 n_d fixed at round(N*E[prev])",
         v_se & (w_se <= dse), v_sp & (w_sp <= dsp), dict(fixed=True)),
    ]
    for name, oks, okp, kw in rows:
        N, Pse, Psp = run(oks, okp, **kw)
        if base is None:
            base = N
        sh = "-" if (N is None or base is None) else f"{N-base:+d}"
        print(f"{tag:5} {name:44} {str(N):>7} {sh:>8}")
    # V3: delta as half width -- report assurance at the claimed N instead of
    # searching, because N* is far outside the cache.
    ok3s = v_se & (w_se <= dse / 2); ok3s[0, :] = False
    ok3p = v_sp & (w_sp <= dsp / 2); ok3p[0, :] = False
    P3s = arm_prob_from_ok(ok3s, pse[0], pse[1], NMAX)
    P3p = arm_prob_from_ok(ok3p, psp[0], psp[1], NMAX)
    a3 = assur(Nclaim, P3s, P3p, pprev[0], pprev[1])
    print(f"{tag:5} {'V3 delta treated as HALF width':44} {'>'+str(NMAX):>7} "
          f"{'':>8}  assurance at N={Nclaim} is only {a3:.6f}")
    print()
