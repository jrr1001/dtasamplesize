"""
Extension of the exact BAM assurance reference to SEPARATE design and analysis
priors, to the one-sided width criterion and to marginal assurance objectives.

WHAT THIS DOES
--------------
Computes, in closed form (Beta-Binomial enumeration, no simulation), the
assurance A(N) of a cohort DTA study in which

  * the DESIGN priors (prior_se, prior_sp, prior_prev) generate the data:
        n_d ~ BetaBinomial(N, prior_prev),
        x | n_d ~ BetaBinomial(n_d, prior_se),  y | n_nd ~ BetaBinomial(n_nd, prior_sp);
  * the ANALYSIS priors (analysis_prior_se, analysis_prior_sp; default = design
    priors) define the posterior of each arm, Beta(a_A + x, b_A + n - x),
    whose width is compared with the target (width <= delta, no tolerance);
  * prevalence has no analysis prior: prior_prev acts only in the design;
  * width_type = "full_width"        : equal-tailed two-sided width at `level`;
    width_type = "lower_half_width"  : posterior median minus the ONE-sided
                                       lower bound qbeta(1 - level);
  * assurance_for = "joint": A(N) = sum_k P(n_d=k) P_se(k) P_sp(N-k)
                    "se"   : A(N) = sum_k P(n_d=k) P_se(k)
                    "sp"   : A(N) = sum_k P(n_d=k) P_sp(N-k)
    with the empty-arm convention P_se(0) = P_sp(0) = 0 (an empty arm is a
    failure, whatever the prior-only width).

REUSE OF ref_bam
----------------
`ref_bam` is imported unmodified (its folder is put on sys.path).  Reused:
`ref_bam.ci_width` (full-width quantile function) and
`ref_bam.log_betabinom_pmf`.  The per-arm matrix/probability routines of
`ref_bam` assume design prior == analysis prior and the full-width criterion,
so `arm_ok_matrix_sep` / `arm_prob_sep` / `BAMExactSep` below are new code.

PROVENANCE
----------
Written by an LLM-based coding assistant under the supervision of the
authors, from the Rd pages and NEWS of dtasamplesize 0.7.0
(`man/bam_sample_size.Rd`, `man/bam_prior_from_pilot.Rd`, `NEWS.md`) and from
`validation/python/ref_bam.py`, WITHOUT reading the R code of the package.
It is a check by the same team, not an implementation by an external group.
"""

import os
import sys

import numpy as np
from scipy.special import betaincinv

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import ref_bam  # noqa: E402  (imported unmodified)


# --------------------------------------------------------------------------
# 1. width criteria
# --------------------------------------------------------------------------
def width_sep(x, n, a_A, b_A, level=0.95, width_type="full_width"):
    """Width of the posterior Beta(a_A + x, b_A + n - x) under the chosen criterion."""
    if width_type == "full_width":
        return ref_bam.ci_width(x, n, a_A, b_A, level)
    if width_type == "lower_half_width":
        x = np.asarray(x, dtype=float)
        n = np.asarray(n, dtype=float)
        ap = a_A + x
        bp = b_A + n - x
        return betaincinv(ap, bp, 0.5) - betaincinv(ap, bp, 1.0 - level)
    raise ValueError("width_type must be 'full_width' or 'lower_half_width'")


# --------------------------------------------------------------------------
# 2. per-arm success probability with separate priors
# --------------------------------------------------------------------------
def arm_ok_matrix_sep(nmax, a_A, b_A, delta, level=0.95, width_type="full_width"):
    """M[m, x] = 1[w(x, m; analysis prior) <= delta]; row m = 0 all False."""
    m = np.arange(nmax + 1)[:, None].astype(float)
    x = np.arange(nmax + 1)[None, :].astype(float)
    valid = x <= m
    xs = np.where(valid, x, 0.0)
    w = width_sep(xs, m, a_A, b_A, level, width_type)
    ok = valid & (w <= delta)
    ok[0, :] = False           # empty arm = failure
    return ok


def arm_prob_sep(nmax, design, analysis, delta, level=0.95, width_type="full_width"):
    """P_arm(m), m = 0..nmax: Beta-Binomial pmf under the DESIGN prior summed
    over counts whose ANALYSIS-prior posterior width meets delta. P_arm(0) = 0."""
    a_D, b_D = float(design[0]), float(design[1])
    a_A, b_A = float(analysis[0]), float(analysis[1])
    ok = arm_ok_matrix_sep(nmax, a_A, b_A, delta, level, width_type)
    out = np.zeros(nmax + 1)
    for m in range(1, nmax + 1):
        idx = np.nonzero(ok[m, : m + 1])[0]
        if idx.size:
            out[m] = np.exp(ref_bam.log_betabinom_pmf(idx, m, a_D, b_D)).sum()
    return out


# --------------------------------------------------------------------------
# 3. assurance
# --------------------------------------------------------------------------
class BAMExactSep:
    def __init__(self, prior_se, prior_sp, prior_prev, delta_se, delta_sp,
                 level=0.95, nmax=1000, analysis_prior_se=None,
                 analysis_prior_sp=None, width_type="full_width",
                 assurance_for="joint"):
        if assurance_for not in ("joint", "se", "sp"):
            raise ValueError("assurance_for must be 'joint', 'se' or 'sp'")
        self.prior_se = tuple(float(v) for v in prior_se)
        self.prior_sp = tuple(float(v) for v in prior_sp)
        self.prior_prev = tuple(float(v) for v in prior_prev)
        self.an_se = (tuple(float(v) for v in analysis_prior_se)
                      if analysis_prior_se is not None else self.prior_se)
        self.an_sp = (tuple(float(v) for v in analysis_prior_sp)
                      if analysis_prior_sp is not None else self.prior_sp)
        self.delta_se, self.delta_sp = float(delta_se), float(delta_sp)
        self.level = float(level)
        self.nmax = int(nmax)
        self.width_type = width_type
        self.assurance_for = assurance_for
        # arms that do not enter the objective are not computed
        zero = np.zeros(self.nmax + 1)
        self.Pse = (arm_prob_sep(self.nmax, self.prior_se, self.an_se,
                                 self.delta_se, self.level, width_type)
                    if assurance_for in ("joint", "se") else zero)
        self.Psp = (arm_prob_sep(self.nmax, self.prior_sp, self.an_sp,
                                 self.delta_sp, self.level, width_type)
                    if assurance_for in ("joint", "sp") else zero)

    def assurance(self, N):
        N = int(N)
        if N > self.nmax:
            raise ValueError(f"N={N} exceeds cache nmax={self.nmax}")
        k = np.arange(N + 1)
        a_p, b_p = self.prior_prev
        pk = np.exp(ref_bam.log_betabinom_pmf(k, N, a_p, b_p))
        if self.assurance_for == "joint":
            return float(np.sum(pk * self.Pse[k] * self.Psp[N - k]))
        if self.assurance_for == "se":
            return float(np.sum(pk * self.Pse[k]))
        return float(np.sum(pk * self.Psp[N - k]))

    def search(self, target, N_lo=2, N_hi=None):
        """Exhaustive ascending scan. Returns dict with target_reached,
        N_total, joint_assurance (A at N_total), max_assurance_evaluated and
        N_at_max_assurance (first N attaining the maximum)."""
        if N_hi is None:
            N_hi = self.nmax
        best, best_N = -np.inf, None
        for N in range(int(N_lo), int(N_hi) + 1):
            a = self.assurance(N)
            if a > best:
                best, best_N = a, N
            if a >= target:
                return dict(target_reached=True, N_total=N, assurance=a,
                            max_assurance_evaluated=a, N_at_max=N)
        return dict(target_reached=False, N_total=None, assurance=None,
                    max_assurance_evaluated=best, N_at_max=best_N)


# --------------------------------------------------------------------------
# 4. power-prior update (bam_prior_from_pilot.Rd)
# --------------------------------------------------------------------------
def prior_from_pilot(tp, fn, tn, fp, base_se, base_sp, base_prev, a0=1.0,
                     design="cohort", use_in="both"):
    """Returns dict with design priors, analysis priors (as in the Rd)."""
    if design not in ("cohort", "case_control"):
        raise ValueError("design")
    if use_in not in ("both", "design"):
        raise ValueError("use_in")
    se = (base_se[0] + a0 * tp, base_se[1] + a0 * fn)
    sp = (base_sp[0] + a0 * tn, base_sp[1] + a0 * fp)
    if design == "cohort":
        prev = (base_prev[0] + a0 * (tp + fn), base_prev[1] + a0 * (tn + fp))
    else:
        prev = (float(base_prev[0]), float(base_prev[1]))
    if use_in == "both":
        an_se, an_sp = se, sp
    else:
        an_se, an_sp = (float(base_se[0]), float(base_se[1])), \
                       (float(base_sp[0]), float(base_sp[1]))
    return dict(prior_se=se, prior_sp=sp, prior_prev=prev,
                analysis_prior_se=an_se, analysis_prior_sp=an_sp)
