"""
Independent reference implementation of the BAM joint-assurance estimand.

Written from the mathematical specification only (Wilson et al. 2022 style
Bayesian assurance for a cohort DTA study).  No code from the R package
`dtasamplesize` was read or copied; only its .Rd documentation was consulted
to fix the parameterisation (delta = FULL width, alpha_ci = equal-tailed
credible level, degenerate arm = failure).

MODEL
-----
    prev  ~ Beta(a_p , b_p)
    n_d   ~ Binomial(N, prev)          n_nd = N - n_d
    Se    ~ Beta(a_se, b_se)      x ~ Binomial(n_d , Se)
    Sp    ~ Beta(a_sp, b_sp)      y ~ Binomial(n_nd, Sp)

    posterior Se | x  =  Beta(a_se + x, b_se + n_d - x)
    posterior Sp | y  =  Beta(a_sp + y, b_sp + n_nd - y)

    w(x,n;a,b) = qbeta(1-e/2 ; a+x, b+n-x) - qbeta(e/2 ; a+x, b+n-x),  e = 1-level

    success  <=>  n_d >= 1  AND  n_nd >= 1
                  AND w(x,n_d ;a_se,b_se) <= delta_se
                  AND w(y,n_nd;a_sp,b_sp) <= delta_sp

    assurance(N) = P(success)

Four INDEPENDENT evaluators of assurance(N) are provided:

  closed_form(N)    -- Beta-Binomial marginalisation + conditional
                       independence factorisation  (the analytic route)
  quad_full(N)      -- full 3-D Gauss-Jacobi quadrature over
                       (prev, Se, Sp) of the *raw* generative model.
                       The integrand is a POLYNOMIAL in each variable of
                       degree <= N, so Gauss-Jacobi with >= (N+1)/2 nodes
                       is exact to machine precision.  This route never
                       uses the Beta-Binomial pmf at all.
  enum_triples(N)   -- exhaustive enumeration of every (n_d, x, y) triple
                       with its Beta-Binomial mass, no factorisation.
  mc(N, B, seed)    -- direct Monte Carlo of the generative model.
"""

import numpy as np
from scipy.special import gammaln, roots_jacobi, betaincinv

LOG0 = -np.inf


# --------------------------------------------------------------------------
# 1. credible interval width
# --------------------------------------------------------------------------
def ci_width(x, n, a, b, level=0.95):
    """Equal-tailed posterior credible interval FULL width for Beta(a+x, b+n-x).

    betaincinv is the regularised incomplete beta inverse, i.e. the Beta
    quantile function -- used directly rather than scipy.stats.beta.ppf so
    that the computation is one call into the underlying special function.
    """
    x = np.asarray(x, dtype=float)
    n = np.asarray(n, dtype=float)
    e = 1.0 - level
    ap = a + x
    bp = b + n - x
    return betaincinv(ap, bp, 1.0 - e / 2.0) - betaincinv(ap, bp, e / 2.0)


# --------------------------------------------------------------------------
# 2. Beta-Binomial log pmf
# --------------------------------------------------------------------------
def log_betabinom_pmf(x, n, a, b):
    """log P(X = x), X ~ BetaBinomial(n, a, b)."""
    x = np.asarray(x, dtype=float)
    n = np.asarray(n, dtype=float)
    return (
        gammaln(n + 1.0) - gammaln(x + 1.0) - gammaln(n - x + 1.0)
        + gammaln(x + a) + gammaln(n - x + b) - gammaln(n + a + b)
        + gammaln(a + b) - gammaln(a) - gammaln(b)
    )


# --------------------------------------------------------------------------
# 3. per-arm success indicator matrix and P_arm(m)
# --------------------------------------------------------------------------
def arm_ok_matrix(nmax, a, b, delta, level=0.95):
    """Boolean (nmax+1) x (nmax+1) matrix M[m, x] = 1[w(x,m) <= delta].

    Entries with x > m are meaningless and set False.
    Row m = 0 is left all-False: an arm of size zero is a DEGENERATE
    replication and counts as a FAILURE regardless of the prior-only width.
    """
    m = np.arange(nmax + 1)[:, None].astype(float)
    x = np.arange(nmax + 1)[None, :].astype(float)
    valid = x <= m
    xs = np.where(valid, x, 0.0)
    w = ci_width(xs, m, a, b, level)
    ok = valid & (w <= delta)
    ok[0, :] = False          # degeneracy convention
    return ok


def arm_prob(nmax, a, b, delta, level=0.95):
    """P_arm(m) for m = 0..nmax.  P_arm(0) = 0."""
    ok = arm_ok_matrix(nmax, a, b, delta, level)
    out = np.zeros(nmax + 1)
    for m in range(1, nmax + 1):
        idx = np.nonzero(ok[m, : m + 1])[0]
        if idx.size:
            out[m] = np.exp(log_betabinom_pmf(idx, m, a, b)).sum()
    return out


# --------------------------------------------------------------------------
# 4. evaluator A: closed form
# --------------------------------------------------------------------------
class BAMExact:
    """Closed-form assurance with cached per-arm probabilities."""

    def __init__(self, prior_se, prior_sp, prior_prev,
                 delta_se, delta_sp, level=0.95, nmax=1000):
        self.a_se, self.b_se = float(prior_se[0]), float(prior_se[1])
        self.a_sp, self.b_sp = float(prior_sp[0]), float(prior_sp[1])
        self.a_p, self.b_p = float(prior_prev[0]), float(prior_prev[1])
        self.delta_se, self.delta_sp = float(delta_se), float(delta_sp)
        self.level = float(level)
        self.nmax = int(nmax)
        self.Pse = arm_prob(self.nmax, self.a_se, self.b_se,
                            self.delta_se, self.level)
        self.Psp = arm_prob(self.nmax, self.a_sp, self.b_sp,
                            self.delta_sp, self.level)

    def assurance(self, N):
        N = int(N)
        if N > self.nmax:
            raise ValueError(f"N={N} exceeds cache nmax={self.nmax}")
        k = np.arange(N + 1)
        pk = np.exp(log_betabinom_pmf(k, N, self.a_p, self.b_p))
        return float(np.sum(pk * self.Pse[k] * self.Psp[N - k]))

    def curve(self, Ns):
        return np.array([self.assurance(N) for N in Ns])

    def search(self, target, N_lo=2, N_hi=None):
        """Smallest N in [N_lo, N_hi] with assurance(N) >= target.

        Evaluated by an exhaustive scan of the grid (no monotonicity
        assumption is made about assurance(N)).
        """
        if N_hi is None:
            N_hi = self.nmax
        for N in range(int(N_lo), int(N_hi) + 1):
            if self.assurance(N) >= target:
                return N, self.assurance(N)
        return None, None


# --------------------------------------------------------------------------
# 5. evaluator B: full 3-D Gauss-Jacobi quadrature of the raw model
# --------------------------------------------------------------------------
def _beta_nodes(a, b, n):
    """Gauss-Jacobi nodes/weights integrating f(p) against the Beta(a,b) pdf.

    Beta(a,b) density on (0,1) is p^(a-1)(1-p)^(b-1)/B(a,b).
    Substituting p = (1+t)/2 maps to the Jacobi weight
    (1-t)^(b-1) (1+t)^(a-1) on (-1,1).  roots_jacobi(n, alpha, beta) uses
    weight (1-t)^alpha (1+t)^beta, so alpha = b-1, beta = a-1.
    Weights are renormalised to sum to 1, which is exactly the Beta
    normalising constant.  Exact for polynomial f of degree <= 2n-1.
    """
    t, w = roots_jacobi(n, b - 1.0, a - 1.0)
    p = (1.0 + t) / 2.0
    w = w / w.sum()
    return p, w


def _log_binom_coef(n, k):
    return gammaln(n + 1.0) - gammaln(k + 1.0) - gammaln(n - k + 1.0)


def quad_full(N, prior_se, prior_sp, prior_prev, delta_se, delta_sp,
              level=0.95, extra_nodes=8, ok_se=None, ok_sp=None):
    """assurance(N) by exact 3-D quadrature of the generative model.

    No Beta-Binomial pmf is used anywhere in this routine.
    """
    a_se, b_se = float(prior_se[0]), float(prior_se[1])
    a_sp, b_sp = float(prior_sp[0]), float(prior_sp[1])
    a_p, b_p = float(prior_prev[0]), float(prior_prev[1])
    nq = int(np.ceil((N + 1) / 2.0)) + int(extra_nodes)

    p, wp = _beta_nodes(a_p, b_p, nq)
    s, ws = _beta_nodes(a_se, b_se, nq)
    q, wq = _beta_nodes(a_sp, b_sp, nq)

    if ok_se is None:
        ok_se = arm_ok_matrix(N, a_se, b_se, delta_se, level)
    if ok_sp is None:
        ok_sp = arm_ok_matrix(N, a_sp, b_sp, delta_sp, level)

    k = np.arange(N + 1)

    # Binom(k; N, p)  ->  (N+1, nq)
    with np.errstate(divide="ignore", invalid="ignore"):
        logp = np.where(p > 0, np.log(p), LOG0)
        log1p_ = np.where(p < 1, np.log1p(-p), LOG0)
    lb = (_log_binom_coef(float(N), k.astype(float))[:, None]
          + k[:, None] * logp[None, :]
          + (N - k)[:, None] * log1p_[None, :])
    Bk = np.exp(lb)                                        # (N+1, nq)

    # g_se[k, j] = sum_x Binom(x; k, s_j) 1[w(x,k) <= delta_se]
    def g_matrix(ok, theta):
        with np.errstate(divide="ignore", invalid="ignore"):
            lt = np.where(theta > 0, np.log(theta), LOG0)
            l1t = np.where(theta < 1, np.log1p(-theta), LOG0)
        G = np.zeros((N + 1, theta.size))
        for m in range(0, N + 1):
            idx = np.nonzero(ok[m, : m + 1])[0]
            if idx.size == 0:
                continue
            xf = idx.astype(float)
            lg = (_log_binom_coef(float(m), xf)[:, None]
                  + xf[:, None] * lt[None, :]
                  + (m - xf)[:, None] * l1t[None, :])
            G[m] = np.exp(lg).sum(axis=0)
        return G

    Gse = g_matrix(ok_se, s)          # (N+1, nq) indexed by arm size k
    Gsp = g_matrix(ok_sp, q)          # (N+1, nq) indexed by arm size N-k

    # full 3-D integral, NOT factorised across k:
    #   sum_j sum_l sum_m wp_j ws_l wq_m * sum_k Bk[k,j] Gse[k,l] Gsp[N-k,m]
    # done as an einsum over (k, j, l, m)
    total = np.einsum("kj,kl,km,j,l,m->", Bk, Gse, Gsp[::-1], wp, ws, wq,
                      optimize=True)
    return float(total)


# --------------------------------------------------------------------------
# 6. evaluator C: exhaustive enumeration of (n_d, x, y) triples
# --------------------------------------------------------------------------
def enum_triples(N, prior_se, prior_sp, prior_prev, delta_se, delta_sp,
                 level=0.95):
    """assurance(N) by summing the mass of every successful (n_d, x, y)."""
    a_se, b_se = float(prior_se[0]), float(prior_se[1])
    a_sp, b_sp = float(prior_sp[0]), float(prior_sp[1])
    a_p, b_p = float(prior_prev[0]), float(prior_prev[1])
    total = 0.0
    ntrip = 0
    for k in range(0, N + 1):
        m = N - k
        pk = np.exp(log_betabinom_pmf(k, N, a_p, b_p))
        if k == 0 or m == 0:
            continue                       # degenerate -> failure
        for x in range(0, k + 1):
            if ci_width(x, k, a_se, b_se, level) > delta_se:
                continue
            px = np.exp(log_betabinom_pmf(x, k, a_se, b_se))
            for y in range(0, m + 1):
                if ci_width(y, m, a_sp, b_sp, level) > delta_sp:
                    continue
                py = np.exp(log_betabinom_pmf(y, m, a_sp, b_sp))
                total += pk * px * py
                ntrip += 1
    return float(total), ntrip


# --------------------------------------------------------------------------
# 7. evaluator D: direct Monte Carlo of the generative model
# --------------------------------------------------------------------------
def mc(N, prior_se, prior_sp, prior_prev, delta_se, delta_sp,
       level=0.95, B=40000, seed=1, chunk=200000):
    """Monte Carlo assurance.  Returns (phat, mcse)."""
    a_se, b_se = float(prior_se[0]), float(prior_se[1])
    a_sp, b_sp = float(prior_sp[0]), float(prior_sp[1])
    a_p, b_p = float(prior_prev[0]), float(prior_prev[1])
    rng = np.random.default_rng(seed)
    ok_se = arm_ok_matrix(N, a_se, b_se, delta_se, level)
    ok_sp = arm_ok_matrix(N, a_sp, b_sp, delta_sp, level)
    hits = 0
    done = 0
    while done < B:
        nb = min(chunk, B - done)
        prev = rng.beta(a_p, b_p, nb)
        se = rng.beta(a_se, b_se, nb)
        sp = rng.beta(a_sp, b_sp, nb)
        nd = rng.binomial(N, prev)
        nnd = N - nd
        x = rng.binomial(nd, se)
        y = rng.binomial(nnd, sp)
        good = ok_se[nd, x] & ok_sp[nnd, y]
        hits += int(good.sum())
        done += nb
    phat = hits / B
    return phat, float(np.sqrt(phat * (1 - phat) / B))
