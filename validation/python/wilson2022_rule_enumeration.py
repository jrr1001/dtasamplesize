# Rule enumeration and direct simulation for the comparison with Wilson et al. (2022).
# A cross-check by the same team, written in Python with SciPy (beta quantiles from Boost,
# not Rmath), from the definitions in the article (sections 3.1.1, 3.1.3, 5.1 and 6.1 of the
# original text, PMC9325402) and, for the rule implemented in the code, from a direct reading
# of Data S2 `Paper Functions.R` (L3-44, L48-119). It is not an implementation by an
# external group and it does not run, copy or load the Data S2 code.
#
# Definitions in the article (original text):
#   3.1.1: symmetric 100(1-a)% interval (lambda_L, lambda_U) with Pr(<= L) = a/2 and Pr(>= U) = a/2;
#          success: lambda_U - lambda_L <= w*.
#   3.1.3: c1 = argmin{W(n11) >= w*} - 1, c2 = argmax{W(n11) >= w*} + 1; N = {n11 <= c1 or n11 >= c2};
#          A = 0 if success is impossible and 1 if every outcome is a success; eq. (2) sums n_T,1 = 0..n_T.
#   Search: the smallest n_T with assurance above the target (step 4 and the example in 3.1.3).
# Rules compared per arm (k = n_T,1 diseased, x = n11):
#   "def"    : success if W(k, x) <= w*  (direct definition of 3.1.1; uses no cut-offs).
#   "printed": N of 3.1.3 with F = {x : W(k, x) >= w*} (complement of the hull [min F, max F]).
#   "code"   : cut-offs and set as programmed by `cutoffs()` and `assurance()` in Data S2.
# Usage: python -I wilson2022_rule_enumeration.py <output_folder> <part: enum | mc>
# Output: wilson2022_rule_enumeration_<part>.log in <output_folder>.
import sys, os, time
import numpy as np
from scipy import stats

sys.stdout.reconfigure(encoding='utf-8')
out = sys.argv[1]; parte = sys.argv[2]
os.makedirs(out, exist_ok=True)
logf = open(os.path.join(out, f'wilson2022_rule_enumeration_{parte}.log'), 'w', encoding='utf-8')


def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s); logf.write(s + '\n'); logf.flush()


say('Rule enumeration and simulation, Wilson et al. 2022; part', parte, '; scipy', __import__('scipy').__version__,
    '; numpy', np.__version__, '; start', time.strftime('%Y-%m-%d %H:%M:%S'))


def widths(a, b, k, alpha=0.05):
    x = np.arange(k + 1)
    lo = stats.beta.ppf(alpha / 2, a + x, b + k - x)
    hi = stats.beta.ppf(1 - alpha / 2, a + x, b + k - x)
    return hi - lo


def sets_for_k(W, w, k):
    """Return three boolean vectors (def, printed, code) of length k + 1."""
    x = np.arange(k + 1)
    s_def = W <= w
    F = np.where(W >= w)[0]
    if F.size == 0:
        s_imp = np.ones(k + 1, bool)
    else:
        s_imp = (x <= F.min() - 1) | (x >= F.max() + 1)
    # Rule in the code (Data S2, Paper Functions.R): k = 0 is not summed (loop k in 1:N).
    s_cod = np.zeros(k + 1, bool)
    if k >= 1:
        fail = W > w                      # the code marks failure when width > w
        c1 = None
        for i in range(1, k + 1):         # for (i in 1:N): first i that fails
            if fail[i]:
                c1 = i - 1; break
        c2 = k                            # c2 = N
        for i in range(1, k):             # for (i in 1:(N-1)): the last i that fails
            if fail[i]:
                c2 = i + 1
        if k == 1:                        # 1:(N-1) in R is c(1, 0) when N = 1
            for i in (1, 0):
                if fail[i]:
                    c2 = i + 1
        if c1 is None:
            s_cod[:] = True               # index = 0:k
        elif c1 > 1:
            if c2 < k:
                s_cod[:c1 + 1] = True; s_cod[c2:] = True
            else:
                s_cod[:c1 + 1] = True
        elif c2 < k:
            s_cod[c2:] = True
        # otherwise, index = NA: nothing
    return s_def, s_imp, s_cod


def tables(a, b, w, kmax, alpha=0.05):
    """Conditional probability of success by k (three rules), contiguity and minimum margin."""
    P = np.zeros((3, kmax + 1)); contig = 0; noncontig = []; margin = np.inf; witnesses = {}
    for k in range(kmax + 1):
        W = widths(a, b, k, alpha)
        sd, si, sc = sets_for_k(W, w, k)
        pmf = stats.betabinom.pmf(np.arange(k + 1), k, a, b)
        P[0, k] = pmf[sd].sum(); P[1, k] = pmf[si].sum(); P[2, k] = pmf[sc].sum()
        F = np.where(~sd)[0]
        if F.size and (F.max() - F.min() + 1) != F.size:
            noncontig.append(k)
        margin = min(margin, np.min(np.abs(W - w)))
        diff = np.where(sd != sc)[0]
        if diff.size:
            witnesses[k] = [(int(x), float(W[x]), bool(sd[x]), bool(sc[x]), float(pmf[x])) for x in diff]
    return P, noncontig, margin, witnesses


def assurance(P, N, ap, bp):
    k = np.arange(N + 1)
    wk = stats.betabinom.pmf(k, N, ap, bp)
    return (wk * P[:, :N + 1]).sum(axis=1)


if parte == 'enum':
    # ---- VAP (section 5.1): Se Beta(25.9, 2.1), prevalence Beta(29, 98), w* = 0.16, 95% ----
    a, b, ap, bp, w = 25.9, 2.1, 29.0, 98.0, 0.16
    t0 = time.time(); P, nc, mg, wit = tables(a, b, w, 160)
    say(f'\n[VAP] tables k = 0..160 in {time.time() - t0:.1f} s; k with non-contiguous failure set: {nc}; minimum margin |W - w*| = {mg:.3e}')
    say('[VAP] prior-only width (k = 0):', f'{widths(a, b, 0)[0]:.6f}')
    say('  N    A_def        A_printed    A_code')
    for N in list(range(95, 116)) + [150]:
        A = assurance(P, N, ap, bp)
        say(f'  {N:3d}  {A[0]:.10f} {A[1]:.10f} {A[2]:.10f}')
    first = {r: next(N for N in range(2, 400) if assurance(P, N, ap, bp)[i] > 0.8) for i, r in enumerate(['def', 'printed', 'code'])}
    say('[VAP] smallest N with A > 0.80:', first)
    # Witnesses at N = 104: outcomes (k, x) where the definition and the code differ, weighted.
    N = 104; wk = stats.betabinom.pmf(np.arange(N + 1), N, ap, bp)
    tot = 0.0; rows = []
    for k, lst in wit.items():
        if k > N: continue
        for (x, Wx, sd, sc, px) in lst:
            contrib = wk[k] * px * (1 if sd else -1)
            tot += contrib; rows.append((abs(contrib), k, x, Wx, sd, sc, contrib))
    rows.sort(reverse=True)
    A = assurance(P, N, ap, bp)
    say(f'[VAP] N = 104: A_def - A_code = {A[0] - A[2]:.12f}; sum of the witness contributions = {tot:.12f}')
    say('  ten witnesses of largest weight (k, x, W(k, x), success by definition, success by the code, P(k) P(x | k)):')
    for r in rows[:10]:
        say(f'    k = {r[1]:3d}, x = {r[2]:3d}, W = {r[3]:.6f}, def = {r[4]}, code = {r[5]}, weight = {r[6]:+.6e}')
    cases = {'x = k omitted with c2 = k': 0.0, 'lower tail omitted (c1 <= 1)': 0.0, 'other': 0.0}
    for r in rows:
        k, x = r[1], r[2]
        if x == k and r[4] and not r[5]:
            cases['x = k omitted with c2 = k'] += r[6]
        elif x <= 1 and r[4] and not r[5]:
            cases['lower tail omitted (c1 <= 1)'] += r[6]
        else:
            cases['other'] += r[6]
    say('  split of the difference by cause:', {kk: f'{v:.6e}' for kk, v in cases.items()})

    # ---- Section 6.1: lambda = 0.9, n0 = 50; Se Beta(n0 l rho, n0 rho (1 - l)), prev Beta(n0 rho, n0 (1 - rho)), w* = 0.18 ----
    for rho, Ns, kmax in [(0.10, list(range(655, 686)), 700), (0.05, [1405, 1406, 1407, 1408, 1409, 1410, 1640, 1641, 1642, 1643, 1644], 1700)]:
        a, b, ap, bp, w = 50 * 0.9 * rho, 50 * rho * 0.1, 50 * rho, 50 * (1 - rho), 0.18
        t0 = time.time(); P, nc, mg, wit = tables(a, b, w, kmax)
        say(f'\n[6.1 rho = {rho}] Se Beta({a:g}, {b:g}), prev Beta({ap:g}, {bp:g}); tables k = 0..{kmax} in {time.time() - t0:.1f} s; '
            f'non-contiguous: {len(nc)}; minimum margin = {mg:.3e}')
        say('  N     A_def        A_printed    A_code')
        for N in Ns:
            A = assurance(P, N, ap, bp); say(f'  {N:4d}  {A[0]:.10f} {A[1]:.10f} {A[2]:.10f}')
        lim = kmax
        fd = next((N for N in range(2, lim + 1) if assurance(P, N, ap, bp)[0] > 0.8), None)
        fc = next((N for N in range(2, lim + 1) if assurance(P, N, ap, bp)[2] > 0.8), None)
        say(f'  smallest N with A > 0.80: definition {fd}; code {fc} (search up to {lim})')
    # rho = 0.01: only values on a coarse grid (2770 published)
    rho = 0.01; a, b, ap, bp, w = 50 * 0.9 * rho, 50 * rho * 0.1, 50 * rho, 50 * (1 - rho), 0.18
    t0 = time.time(); P, nc, mg, wit = tables(a, b, w, 4000)
    say(f'\n[6.1 rho = 0.01] Se Beta({a:g}, {b:g}), prev Beta({ap:g}, {bp:g}); tables k = 0..4000 in {time.time() - t0:.1f} s; non-contiguous: {len(nc)}; minimum margin = {mg:.3e}')
    for N in [1000, 2000, 2500, 2769, 2770, 2771, 3000, 3500, 4000]:
        A = assurance(P, N, ap, bp); say(f'  {N:4d}  {A[0]:.10f} {A[1]:.10f} {A[2]:.10f}')

elif parte == 'mc':
    # Direct simulation of the generative model: no enumeration, no sets, no cut-offs.
    # rho ~ Beta(ap, bp); K ~ Bin(N, rho); lambda ~ Beta(a, b); X ~ Bin(K, lambda);
    # success if the width of the symmetric 95% interval of Beta(a + X, b + K - X) is <= w*.
    rng = np.random.default_rng(20261115)

    def mc(a, b, ap, bp, w, N, R, chunk=1_000_000):
        s = 0
        for start in range(0, R, chunk):
            m = min(chunk, R - start)
            rho = rng.beta(ap, bp, m); K = rng.binomial(N, rho)
            lam = rng.beta(a, b, m); X = rng.binomial(K, lam)
            W = stats.beta.ppf(0.975, a + X, b + K - X) - stats.beta.ppf(0.025, a + X, b + K - X)
            s += int(np.count_nonzero(W <= w))
        p = s / R
        return p, np.sqrt(p * (1 - p) / R)

    for (lab, a, b, ap, bp, w, N, R) in [
            ('VAP', 25.9, 2.1, 29, 98, 0.16, 104, 20_000_000),
            ('VAP', 25.9, 2.1, 29, 98, 0.16, 105, 20_000_000),
            ('VAP', 25.9, 2.1, 29, 98, 0.16, 106, 20_000_000),
            ('6.1 rho 0.10', 4.5, 0.5, 5, 45, 0.18, 661, 10_000_000),
            ('6.1 rho 0.10', 4.5, 0.5, 5, 45, 0.18, 681, 10_000_000),
            ('6.1 rho 0.05', 2.25, 0.25, 2.5, 47.5, 0.18, 1408, 4_000_000),
            ('6.1 rho 0.05', 2.25, 0.25, 2.5, 47.5, 0.18, 1643, 4_000_000)]:
        t0 = time.time(); p, se = mc(a, b, ap, bp, w, N, R)
        say(f'[MC] {lab}, N = {N}, R = {R}: A = {p:.6f} (MCSE {se:.6f}) [{time.time() - t0:.0f} s]')

say('end', time.strftime('%Y-%m-%d %H:%M:%S'))
