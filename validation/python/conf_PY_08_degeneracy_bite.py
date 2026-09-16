"""
Is the "degenerate replication = failure" convention actually EXERCISED by the
confirmatory grid?  At N ~ 670 it is not (P(n_d in {0,N}) is ~1e-30), so the
published scenarios cannot test it.  The small-N confirmatory scenarios can.
"""
import json
import numpy as np
import ref_bam as R

grid = json.load(open("LOCKED_confirmatory_grid.json"))
SC = {s["id"]: s for s in grid["scenarios"]}
pkg = json.load(open("conf_package.json"))

print(f"{'id':5} {'N':>5} {'P(n_d=0 or N)':>15} {'A (degen=fail)':>16} "
      f"{'A (degen=ok)':>14} {'shift':>10}")
for sid in ["C09", "C10", "C11", "C05", "C07"]:
    s = SC[sid]
    N = pkg[sid]["main"]["N_total"]
    ap, bp = s["prior_prev"]
    p0 = float(np.exp(R.log_betabinom_pmf(0, N, ap, bp)))
    pN = float(np.exp(R.log_betabinom_pmf(N, N, ap, bp)))
    ex = R.BAMExact(s["prior_se"], s["prior_sp"], s["prior_prev"],
                    s["delta_se"], s["delta_sp"], s["level"], nmax=N)
    a_fail = ex.assurance(N)
    # variant: an arm of size 0 succeeds iff the PRIOR-only width already meets
    # the target (this is what "not treating degeneracy specially" would give)
    Pse2 = ex.Pse.copy(); Psp2 = ex.Psp.copy()
    Pse2[0] = 1.0 if R.ci_width(0, 0, s["prior_se"][0], s["prior_se"][1],
                                s["level"]) <= s["delta_se"] else 0.0
    Psp2[0] = 1.0 if R.ci_width(0, 0, s["prior_sp"][0], s["prior_sp"][1],
                                s["level"]) <= s["delta_sp"] else 0.0
    k = np.arange(N + 1)
    pk = np.exp(R.log_betabinom_pmf(k, N, ap, bp))
    a_ok = float(np.sum(pk * Pse2[k] * Psp2[N - k]))
    print(f"{sid:5} {N:>5} {p0+pN:>15.3e} {a_fail:>16.8f} {a_ok:>14.8f} "
          f"{a_ok-a_fail:>+10.2e}")

print("\nPublished scenarios at their selected N:")
for tag, pse, psp, N in [("P1", (17, 3), (2, 2), 678),
                         ("P2", (17, 3), (18, 2), 672)]:
    p0 = float(np.exp(R.log_betabinom_pmf(0, N, 4, 16)))
    pN = float(np.exp(R.log_betabinom_pmf(N, N, 4, 16)))
    print(f"  {tag}  N={N}  P(degenerate split) = {p0+pN:.3e}"
          "   -> convention is untestable at this N")
