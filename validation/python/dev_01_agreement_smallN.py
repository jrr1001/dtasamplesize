"""
DEVELOPMENT grid, step 1: do the four independent evaluators agree?

Small N (2..30) where exhaustive enumeration and exact quadrature are cheap.
This is the DEBUGGING grid; it is not the confirmatory grid.
"""
import sys, time, json
import numpy as np
import ref_bam as R

DEV_SCENARIOS = [
    # name, prior_se, prior_sp, prior_prev, delta_se, delta_sp
    ("D1 vague/vague wide-delta",  (2, 2),  (2, 2),  (6, 14),  0.55, 0.55),
    ("D2 informative Se",          (17, 3), (2, 2),  (6, 14),  0.40, 0.60),
    ("D3 skewed prevalence",       (17, 3), (18, 2), (4, 16),  0.45, 0.45),
    ("D4 near-uniform prev",       (5, 5),  (3, 7),  (1, 1),   0.50, 0.50),
    ("D5 asymmetric deltas",       (9, 1),  (2, 8),  (2, 3),   0.35, 0.65),
]

LEVEL = 0.95
NS = list(range(2, 31))

rows = []
t0 = time.time()
for name, pse, psp, pprev, dse, dsp in DEV_SCENARIOS:
    ex = R.BAMExact(pse, psp, pprev, dse, dsp, LEVEL, nmax=max(NS))
    for N in NS:
        a_cf = ex.assurance(N)
        a_q = R.quad_full(N, pse, psp, pprev, dse, dsp, LEVEL)
        a_e, ntrip = R.enum_triples(N, pse, psp, pprev, dse, dsp, LEVEL)
        rows.append(dict(scenario=name, N=N, closed=a_cf, quad=a_q,
                         enum=a_e, ntrip=ntrip))

d_qc = max(abs(r["quad"] - r["closed"]) for r in rows)
d_ec = max(abs(r["enum"] - r["closed"]) for r in rows)
d_qe = max(abs(r["quad"] - r["enum"]) for r in rows)

print(f"rows = {len(rows)}   elapsed = {time.time()-t0:.1f}s")
print(f"max |quad_full  - closed_form| = {d_qc:.3e}")
print(f"max |enum_trip  - closed_form| = {d_ec:.3e}")
print(f"max |quad_full  - enum_trip  | = {d_qe:.3e}")

# per-scenario worst
print("\nper-scenario max |quad - closed| and assurance range:")
for name, *_ in DEV_SCENARIOS:
    sub = [r for r in rows if r["scenario"] == name]
    m = max(abs(r["quad"] - r["closed"]) for r in sub)
    print(f"  {name:28s} maxdiff={m:.3e}  A(2)={sub[0]['closed']:.6f} "
          f"A(30)={sub[-1]['closed']:.6f}")

# Monte Carlo cross-check on a handful of points (MCSE-budgeted)
print("\nMonte Carlo cross-check (B=160000, MCSE~0.0012):")
for name, pse, psp, pprev, dse, dsp in DEV_SCENARIOS:
    ex = R.BAMExact(pse, psp, pprev, dse, dsp, LEVEL, nmax=30)
    for N in (10, 20, 30):
        a_cf = ex.assurance(N)
        ph, se = R.mc(N, pse, psp, pprev, dse, dsp, LEVEL, B=160000, seed=101)
        z = (ph - a_cf) / se if se > 0 else 0.0
        print(f"  {name:28s} N={N:3d} exact={a_cf:.6f} mc={ph:.6f} "
              f"mcse={se:.5f} z={z:+.2f}")

json.dump(rows, open("dev_01_agreement_smallN.json", "w"))
