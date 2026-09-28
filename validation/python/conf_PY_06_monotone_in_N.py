"""
EXPLORATORY (not pre-registered): is assurance(N) monotone in N?

Relevance: "the smallest N reaching the target" is only equal to "the first N
found by a bisection" if the curve is monotone.  My reference deliberately
scans the whole grid, so it finds the true minimum either way; this script
asks whether the two definitions can ever come apart, and whether the package
(which agreed with my full scan on all 14 scenarios) could be relying on an
assumption that fails somewhere.

Also reports P_arm(m): is the per-arm probability itself monotone in m?
"""
import json
import numpy as np
import ref_bam as R

# Declared tolerance for counting a step as a real reversal (not roundoff
# noise near P == 1, e.g. P[800] = 1.0000000000003 for prior_se=[17,3]).
# A "reversal" is a step m -> m+1 with delta = P[m+1] - P[m] < -TOL.
# Per-arm probabilities are clipped to [0, 1] before differencing so that
# floating-point overshoot above 1.0 cannot itself register as a drop.
TOL = 1e-9

grid = json.load(open("LOCKED_confirmatory_grid.json"))
SC = grid["scenarios"]
PUB = [
    {"id": "P1", "prior_se": [17, 3], "prior_sp": [2, 2], "prior_prev": [4, 16],
     "delta_se": 0.14, "delta_sp": 0.10, "level": 0.95, "N_lo": 2, "N_hi": 800},
    {"id": "P2", "prior_se": [17, 3], "prior_sp": [18, 2], "prior_prev": [4, 16],
     "delta_se": 0.14, "delta_sp": 0.10, "level": 0.95, "N_lo": 2, "N_hi": 800},
]

print(f"Tolerance TOL = {TOL:.0e}: a step m -> m+1 counts as a reversal only "
      f"if delta = P[m+1] - P[m] < -TOL. Per-arm P is clipped to [0, 1] "
      f"before differencing. m-range used per scenario is NH = min(N_hi, 900) "
      f"(the published P1/P2 cases have N_hi = 800, so m <= 800 there).")
print(f"{'id':5} {'range':>12} {'A(N) non-decr?':>15} {'#dips':>6} "
      f"{'worst dip':>12} {'first N':>8} {'P_se mono':>10} {'P_sp mono':>10}")
for s in SC + PUB:
    if s["id"] == "C12":
        continue
    NH = min(s["N_hi"], 900)
    ex = R.BAMExact(s["prior_se"], s["prior_sp"], s["prior_prev"],
                    s["delta_se"], s["delta_sp"], s["level"], nmax=NH)
    Ns = np.arange(2, NH + 1)
    A = np.clip(ex.curve(Ns), 0.0, 1.0)
    d = np.diff(A)
    dips = int((d < -TOL).sum())
    worst = float(d.min())
    # first N at which the curve dips, if any
    first = int(Ns[np.argmax(d < -TOL)] ) if dips else -1
    Pse_c = np.clip(ex.Pse, 0.0, 1.0)
    Psp_c = np.clip(ex.Psp, 0.0, 1.0)
    mse = bool(np.all(np.diff(Pse_c[1:]) >= -TOL))
    msp = bool(np.all(np.diff(Psp_c[1:]) >= -TOL))
    print(f"{s['id']:5} {f'2..{NH}':>12} {str(dips == 0):>15} {dips:>6} "
          f"{worst:>12.3e} {str(first if dips else '-'):>8} "
          f"{str(mse):>10} {str(msp):>10}")

print("\nDetail: where P_arm(m) is NOT monotone, list the first few reversals")
for s in SC[:3] + PUB:
    NH = min(s["N_hi"], 900)
    ex = R.BAMExact(s["prior_se"], s["prior_sp"], s["prior_prev"],
                    s["delta_se"], s["delta_sp"], s["level"], nmax=NH)
    for nm, P0 in (("Se", ex.Pse), ("Sp", ex.Psp)):
        P = np.clip(P0, 0.0, 1.0)
        bad = np.nonzero(np.diff(P[1:]) < -TOL)[0] + 1
        if bad.size:
            ex3 = [f"m={m}->{m+1}: {P[m]:.6f}->{P[m+1]:.6f}" for m in bad[:3]]
            print(f"  {s['id']} P_{nm}: {bad.size} reversals; " + "; ".join(ex3))
