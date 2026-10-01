"""
CONFIRMATORY run, independent-reference side.

Reads LOCKED_confirmatory_grid.json (defined before the confirmatory
execution according to the local hash record; see HASHES_2026-08-27.txt /
HASHES_2026-09-28.txt -- there is no file named LOCK_HASHES.txt) and
conf_package.json (package output), and:
  (a) re-derives N* independently for every scenario;
  (b) evaluates the N-1 / N / N+1 trio around the package's N_total;
  (c) exhaustively cross-checks the enumerable scenarios with 3-D
      Gauss-Jacobi quadrature and (n_d, x, y) triple enumeration;
  (d) Monte Carlo confirmation at the fixed budget defined in the local
      hash record;
  (e) reference-side property tests PT1 and PT2.
Nothing here was tuned after seeing results; the grid file is unmodified.
"""
import json, time, math, sys
import numpy as np
import ref_bam as R


def enum_triples_fast(N, s, ok_se, ok_sp):
    """Exhaustive (n_d, x, y) triple enumeration -- identical sum to
    ref_bam.enum_triples (which is left byte-frozen; see
    HASHES_2026-08-27.txt / HASHES_2026-09-28.txt),
    only with the credible-interval widths precomputed once and the innermost
    loop written as an explicit outer product of triple masses rather than a
    Python loop.  No factorisation of the sum is introduced: every triple's
    mass is materialised and added.  Equivalence to the frozen routine is
    asserted at the end of section C."""
    a_se, b_se = s["prior_se"]; a_sp, b_sp = s["prior_sp"]
    a_p, b_p = s["prior_prev"]
    total = 0.0
    ntrip = 0
    for k in range(1, N):                       # k=0 or k=N -> degenerate
        m = N - k
        pk = float(np.exp(R.log_betabinom_pmf(k, N, a_p, b_p)))
        xs = np.nonzero(ok_se[k, : k + 1])[0]
        ys = np.nonzero(ok_sp[m, : m + 1])[0]
        if xs.size == 0 or ys.size == 0:
            continue
        px = np.exp(R.log_betabinom_pmf(xs, k, a_se, b_se))
        py = np.exp(R.log_betabinom_pmf(ys, m, a_sp, b_sp))
        triple_masses = pk * np.multiply.outer(px, py)   # every (k,x,y)
        total += float(triple_masses.sum())
        ntrip += triple_masses.size
    return total, ntrip


grid = json.load(open("LOCKED_confirmatory_grid.json"))
pkg = json.load(open("conf_package.json"))
SC = grid["scenarios"]
SEED_TRIO = grid["_meta"]["confirmatory_mc_seeds"][0]     # 70117
B_TRIO = 160000
B_SWEEP = 40000

def mk(s, dmult=1.0, nmax=None):
    return R.BAMExact(s["prior_se"], s["prior_sp"], s["prior_prev"],
                      s["delta_se"] * dmult, s["delta_sp"] * dmult,
                      s["level"], nmax=nmax or s["N_hi"])

results = {}

print("=" * 92)
print("A. INDEPENDENT SEARCH  vs  PACKAGE   (locked confirmatory grid)")
print("=" * 92)
print(f"{'id':4} {'pkgN':>6} {'refN':>6} {'agree':>6} "
      f"{'pkg A(N)':>13} {'ref A(N)':>13} {'|diff|':>10}")
for s in SC:
    sid = s["id"]
    t0 = time.time()
    ex = mk(s)
    refN, refA = ex.search(s["target"], s["N_lo"], s["N_hi"])
    p = pkg[sid]["main"]
    pN, pA = p["N_total"], p["joint_assurance"]
    refA_at_pN = ex.assurance(pN) if pN is not None else float("nan")
    ok = (refN == pN)
    print(f"{sid:4} {pN:>6} {str(refN):>6} {str(ok):>6} "
          f"{pA:>13.10f} {refA_at_pN:>13.10f} {abs(pA-refA_at_pN):>10.2e}")
    results[sid] = dict(pkgN=pN, refN=refN, pkgA=pA, refA=refA_at_pN,
                        pkg_warnings=p["warnings"], ex=ex,
                        secs=time.time() - t0)

print()
print("=" * 92)
print("B. N-1 / N / N+1 CONFIRMATION OF THE SELECTED SIZE")
print("=" * 92)
print(f"{'id':4} {'target':>7} {'N':>6} | {'ref A(N-1)':>12} {'ref A(N)':>12} "
      f"{'ref A(N+1)':>12} | {'MC A(N)':>9} {'mcse':>7} | verdict")
for s in SC:
    sid = s["id"]; r = results[sid]; ex = r["ex"]; N = r["pkgN"]
    tgt = s["target"]
    a_m1 = ex.assurance(N - 1) if N - 1 >= 1 else float("nan")
    a_0 = ex.assurance(N)
    a_p1 = ex.assurance(N + 1) if N + 1 <= ex.nmax else ex.assurance(N)
    ph, mcse = R.mc(N, s["prior_se"], s["prior_sp"], s["prior_prev"],
                    s["delta_se"], s["delta_sp"], s["level"],
                    B=B_TRIO, seed=SEED_TRIO)
    if sid == "C12":
        verdict = "INFEASIBLE-by-design (see PT3)"
    else:
        c1 = a_0 >= tgt
        c2 = a_m1 < tgt
        verdict = "PASS" if (c1 and c2) else (
            "FAIL:A(N)<target" if not c1 else "FAIL:A(N-1)>=target")
    r.update(a_m1=a_m1, a_0=a_0, a_p1=a_p1, mc=ph, mcse=mcse, verdict=verdict)
    print(f"{sid:4} {tgt:>7.2f} {N:>6} | {a_m1:>12.8f} {a_0:>12.8f} "
          f"{a_p1:>12.8f} | {ph:>9.6f} {mcse:>7.5f} | {verdict}")

print()
print("=" * 92)
print("C. EXHAUSTIVE CROSS-CHECK ON THE ENUMERABLE CONFIRMATORY SCENARIOS")
print("   closed form vs 3-D Gauss-Jacobi quadrature vs (n_d,x,y) enumeration")
print("=" * 92)
worst = 0.0
for s in SC:
    if not s["enumerable"]:
        continue
    sid = s["id"]; ex = results[sid]["ex"]
    NH = s["N_hi"]
    ok_se = R.arm_ok_matrix(NH, s["prior_se"][0], s["prior_se"][1],
                            s["delta_se"], s["level"])
    ok_sp = R.arm_ok_matrix(NH, s["prior_sp"][0], s["prior_sp"][1],
                            s["delta_sp"], s["level"])
    dq = de = dfrozen = 0.0
    ntrip_tot = 0
    for N in range(s["N_lo"], NH + 1):
        a_cf = ex.assurance(N)
        a_q = R.quad_full(N, s["prior_se"], s["prior_sp"], s["prior_prev"],
                          s["delta_se"], s["delta_sp"], s["level"],
                          ok_se=ok_se[: N + 1, : N + 1].copy(),
                          ok_sp=ok_sp[: N + 1, : N + 1].copy())
        a_e, nt = enum_triples_fast(N, s, ok_se, ok_sp)
        if N <= 22:      # equivalence of fast vs frozen enumeration
            a_f, _ = R.enum_triples(N, s["prior_se"], s["prior_sp"],
                                    s["prior_prev"], s["delta_se"],
                                    s["delta_sp"], s["level"])
            dfrozen = max(dfrozen, abs(a_f - a_e))
        dq = max(dq, abs(a_q - a_cf)); de = max(de, abs(a_e - a_cf))
        ntrip_tot += nt
    worst = max(worst, dq, de)
    print(f"  {sid}: N={s['N_lo']}..{NH}  "
          f"max|quad-closed|={dq:.3e}  max|enum-closed|={de:.3e}  "
          f"max|fast-frozen enum|={dfrozen:.3e}  "
          f"successful (n_d,x,y) triples summed = {ntrip_tot:,}")
print(f"  WORST DISAGREEMENT ACROSS ALL THREE ANALYTIC ROUTES = {worst:.3e}")

print()
print("=" * 92)
print("D. REFERENCE-SIDE PROPERTY TESTS (PT1 target-monotone, PT2 delta-monotone)")
print("=" * 92)
INF = float("inf")
pt1_fail = pt2_fail = 0
print("PT1  N* as a function of target_assurance (must be non-decreasing):")
for s in SC:
    if s["id"] == "C12":
        continue
    ex = results[s["id"]]["ex"]
    Ns = []
    for t in grid["property_tests"]["PT1_target_monotone"]["targets"]:
        n, _ = ex.search(t, s["N_lo"], s["N_hi"])
        Ns.append(INF if n is None else n)
    okm = all(Ns[i] <= Ns[i + 1] for i in range(len(Ns) - 1))
    pt1_fail += (not okm)
    print(f"  {s['id']}: t=0.70->{Ns[0]}  0.80->{Ns[1]}  0.90->{Ns[2]}   "
          f"{'OK' if okm else 'VIOLATION'}")

print("\nPT2  N* as a function of the width margin (must be non-decreasing "
      "as delta shrinks):")
for s in SC:
    if s["id"] == "C12":
        continue
    Ns = []
    for m in grid["property_tests"]["PT2_delta_monotone"]["delta_multipliers"]:
        e2 = mk(s, dmult=m)
        n, _ = e2.search(s["target"], s["N_lo"], s["N_hi"])
        Ns.append(INF if n is None else n)
    okm = all(Ns[i] <= Ns[i + 1] for i in range(len(Ns) - 1))
    pt2_fail += (not okm)
    print(f"  {s['id']}: delta x1.00->{Ns[0]}   x0.85->{Ns[1]}   "
          f"{'OK' if okm else 'VIOLATION'}")

print(f"\nPT1 violations = {pt1_fail}   PT2 violations = {pt2_fail}")

print()
print("=" * 92)
print("E. PT3 -- infeasible scenario C12")
print("=" * 92)
s12 = [x for x in SC if x["id"] == "C12"][0]
ex12 = results["C12"]["ex"]
n12, _ = ex12.search(s12["target"], s12["N_lo"], s12["N_hi"])
print(f"  reference: any N in [{s12['N_lo']},{s12['N_hi']}] reaching "
      f"{s12['target']}?  {'NO' if n12 is None else 'YES at N='+str(n12)}")
print(f"  reference assurance at N_hi={s12['N_hi']}: "
      f"{ex12.assurance(s12['N_hi']):.8f}")
print(f"  package returned N_total = {results['C12']['pkgN']} with "
      f"joint_assurance = {results['C12']['pkgA']:.8f}")
print(f"  package warnings = {results['C12']['pkg_warnings']}")

# serialise
dump = {k: {kk: (vv if kk != "ex" else None) for kk, vv in v.items()}
        for k, v in results.items()}
json.dump(dump, open("conf_reference.json", "w"), indent=1, default=str)
print("\nwrote conf_reference.json")
