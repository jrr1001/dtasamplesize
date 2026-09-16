"""
MCSE-budgeted Monte Carlo confirmation of the two PUBLISHED decision points.
B = 160000 -> MCSE ~ 0.0010 at p = 0.80 (pre-registered budget; no 1e8 runs).
Three independent seeds from the locked grid.
The simulation draws the generative model directly (prev, Se, Sp, n_d, x, y)
and uses NO closed form, so it is an independent check of the analytic route.
"""
import numpy as np
import ref_bam as R

B = 160000
SEEDS = [70117, 70118, 70119]
CASES = [("P1", (17, 3), (2, 2), 678, 0.8003489948, 677, 0.7996848824),
         ("P2", (17, 3), (18, 2), 672, 0.8002692084, 671, 0.7996594780)]

print(f"{'case':5} {'N':>5} {'seed':>7} {'exact':>13} {'MC':>11} "
      f"{'mcse':>8} {'z':>7}")
for tag, pse, psp, N, Aex, Nm1, Am1 in CASES:
    for Nq, Aq, lab in ((Nm1, Am1, "N-1"), (N, Aex, "N  "), (N + 1, None, "N+1")):
        ex = R.BAMExact(pse, psp, (4, 16), 0.14, 0.10, 0.95, nmax=Nq)
        a = ex.assurance(Nq)
        if Aq is not None:
            assert abs(a - Aq) < 1e-9, (Nq, a, Aq)
        for sd in SEEDS:
            ph, mcse = R.mc(Nq, pse, psp, (4, 16), 0.14, 0.10, 0.95,
                            B=B, seed=sd)
            z = (ph - a) / mcse
            print(f"{tag:5} {Nq:>5} {sd:>7} {a:>13.10f} {ph:>11.6f} "
                  f"{mcse:>8.5f} {z:>+7.2f}   [{lab}]")
    print()
print("All |z| below 3 means the closed form and a direct simulation of the "
      "generative model are indistinguishable at MCSE ~ 0.001.")
