"""
DEVELOPMENT step 2: reproduce the two PUBLISHED values with the independent
reference implementation.

P1: prior_se=(17,3) prior_sp=(2,2)  prior_prev=(4,16) d_se=.14 d_sp=.10 t=.80
    package claims N = 678, assurance = 0.800349, A(677) = 0.799685
P2: prior_sp=(18,2), rest identical
    package claims N = 672, assurance = 0.800269
"""
import time
import numpy as np
import ref_bam as R

LEVEL = 0.95
TARGET = 0.80
NMAX = 900

CASES = [
    ("P1", (17, 3), (2, 2), (4, 16), 0.14, 0.10, 678, 0.800349, 0.799685),
    ("P2", (17, 3), (18, 2), (4, 16), 0.14, 0.10, 672, 0.800269, None),
]

for tag, pse, psp, pprev, dse, dsp, N_claim, A_claim, Am1_claim in CASES:
    t0 = time.time()
    ex = R.BAMExact(pse, psp, pprev, dse, dsp, LEVEL, nmax=NMAX)
    Nfound, Afound = ex.search(TARGET, N_lo=2, N_hi=NMAX)
    print(f"\n=== {tag} ===  (cache built + full scan in {time.time()-t0:.1f}s)")
    print(f"  reference N*      = {Nfound}      package claim = {N_claim}")
    print(f"  reference A(N*)   = {Afound:.10f}  package claim = {A_claim}")
    for d in (-2, -1, 0, 1, 2):
        N = N_claim + d
        print(f"     A({N}) = {ex.assurance(N):.10f}"
              + ("   <-- claimed N" if d == 0 else "")
              + (f"   [package claim {Am1_claim}]"
                 if (d == -1 and Am1_claim is not None) else ""))
    # monotonicity of the assurance curve near the crossing
    Ns = np.arange(N_claim - 30, N_claim + 31)
    A = ex.curve(Ns)
    print(f"  assurance strictly increasing on [{Ns[0]},{Ns[-1]}]? "
          f"{bool(np.all(np.diff(A) > 0))}")
    print(f"  first N in that window with A>=target: "
          f"{int(Ns[np.argmax(A >= TARGET)])}")
