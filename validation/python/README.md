# Python validation of `dtasamplesize::bam_sample_size(method = "exact")`

## Purpose

This directory holds a **separate reference implementation**
(`ref_bam.py`) of the joint sensitivity/specificity ("Se/Sp") assurance
calculation performed by `dtasamplesize::bam_sample_size(method = "exact")`,
plus development scripts (`dev_*`), confirmatory scripts (`conf_*`), a locked
scenario grid (`LOCKED_confirmatory_grid.json`), and the logs/JSON outputs of
one completed run. It checks the package's exact computation against a
from-scratch Python implementation; it does not re-derive or critique the
statistical specification itself (see Limitations).

## Provenance

`ref_bam.py` and the confirmatory scripts were written on **2026-08-27** by an
LLM-based coding agent (Claude, Anthropic) under the authors' supervision,
working **only from the package help page** (`man/bam_sample_size.Rd`) and
`NAMESPACE` — it did **not** read the R source under `R/`. The authors
reviewed the results. This separation is the point: agreement between a
description-only reimplementation and the package's own source is more
informative than agreement between two readings of the same source code.

## Reconstruction (2026-09-16)

The original working directory was deleted; the files here were reconstructed
byte-for-byte on **2026-09-16** from the session record. `ref_bam.py` and
`LOCKED_confirmatory_grid.json` match the SHA-256 hashes recorded at lock
time; `conf_PY_02`, `dev_01`, `dev_02`, `conf_PY_06`, and the R scripts match
the hashes recorded at 18:56Z (see `HASHES_2026-08-27.txt`), and so does
`conf_PY_05_shared_assumptions.py` once a one-line shell correction made at
18:32Z (HPD bisection direction) is re-applied; a first reconstruction that
missed it is kept, labelled SUPERSEDED, in `logs/`.
`conf_PY_08_degeneracy_bite.py`/`conf_PY_09_published_mc.py` were written
later and have no recorded hash.

## Correction (2026-09-28, H-03)

`conf_PY_06_monotone_in_N.py` was edited on 2026-09-28 to declare a tolerance
(`TOL = 1e-9`) for counting a step as a real reversal, clip per-arm P (and the
joint assurance curve) to `[0, 1]` before differencing, and print the
tolerance and m-range used. This intentionally breaks its 2026-08-27 locked
hash (`57fa5ec9f...`); see `HASHES_2026-09-28.txt` for the reason and the new
hash (`5e6a6975...`). The re-run shows **73** real reversals for the
published P1/P2 Se prior over m in 2..800 (previously reported as 395, which
counted floating-point roundoff near `P == 1` as reversals); joint assurance
`A(N)` remains non-decreasing for all 13 scenarios. `ref_bam.py` and
`LOCKED_confirmatory_grid.json` were not touched and keep their original
hashes.

## Reverification (2026-09-28, H-05)

`conf_R_01_package.R` and `conf_R_07_reverify_version.R` were re-run against
**dtasamplesize 0.6.5** (built and installed from this release's tarball into
a private library, not the system library). All 12 confirmatory scenarios
(C01-C12) and both published scenarios (P1: N=678, A=0.8003489948; P2:
N=672, A=0.8002692084) reproduce exactly (`ALL CONFIRMATORY RESULTS UNCHANGED
ACROSS THE REBUILD: TRUE`) -- expected, since H-06 (the only 0.6.5 code
change) touches `joint_sample_size()`, not `bam_sample_size()`, which is all
these scripts exercise. Logs archived to
`logs/conf_R_01_package.R.log` and `logs/conf_R_07_reverify_version.R.log`.
The rest of the confirmatory/development grid (the Python scripts) was not
re-run under 0.6.5: their computational surface is unaffected by H-06, and
the 2026-09-16 logs already reflect the H-03 correction above.

**What the lock shows and does not show.** The hashes were recorded locally
in the agent's own session log, not with any external timestamping service or
registry. They show the grid and reference implementation were unchanged
between locking and running (internal consistency), not that they were
registered with a third party in advance. This is correctly described as
**pre-specified and hash-locked**, not **pre-registered**.

## Design

`dev_*` scripts probe 5 development scenarios (D1-D5), chosen freely while
building the reference implementation. `conf_*` scripts run only against the
locked confirmatory grid: 12 scenarios (C01-C12), 3 Monte Carlo seeds, and 6
property tests (PT1-PT6). Development and confirmatory scenarios use disjoint
parameter sets by construction.

## How to run

```bash
pip install -r requirements.txt   # numpy, scipy (R/dtasamplesize installed separately)
./run_all.sh    # POSIX, or ./run_all.ps1 for PowerShell
```

The rerun documented here (`run_all_top.log`, `logs/*.log`, 2026-09-16) used
**dtasamplesize 0.6.3**, whose computational code is identical to 0.6.4
(which changed only documentation/tests/validation tooling) and, for
`bam_sample_size()`, to 0.6.5 (which changed only `joint_sample_size()`; see
"Reverification (2026-09-28, H-05)" above for the 0.6.5 re-run of the two
R-side confirmatory scripts).

## Results

Manuscript claims checked against today's (2026-09-16) logs/JSON. Only MATCH
rows carry reproduced values.

| Claim | Verdict | Today's value | Source |
|---|---|---|---|
| Analytic routes agree within 4.086e-14 | MATCH | 4.086e-14 | `logs/conf_PY_02_reference.py.log:42` |
| Reference vs package agree within 7.20e-13 (max diff) | MATCH | 7.20e-13 (C06) | `logs/conf_PY_02_reference.py.log:10` |
| 12 scenarios; 11 feasible (C12 infeasible) | MATCH | C01-C12; C12 infeasible-by-design | `logs/conf_PY_02_reference.py.log:5-16,33,76-81` |
| N-1/N/N+1 target reached at N not N-1, all feasible scenarios | MATCH | PASS on all 11 | `logs/conf_PY_02_reference.py.log:22-32` |
| 41 property-based tests total | MATCH | 19 (R) + 22 (Python PT1 x11 + PT2 x11) = 41 | `logs/conf_R_03_properties.R.log:28`; `logs/conf_PY_02_reference.py.log:47-58,60-71,73` |
| 1,846,401 enumerated triplets | MATCH | 604,173+631,595+610,633=1,846,401 | `logs/conf_PY_02_reference.py.log:39-41` |
| N=678 (0.8003489948) and N=672 (0.8002692084) | MATCH | reproduced exactly | `logs/conf_R_07_reverify_version.R.log:4-5` |
| Half-width misreading -> assurance 0.0032 at N=678 | MATCH | 0.003178 | `logs/conf_PY_05_shared_assumptions.py.log:7` |
| Fixed diseased count -> N=533 | MATCH | N=533 (P1 and P2) | `logs/conf_PY_05_shared_assumptions.py.log:6,13` |
| HPD intervals -> 666 and 658 | MATCH | 666 (P1), 658 (P2) | `logs/conf_PY_05_shared_assumptions.py.log:4,11` |
| 73 real non-monotone per-arm reversals for m<=800 (TOL=1e-9, P clipped to [0,1]); joint assurance non-decreasing | MATCH | 73 (P1/P2 P_Se, m 2..800); A(N) non-decreasing, 0 dips, 13/13 rows | `logs/conf_PY_06_monotone_in_N.py.log:1-15,24,26` |
| MC budgets 40,000 and 160,000; combined total under 1e7 | MATCH | both budgets used; no printed grand total, but B x calls sums to ~6-8M | `conf_PY_02_reference.py:50-51`; `conf_R_03_properties.R:68-70`; `conf_R_04_mc_vs_exact.R:10` |
| RNG state unchanged after a call | MATCH | PT5 PASS, seed and `runif(3)` unaffected | `logs/conf_R_03_properties.R.log:20-22` |
| Grid = 12 scenarios x 3 seeds x 6 property tests; dev = 5, disjoint | MATCH | ids/seeds/PT1-PT6 confirmed; D1-D5 distinct priors | `LOCKED_confirmatory_grid.json:16,22-33`; `dev_01_agreement_smallN.py:13-17` |

**vs. `report.html` (2026-08-27):** all rows above agree exactly with the
historical report. `conf_PY_05` was first rerun from an incomplete
reconstruction (without the 18:32Z correction) and returned no HPD sample
size; that log is kept as
`logs/conf_PY_05_shared_assumptions.py.SUPERSEDED_incomplete_reconstruction.log`.
The current log comes from the file whose hash matches the 18:56Z record.

## Limitations

- **Agreement between two implementations cannot validate the specification
  itself.** Both `ref_bam.py` and the package share the same reading of
  ambiguous points (equal-tailed vs. HPD interval, random vs. fixed diseased
  count, delta as full vs. half width). `conf_PY_05_shared_assumptions.py`
  quantifies how far the answer moves under alternative readings; it does not
  prove the shared reading is correct.
- **Correlated implementers.** Both implementations were built with LLM
  assistance from the same model family, so agreement does not exclude a
  correlated misreading of the specification.
- **Scope.** Covers `method = "exact"` only, on 12 confirmatory + 2 published
  scenarios; not exhaustive over the package's parameter space.

## What is kept

Every script, the locked grid, the JSON outputs and the per-script logs of the
2026-09-16 run are retained as evidence (`logs/`, `run_all_top.log`,
`conf_*.json`, `dev_01_agreement_smallN.json`). `__pycache__/` is excluded via
`.gitignore`.
