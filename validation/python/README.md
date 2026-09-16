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
(which changed only documentation/tests/validation tooling).

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
| 395 non-monotone per-arm reversals <900 subjects; joint assurance non-decreasing | MATCH | 395 (P1/P2 P_Se); A(N) non-decreasing, 0 dips, 13/13 rows | `logs/conf_PY_06_monotone_in_N.py.log:1-14,23,25` |
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
