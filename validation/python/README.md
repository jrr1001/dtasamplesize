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
LLM-based coding assistant under the authors' supervision,
working **only from the package help page** (`man/bam_sample_size.Rd`) and
`NAMESPACE` — it did **not** read the R source under `R/`. The authors
reviewed the results. This separation is the point: agreement between a
description-only reimplementation and the package's own source is more
informative than agreement between two readings of the same source code.
It is nevertheless a same-team, LLM-assisted check, not external software
and not an independent implementation in the sense of an unrelated group
(see Limitations).

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
N=672, A=0.8002692084) reproduce exactly, as printed by the two scripts;
this is expected, since H-06 (the only 0.6.5 code change) touches
`joint_sample_size()`, not `bam_sample_size()`, which is all these scripts
exercise. The line `ALL CONFIRMATORY RESULTS UNCHANGED ACROSS THE REBUILD:
TRUE` that `conf_R_07` prints is not evidence of unchanged results across a
rebuild; see "What the line `ALL CONFIRMATORY RESULTS UNCHANGED ACROSS THE
REBUILD` does and does not show" below. Logs archived to
`logs/conf_R_01_package.R.log` and `logs/conf_R_07_reverify_version.R.log`.
The rest of the confirmatory/development grid (the Python scripts) was not
re-run under 0.6.5: their computational surface is unaffected by H-06, and
the 2026-09-16 logs already reflect the H-03 correction above.

**What the lock shows and does not show.** The hashes were recorded locally
in the assistant's own session log (`HASHES_2026-08-27.txt`, `HASHES_2026-09-28.txt`
-- the actual hash-record files kept in this directory; there is no file
named `LOCK_HASHES.txt`), not with any external timestamping service or
registry. They show the grid and reference implementation were unchanged
between recording the hash and running the confirmatory scripts (internal
consistency), not that they were registered with a third party in advance.
The confirmatory grid is correctly described as **defined before the
confirmatory execution according to the local hash record**, never as
"pre-registered" and never cited as proof of preespecification by itself.

## Regeneration under 0.7.0 (current)

The current evidence is the full `run_all.sh` pipeline (dev_01-03, conf_R_01,
conf_PY_02, conf_R_03, conf_R_04, conf_PY_05, conf_PY_06, conf_R_07,
conf_PY_08, conf_PY_09) plus the top-level `validation/*.R` scripts, re-run on
2026-11-23 against **dtasamplesize 0.7.0** installed in a private library (not
the system library; `R_LIBS_PRIVATE` pointed at it), with R 4.6.2. The registered
scripts ran without edits (the SHA-256 of the 14 registered files is the one
recorded when 0.7.0 was released and was checked before and after the run). The
logs of the R scripts that print the package version (conf_R_01, conf_R_07,
dev_03, compare_confirmatory_baselines and sep_priors_check) state it; the Python
scripts read the JSON written by conf_R_01 or do not call the package, and
conf_R_03 and conf_R_04 do not print it (the installed library is recorded in
`../logs/sessionInfo.txt`). The pipeline ended with exit
code 0 (`run_all_top.log`), and `compare_confirmatory_baselines.R` was then run
by hand as an additional step (its log is `logs/compare_confirmatory_baselines.R.log`).

Compared with the 0.6.9 logs, after setting aside versions, dates, times and
paths, the only differences are the following. In `conf_R_01_package.R.log`,
scenario C10 reports `warns=2` instead of `warns=0`: version 0.7.0 warns that
the analysis prior (equal to the design prior) already has a prior-only interval
width within the target; `N_total` and `joint_assurance` are unchanged. In
`dev_03_probe_pkg.R.log`, the list of fields of the object returned by
`bam_sample_size()` has the new fields of 0.7.0 (`alpha_ci`, `width_type`,
`assurance_for` and the analysis-prior fields). Every other number in the twelve
logs is identical to 0.6.9, including the N, the joint assurances, the
agreement of 7.20e-13 and the 41 property checks.

`conf_package.json`, `conf_reference.json`, `conf_properties.json` and
`dev_01_agreement_smallN.json` were regenerated; they are outputs, not
registered files. `report.html` was regenerated from these logs by
`generate_report.R`.

### Scripts added after the registered ones

- `compare_confirmatory_baselines.R` compares C01-C12 of the regenerated
  `conf_package.json` with two earlier files kept in
  `historical_not_current/`: `conf_package_38f12a8_v0.6.9.json` (the file in
  commit `38f12a8`, 0.6.9) and
  `conf_package_bc943e0_rerun_2026-09-16_v0.6.3.json` (the file in commit
  `bc943e0`, the re-run of 2026-09-16 with 0.6.3). Both were extracted byte for byte
  with `git show`. For every field the three files have in common it prints
  `identical()` and the absolute difference, and for every difference both
  values; it describes the differences and does not judge them. In this run, the
  regenerated file and the 0.6.9 file differ only in the new warnings of C10;
  the regenerated file and the 0.6.3 re-run differ in the assurance of the
  single-N probe at N-1 (`NA` in the regenerated file, a number in the 0.6.3
  re-run; this is the single-N behaviour of the anti-noise rule described
  under `conf_R_04_mc_vs_exact.R` below, which is not examined further here), in the same C10 warnings and in C12 (N = 400, A = 0 in the 0.6.3
  re-run; `NA` and `target_reached = FALSE` now). All `N_total`, `joint_assurance`
  and `mcse` values of C01-C11 are identical in the three files. The script is
  run as the last step of `run_all.sh` and `run_all.ps1`.
- `sep_priors_check.py` (with `ref_bam_sep_priors.py`) compares
  `bam_sample_size(method = "exact")` of 0.7.0 with a Python reimplementation
  for separate design and analysis priors, `width_type` and `assurance_for`. It
  was written by an LLM-based coding assistant under the authors' supervision,
  working from the help pages, and is a check by the same team. It is run as
  `python -I sep_priors_check.py <out>` with `RSCRIPT` set to the Rscript
  executable; its console log is `logs/sep_priors_check.py.log` and its own log
  and CSV are in `../logs/sep_priors_check/`. Result of the run of 2026-11-23: 290
  comparisons, 0 failures at the tolerance of 1e-10 in A(N) with identical N;
  largest difference 1.215e-12.
- `wilson2022_rule_enumeration.py` is described in `../README.md`.

## Regeneration under 0.6.9 (historical)

The "Reverification (2026-09-28, H-05)" section and the "Results
(2026-09-16, historical)" section describe earlier runs against
**dtasamplesize 0.6.5 / 0.6.3**; they are kept, explicitly labelled by
version, as the historical record of those runs and are not current validation
evidence. The run of 2026-10-06 was the full
`run_all.sh`/`run_all.ps1` pipeline (dev_01-03, conf_R_01, conf_PY_02,
conf_R_03, conf_R_04, conf_PY_05, conf_PY_06, conf_R_07, conf_PY_08,
conf_PY_09) plus the top-level `validation/*.R` scripts, re-run against a
privately installed **dtasamplesize 0.6.9** (not the system library;
`R_LIBS_PRIVATE` pointed at a private library) on 2026-10-06, exit code 0
end to end. `validation/python/generate_report.R` (the public, versioned
replacement for the undocumented machine-local
`tools/generar_report_html.R` referenced by earlier revisions of this
file, never committed to this repository) regenerated `report.html` from
those logs. The "Results" table below is the 0.6.9 confirmation.

This run also exercised the M-03 fix (see `NEWS.md`, "dtasamplesize
0.6.8"): `conf_R_01_package.R` no longer crashes probing `N_total +/- 1`
around the infeasible scenario C12, and reports it explicitly
(`target_reached = FALSE`, `N_total = NA`, `joint_assurance = NA`,
`max_assurance_evaluated = 0` at `N = 10`); `conf_PY_02_reference.py` and
`conf_R_03_properties.R`'s PT3 agree C12 is infeasible rather than
crashing on a `None`/`NA` `N_total`; `conf_R_07_reverify_version.R` prints
`TRUE` for its comparison of all 12 scenarios (C01-C12, including C12) with
`conf_package.json`, but that comparison is of a regenerated output with
itself (see the next section).

The logs of that run were replaced by those of 0.7.0 (see above); they remain in
the git history. Version 0.7.0 returns identical results for the default
arguments (`NEWS.md`, "Compatibility with 0.6.9"; `../logs/compat_v069_check/`).

## What the line `ALL CONFIRMATORY RESULTS UNCHANGED ACROSS THE REBUILD` does and does not show

`conf_R_07_reverify_version.R` reads `conf_package.json`, calls the installed
package for each scenario and compares `N_total` (with `identical`) and
`joint_assurance` (tolerance 1e-12) with the JSON. In the `run_all` pipeline
`conf_R_01_package.R` runs first and rewrites `conf_package.json` with the
same installed package and the same arguments, so inside the pipeline the
comparison is between an output and its own immediate regeneration. It does
not compare with the run of 2026-08-27, whose `conf_package.json` and package
build are not preserved (the oldest `conf_package.json` in the repository
history is the output of the 2026-09-16 re-run with 0.6.3, and it was
overwritten in later regenerations). The line therefore does not show that
results are unchanged since the original confirmatory run.

The comparison that the line cannot make is made by
`compare_confirmatory_baselines.R` (see "Regeneration under 0.7.0"), which
compares the regenerated file with the two earlier files kept in the git
history; this document does not cite the line as evidence.

What was established, in a separate re-execution carried out for the 2026
revision of the article: the code of `conf_R_01_package.R` run with three
package builds (0.6.2, 0.6.3 and 0.6.9) gives identical `N_total` and
`joint_assurance` in scenarios C01-C11 (maximum difference 0 in
`joint_assurance`). Scenario C12 (infeasible by design) is `N_total = 400`,
`joint_assurance = 0` under the 0.6.2 and 0.6.3 convention and `NA` under
0.6.9; no published number depends on it. This does not show that the build
used on 2026-08-27 gave the same values, because that build cannot be
recovered. The header of `conf_R_07_reverify_version.R` (lines 5-6) says the
installed package is "a privately installed dtasamplesize 0.6.9"; that
mention is out of date for later releases. The script is a registered file
and is not edited; it prints the version actually installed
(`dtasamplesize now installed:`).

## Design

`dev_*` scripts probe 5 development scenarios (D1-D5), chosen freely while
building the reference implementation. `conf_*` scripts run only against the
locked confirmatory grid: 12 scenarios (C01-C12), 3 Monte Carlo seeds, and 6
property tests (PT1-PT6). Development and confirmatory scenarios use disjoint
parameter sets by construction.

## How to run

```bash
pip install -r requirements.txt   # numpy, scipy (R/dtasamplesize installed separately)
R_LIBS_PRIVATE="<library with dtasamplesize 0.7.0>" ./run_all.sh    # POSIX, or ./run_all.ps1 for PowerShell
```

`run_all.sh` ends with `compare_confirmatory_baselines.R`. The registered Python
scripts import `ref_bam` from their own folder, so they are run without the `-I`
flag (`python script.py`, as `run_all.sh` does).

The rerun of 2026-09-16 (no longer in the working tree) used
**dtasamplesize 0.6.3**, whose computational code is identical to 0.6.4
(which changed only documentation/tests/validation tooling) and, for
`bam_sample_size()`, to 0.6.5 (which changed only `joint_sample_size()`; see
"Reverification (2026-09-28, H-05)" above for the 0.6.5 re-run of the two
R-side confirmatory scripts).

## Results (0.7.0, current)

Manuscript claims checked against the 2026-11-23 logs/JSON produced by the
0.7.0 pipeline run described above (the 0.6.9 run of 2026-10-06 gave the same
values). Only MATCH rows carry reproduced values.

| Claim | Verdict | 0.7.0 value | Source |
|---|---|---|---|
| Analytic routes agree within 4.086e-14 | MATCH | 4.086e-14 | `logs/conf_PY_02_reference.py.log:40-42` |
| Reference vs package agree within 7.20e-13 (max diff) | MATCH | 7.20e-13 (C06) | `logs/conf_PY_02_reference.py.log` section A |
| 12 scenarios; 11 feasible (C12 infeasible) | MATCH | C01-C12; C12 infeasible-by-design (`target_reached = FALSE`, `N_total`/`joint_assurance` = NA) | `logs/conf_PY_02_reference.py.log` sections A, B, E; `logs/conf_R_01_package.R.log` |
| N-1/N/N+1 target reached at N not N-1, all feasible scenarios | MATCH | PASS on all 11; C12 "INFEASIBLE-by-design" | `logs/conf_PY_02_reference.py.log` section B |
| 41 property-based tests total | MATCH | 19/19 R + 22/22 Python (PT1 x11 + PT2 x11) = 41, 0 violations | `logs/conf_R_03_properties.R.log` ("19/19 property checks passed"); `logs/conf_PY_02_reference.py.log` ("PT1 violations = 0   PT2 violations = 0") |
| 1,846,401 enumerated triplets | MATCH | 604,173+631,595+610,633=1,846,401 (C09/C10/C11) | `logs/conf_PY_02_reference.py.log` section C |
| N=678 (0.8003489948) and N=672 (0.8002692084) | MATCH | reproduced exactly by the installed 0.7.0 (the C01-C12 `TRUE` line of `conf_R_07` compares a regenerated output with itself; see the section above) | `logs/conf_R_07_reverify_version.R.log` |
| A(677) = 0.7996848824 | MATCH | 0.7996848824 | `logs/reproduce_manuscript.log` |
| Negative control: treating the half-width delta/2 as the total-width target (a criterion twice as strict as the published one) -> assurance 0.0032 at N=678 | MATCH | 0.003178 (P1) | `logs/conf_PY_05_shared_assumptions.py.log` (P1, V3) |
| Fixed diseased count -> N=533 | MATCH | N=533 (P1 and P2) | `logs/conf_PY_05_shared_assumptions.py.log` (V5) |
| HPD intervals -> 666 and 658 | MATCH | 666 (P1), 658 (P2) | `logs/conf_PY_05_shared_assumptions.py.log` (V2) |
| 73 real non-monotone per-arm reversals for m<=800 (TOL=1e-9, P clipped to [0,1]); joint assurance non-decreasing | MATCH | 73 (P1/P2 P_Se) | `logs/conf_PY_06_monotone_in_N.py.log` |
| RNG state unchanged after a call | MATCH | PT5 PASS | `logs/conf_R_03_properties.R.log` |
| joint_sample_size N=580, p=0.8041, mcse=0.002806 | MATCH | reproduced exactly | `logs/reproduce_manuscript.log` |
| Buderer P(CI width <= target) = 0.5730 | MATCH | 0.5730 | `logs/reproduce_manuscript.log` |
| Wald/Wilson sweep ranges 0.469-0.872 and 0.000-1.000 | MATCH | Wald 0.4688-0.8721; Wilson 0.0000-1.0000 | `logs/sweep_wald_wilson.log` |
| G04 exact-mode invariance (0.7.0, expected PASS) | MATCH | PASS across 80 B/seed/N_range combinations + no-crossing contract | `logs/exact_invariance.log` |
| RNG-kind invariance (5 kinds) | MATCH | PASS | `logs/rng_invariance.log` |
| Integration check (full pipeline, reproducibility, ordering) | MATCH | PASS | `logs/integration_check.log` |
| Cross-validation (V1-V8: Wilson/Wald/Buderer/Hanley-McNeil/net benefit/Beta posterior/imperfect-ref) | MATCH | all PASS (2 documented expected deviations under V1/V4, not held to tolerance) | `logs/cross_validation.log` |

**Direction of the negative control.** In `conf_PY_05_shared_assumptions.py`
the variant V3 is labelled "delta read as HALF width" (and "delta treated as
HALF width" in its log), but what it computes is `width <= delta/2`: the
half-width delta/2 applied as if it were the total width, which is twice as
strict as the published criterion (`width <= delta`, delta being the full
width). That is the reading that gives 0.0032 at N = 678 for P1. The
opposite misreading, taking 2 x delta as the total width, gives an assurance
of 0.9998 at the same N (`validation/negative_control_halfwidth.R`, run
against `bam_sample_size(method = "exact")` at a single N). The script label
is registered and is not edited; this paragraph is the correction of its
direction.

**Wording of the registered files.** `ref_bam.py` (line 28),
`conf_PY_02_reference.py` (lines 2, 8 and 65, printed in its log as "A.
INDEPENDENT SEARCH vs PACKAGE"), `conf_PY_09_published_mc.py` (lines 5 and 7),
`dev_01_agreement_smallN.py` and `dev_02_published.py` (line 2) and
`conf_R_04_mc_vs_exact.R` (line 5) use "independent" for separate code paths
or seeds written by the same team. They are registered files and are not
edited. None of them is independent validation (see Limitations).

**MC-vs-exact internal consistency (`conf_R_04_mc_vs_exact.R`), NOT a
discrepancy in the manuscript numbers above.** This script probes
`bam_sample_size(method = "monte_carlo")` at the single, fixed published
N (678, 672) for 3 seeds each and compares the returned
`joint_assurance` against the exact value. Under the current API
(the "anti-noise" acceptance rule added in 0.6.6 -- see `NEWS.md`), a
single-N Monte Carlo probe only returns a finite `joint_assurance` when
the lower bound of its own 95% Monte Carlo CI still clears
`target_assurance`; otherwise it returns `NA` by design (never a noisy
point estimate passed off as a validated one). Since the published
`joint_assurance` (~0.8003) sits well within Monte Carlo noise of the
0.80 target at `B = 160000` (`mcse ~= 0.00100`), most single-N probes are
expected to land on the "not accepted" side and return `NA`. Of the 6
case/seed combinations, only one (P1, seed 70117) cleared the bound:
`joint_assurance = 0.8020937500`, `mcse = 0.00100`,
`z = (MC - exact)/mcse = +1.75`. The other 5 returned `NA` and no `z`
could be computed for them. **The only MCSE multiple this run
measures/supports is `|z| = 1.75` (one data point); it does not support
a stronger claim such as "compatible within 3x MCSE" across multiple
replications, because the current API does not expose a raw simulated
`joint_assurance` at a single N independent of the accept/reject
decision.** The git-tracked `logs/conf_R_04_mc_vs_exact.R.log` prior to
this release (last touched in the 0.6.4-era commit, before the 0.6.6
anti-noise fix) showed finite values for all 6 rows; that was a stale
log that predated the current, deliberate design and was never
regenerated across the 0.6.5/0.6.6/0.6.7 releases. This is now
corrected (the log is regenerated here), but the anti-noise design
itself means this particular script cannot usefully report "Monte Carlo
compatible within its MCSE" the way earlier, stale evidence implied;
anyone relying on that claim should re-derive it from a probe that
requests the raw simulated assurance rather than the accept/reject
headline.

## Results (2026-09-16, historical)

Manuscript claims checked against the 2026-09-16 logs/JSON, under
**dtasamplesize 0.6.3** (see "How to run" below). Kept for provenance;
not current evidence -- see "Results (0.7.0, current)" above. Only MATCH
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
| Negative control: treating the half-width delta/2 as the total-width target (twice as strict as published) -> assurance 0.0032 at N=678 | MATCH | 0.003178 (P1) | `logs/conf_PY_05_shared_assumptions.py.log:7` |
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
  correlated misreading of the specification. The Python checks are a
  same-team cross-check, not an external reference implementation, and
  this directory should not be described as independent validation.
- **Scope.** Covers `method = "exact"` only, on 12 confirmatory + 2 published
  scenarios; not exhaustive over the package's parameter space.

## What is kept

Every script, the locked grid, the JSON outputs and the per-script logs of the
2026-11-23 run are retained as evidence (`logs/`, `run_all_top.log`,
`conf_*.json`, `dev_01_agreement_smallN.json`, `report.html`), together with the
two earlier `conf_package.json` files in `historical_not_current/`. `__pycache__/` is excluded via
`.gitignore`.
