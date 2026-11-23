# Historical artefacts (not current validation evidence)

The files in this folder are kept for provenance only. They are **not**
evidence for the 0.7.0 release and must not be cited as such.

- `HASHES_2026-10-01.txt` -- an earlier SHA-256 manifest of
  `validation/`, superseded the same day by `HASHES_2026-10-01_v2.txt`
  (which differs only in the hashes of `python/report.html` and
  `reproduce_manuscript.R`, both regenerated after the v1 manifest was
  recorded). Superseded by `HASHES_2026-10-05_v0.6.8.txt` (0.6.8), itself
  superseded by `validation/HASHES_2026-10-06_v0.6.9.txt` (0.6.9), which is in
  turn superseded by `validation/HASHES_2026-11-23_v0.7.0.txt` (0.7.0).
- `HASHES_2026-10-05_v0.6.8.txt` -- the SHA-256 manifest of `validation/`
  for the 0.6.8 release, superseded by
  `validation/HASHES_2026-10-06_v0.6.9.txt` (0.6.9 is a maintenance
  release: no computational code, API, or results changed; most file
  hashes are unchanged, but the manifest itself is recomputed per
  release).
- `exact_invariance_v065.log` -- the deliberate demonstration run of
  `validation/exact_invariance.R --version-check=off` against the
  (pre-fix) 0.6.5 installation, kept to show the G04 gate correctly
  FAILS against the version with the defect it exists to catch. Not a
  current-version check.

- `R_CMD_check_with_pdf_manual_00check_v0.6.9.log` -- the
  `R CMD check --as-cran` of 0.6.9 that also attempted the PDF manual (1 ERROR,
  1 WARNING: the local TeX installation lacks the Courier font metrics). It cannot
  be regenerated here (`makeindex` is missing) and is not evidence for 0.7.0.
- `HASHES_2026-10-06_v0.6.9.txt` (in `validation/`) is the manifest of 0.6.9;
  the current one is `validation/HASHES_2026-11-23_v0.7.0.txt`.

See `validation/README.md` for what each current script proves.
