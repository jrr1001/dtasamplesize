# Historical artefacts (not current validation evidence)

The files in this folder are kept for provenance only. They are **not**
evidence for the 0.6.8 release and must not be cited as such.

- `HASHES_2026-10-01.txt` -- an earlier SHA-256 manifest of
  `validation/`, superseded the same day by `HASHES_2026-10-01_v2.txt`
  (which differs only in the hashes of `python/report.html` and
  `reproduce_manuscript.R`, both regenerated after the v1 manifest was
  recorded). Superseded again by `validation/HASHES_2026-10-05_v0.6.8.txt`
  for this release.
- `exact_invariance_v065.log` -- the deliberate demonstration run of
  `validation/exact_invariance.R --version-check=off` against the
  (pre-fix) 0.6.5 installation, kept to show the G04 gate correctly
  FAILS against the version with the defect it exists to catch. Not a
  current-version check.

See `validation/README.md` for what each current script proves.
