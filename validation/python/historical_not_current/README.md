# Historical artefacts (not current validation evidence)

The files in this folder document the provenance of the Python
reference-implementation validation but are **not** current evidence for
the 0.6.8 release. Do not cite them as confirming any current-release
number.

- `HASHES_2026-08-27.txt` -- the original lock-time hash record for
  `ref_bam.py`, `LOCKED_confirmatory_grid.json`, and the confirmatory/
  development scripts as they stood on 2026-08-27, before the 2026-09-16
  byte-for-byte reconstruction and the 2026-09-28 (H-03) correction to
  `conf_PY_06_monotone_in_N.py`. `validation/python/README.md`'s
  "Provenance" and "Reconstruction" sections quote the hashes recorded
  here; the grid and reference implementation they describe are
  unchanged today (verified against `LOCKED_confirmatory_grid.json`'s
  current SHA-256 for this release).
- `HASHES_2026-09-28.txt` -- the hash record for the H-03 correction to
  `conf_PY_06_monotone_in_N.py` (declaring `TOL = 1e-9`, clipping P to
  `[0, 1]`), which intentionally broke the 2026-08-27 hash of that one
  file; see `validation/python/README.md`'s "Correction (2026-09-28,
  H-03)" section.
- `conf_PY_05_shared_assumptions.py.SUPERSEDED_incomplete_reconstruction.log`
  -- a first, incomplete reconstruction of this script's log (missing the
  18:32Z HPD-bisection-direction correction), which returned no HPD
  sample size. Superseded by the log matching the 18:56Z hash record,
  kept under `validation/python/logs/`.

See `validation/python/README.md` for the full narrative and
`validation/HASHES_2026-10-05_v0.6.8.txt` for this release's current
manifest.
