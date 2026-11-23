# Historical artefacts (not current validation evidence)

The files in this folder document the provenance of the Python
reference-implementation validation but are **not** current evidence for
the 0.7.0 release. Do not cite them as confirming any current-release
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

- `conf_package_38f12a8_v0.6.9.json` and
  `conf_package_bc943e0_rerun_2026-09-16_v0.6.3.json` -- the `conf_package.json`
  of commit `38f12a8` (0.6.9) and of commit `bc943e0` (the re-run of 2026-09-16
  with 0.6.3), extracted byte for byte with `git show`. They are the baselines that
  `validation/python/compare_confirmatory_baselines.R` compares with the
  regenerated `conf_package.json`. They are not locked or pre-registered files.

See `validation/python/README.md` for the full narrative and
`validation/HASHES_2026-11-23_v0.7.0.txt` for the current manifest.
