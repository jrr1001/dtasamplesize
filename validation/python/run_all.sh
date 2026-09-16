#!/usr/bin/env bash
# Run the full development + confirmatory validation pipeline for the
# independent Python reference implementation of dtasamplesize::bam_sample_size
# (method = "exact") against the installed R package.
#
# Always cd's into its own directory first, so it can be invoked from anywhere.
# Runs dev_* then conf_* scripts in dependency order (see README.md for why
# this order matters: conf_R_01 must run before conf_PY_02 / conf_R_07 /
# conf_PY_08, all of which read its output conf_package.json).
#
# Each script's combined stdout+stderr is teed to logs/<script>.log. The
# script's own exit status (not tee's) determines pass/fail; on the first
# failure this stops and reports which script failed and its exit code.

set -u
cd "$(dirname "${BASH_SOURCE[0]}")"

unset LC_ALL LC_CTYPE LANG
RSCRIPT="/c/Program Files/R/R-4.5.2/bin/Rscript.exe"
PYTHON="C:/Users/DELL/AppData/Local/Python/bin/python.exe"

mkdir -p logs

run_step() {
  local kind="$1" script="$2"
  local log="logs/${script}.log"
  echo "=== [$(date -u +%Y-%m-%dT%H:%M:%SZ)] START $script ($kind) ==="
  if [ "$kind" = "R" ]; then
    "$RSCRIPT" "$script" > "$log" 2>&1
  else
    "$PYTHON" "$script" > "$log" 2>&1
  fi
  local status=$?
  if [ $status -ne 0 ]; then
    echo "=== FAILED: $script (exit $status). Last 40 lines of $log: ==="
    tail -n 40 "$log"
    echo "=== stopping run_all.sh: $script failed, downstream scripts not run ==="
    exit $status
  fi
  echo "=== [$(date -u +%Y-%m-%dT%H:%M:%SZ)] OK    $script ==="
}

echo "Python: $("$PYTHON" --version 2>&1)"
echo "Rscript: $("$RSCRIPT" --version 2>&1 | head -n 1)"
echo

# -- development scripts (debugging grid; not pre-registered) --------------
run_step PY dev_01_agreement_smallN.py
run_step PY dev_02_published.py
run_step R  dev_03_probe_pkg.R

# -- confirmatory scripts (locked grid) -------------------------------------
run_step R  conf_R_01_package.R          # writes conf_package.json
run_step PY conf_PY_02_reference.py      # reads conf_package.json, writes conf_reference.json
run_step R  conf_R_03_properties.R       # writes conf_properties.json
run_step R  conf_R_04_mc_vs_exact.R
run_step PY conf_PY_05_shared_assumptions.py
run_step PY conf_PY_06_monotone_in_N.py
run_step R  conf_R_07_reverify_version.R # reads conf_package.json
run_step PY conf_PY_08_degeneracy_bite.py # reads conf_package.json
run_step PY conf_PY_09_published_mc.py

echo
echo "=== ALL SCRIPTS COMPLETED SUCCESSFULLY ==="
