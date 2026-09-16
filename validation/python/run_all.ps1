# Run the full development + confirmatory validation pipeline for the
# independent Python reference implementation of dtasamplesize::bam_sample_size
# (method = "exact") against the installed R package.
#
# Always cd's into its own directory first, so it can be invoked from anywhere.
# Runs dev_* then conf_* scripts in dependency order (see README.md: conf_R_01
# must run before conf_PY_02 / conf_R_07 / conf_PY_08, all of which read its
# output conf_package.json).
#
# Each script's combined stdout+stderr is teed to logs/<script>.log. On the
# first non-zero exit code this stops and reports which script failed.

Set-StrictMode -Version Latest
Set-Location -Path $PSScriptRoot

Remove-Item Env:\LC_ALL -ErrorAction SilentlyContinue
Remove-Item Env:\LC_CTYPE -ErrorAction SilentlyContinue
Remove-Item Env:\LANG -ErrorAction SilentlyContinue

$RSCRIPT = "C:/Program Files/R/R-4.5.2/bin/Rscript.exe"
$PYTHON  = "C:/Users/DELL/AppData/Local/Python/bin/python.exe"

New-Item -ItemType Directory -Force -Path "logs" | Out-Null

function Run-Step {
    param(
        [string]$Kind,
        [string]$Script
    )
    $log = "logs/$Script.log"
    Write-Host "=== [$(Get-Date -AsUTC -Format o)] START $Script ($Kind) ==="
    if ($Kind -eq "R") {
        & $RSCRIPT $Script *> $log
    } else {
        & $PYTHON $Script *> $log
    }
    $status = $LASTEXITCODE
    if ($status -ne 0) {
        Write-Host "=== FAILED: $Script (exit $status). Last 40 lines of $log : ==="
        Get-Content $log -Tail 40
        Write-Host "=== stopping run_all.ps1: $Script failed, downstream scripts not run ==="
        exit $status
    }
    Write-Host "=== [$(Get-Date -AsUTC -Format o)] OK    $Script ==="
}

Write-Host "Python: $(& $PYTHON --version 2>&1)"
Write-Host "Rscript: $((& $RSCRIPT --version 2>&1) | Select-Object -First 1)"
Write-Host ""

# -- development scripts (debugging grid; not pre-registered) --------------
Run-Step -Kind PY -Script "dev_01_agreement_smallN.py"
Run-Step -Kind PY -Script "dev_02_published.py"
Run-Step -Kind R  -Script "dev_03_probe_pkg.R"

# -- confirmatory scripts (locked grid) -------------------------------------
Run-Step -Kind R  -Script "conf_R_01_package.R"           # writes conf_package.json
Run-Step -Kind PY -Script "conf_PY_02_reference.py"        # reads conf_package.json, writes conf_reference.json
Run-Step -Kind R  -Script "conf_R_03_properties.R"         # writes conf_properties.json
Run-Step -Kind R  -Script "conf_R_04_mc_vs_exact.R"
Run-Step -Kind PY -Script "conf_PY_05_shared_assumptions.py"
Run-Step -Kind PY -Script "conf_PY_06_monotone_in_N.py"
Run-Step -Kind R  -Script "conf_R_07_reverify_version.R"   # reads conf_package.json
Run-Step -Kind PY -Script "conf_PY_08_degeneracy_bite.py"  # reads conf_package.json
Run-Step -Kind PY -Script "conf_PY_09_published_mc.py"

Write-Host ""
Write-Host "=== ALL SCRIPTS COMPLETED SUCCESSFULLY ==="
