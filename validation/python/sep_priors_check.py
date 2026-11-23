"""
sep_priors_check.py -- compares dtasamplesize 0.7.0 `bam_sample_size()`
(method = "exact") with the Python reimplementation in ref_bam_sep_priors.py
for separate design/analysis priors, `width_type` and `assurance_for`.

For each case of a battery (>= 12) the package is called as a BLACK BOX
through `Rscript` in a subprocess (R_LIBS pointing at the library with
dtasamplesize 0.7.0) and returns: N_total, joint_assurance (= A(N_total)),
max_assurance_evaluated / N_at_max_assurance, and A(N) on a grid of fixed N
obtained with `N_range` = single value (field `max_assurance_evaluated`),
plus, for the pilot-prior cases, the priors returned by `bam_prior_from_pilot()`.
Python recomputes every one of them.  Criterion:
|A_pkg - A_py| <= 1e-10 and N_total identical.

Usage:  python -I sep_priors_check.py [output_dir]
        (default output_dir = validation/output/ next to this folder)
        Optional environment variables: RSCRIPT, DTA_PKG_LIB, DTA_EXTRA_LIB
        (see below); without them the installed dtasamplesize and the
        Rscript on the PATH are used.
Outputs: sep_priors_check.log and sep_priors_check.csv in output_dir
         (plus the generated R driver and raw package output for inspection).

PROVENANCE: written by an LLM-based coding assistant under the
supervision of the authors, from the Rd pages and NEWS of dtasamplesize 0.7.0
without reading the R code of the package.  It is a check by the same team,
not an implementation by an external group.
"""

import csv
import math
import os
import subprocess
import sys
import time

sys.stdout.reconfigure(encoding="utf-8")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ref_bam                                   # noqa: E402
import ref_bam_sep_priors as rsp                 # noqa: E402

# Machine-specific settings come from the environment, with portable defaults:
#   RSCRIPT       path to Rscript (default: "Rscript" on the PATH)
#   DTA_PKG_LIB   R library holding the dtasamplesize build to check
#                 (default: none, i.e. R's own library paths)
#   DTA_EXTRA_LIB extra R library for dependencies (default: none)
RSCRIPT = os.environ.get("RSCRIPT", "Rscript")
LIB_PKG = os.environ.get("DTA_PKG_LIB", "")
LIB_CRAN = os.environ.get("DTA_EXTRA_LIB", "")
TOL = 1e-10
NMAX = 1000          # Python cache size and N_max passed to the package

P1 = dict(prior_se=(17.0, 3.0), prior_sp=(2.0, 2.0), prior_prev=(4.0, 16.0),
          delta_se=0.14, delta_sp=0.10)
GRID = [2, 5, 20, 100, 250]


def case(cid, desc, base, **kw):
    d = dict(id=cid, desc=desc, analysis_prior_se=None, analysis_prior_sp=None,
             width_type="full_width", assurance_for="joint", alpha_ci=0.95,
             target=0.8, pilot=None, grid=list(GRID))
    d.update(base)
    d.update(kw)
    return d


VAP_PILOT = dict(tp=16, fn=1, tn=20, fp=35, base_prior_se=(9.9, 1.1),
                 base_prior_sp=(1.0, 1.0), base_prior_prev=(12.0, 43.0))

CASES = [
    case("C01", "P1, analysis = design (matches ref_bam), N=678", P1),
    case("C02", "P1 design, Jeffreys Beta(0.5,0.5) analysis both arms", P1,
         analysis_prior_se=(0.5, 0.5), analysis_prior_sp=(0.5, 0.5)),
    case("C03", "P1 design, uniform Beta(1,1) analysis both arms", P1,
         analysis_prior_se=(1.0, 1.0), analysis_prior_sp=(1.0, 1.0)),
    case("C04", "P1 design, informative analysis Se(30,5) Sp(10,2)", P1,
         analysis_prior_se=(30.0, 5.0), analysis_prior_sp=(10.0, 2.0)),
    case("C05", "P1 design, only Se analysis Jeffreys (Sp analysis NULL)", P1,
         analysis_prior_se=(0.5, 0.5)),
    case("C06", "vague design Beta(1,1)x2, informative analysis = P1 priors",
         dict(P1, prior_se=(1.0, 1.0), prior_sp=(1.0, 1.0)),
         analysis_prior_se=(17.0, 3.0), analysis_prior_sp=(2.0, 2.0)),
    case("C07", "lower_half_width, Williamson example (25,7)(30,2) prev(13.56,122.06) 0.10/0.05",
         dict(prior_se=(25.0, 7.0), prior_sp=(30.0, 2.0),
              prior_prev=(13.56, 122.06), delta_se=0.10, delta_sp=0.05),
         width_type="lower_half_width"),
    case("C08", "lower_half_width + uniform analysis priors (C07 design)",
         dict(prior_se=(25.0, 7.0), prior_sp=(30.0, 2.0),
              prior_prev=(13.56, 122.06), delta_se=0.10, delta_sp=0.05),
         width_type="lower_half_width",
         analysis_prior_se=(1.0, 1.0), analysis_prior_sp=(1.0, 1.0)),
    case("C09", "assurance_for='se', VAP Se(25.9,2.1) prev(29,98) delta_se 0.16",
         dict(prior_se=(25.9, 2.1), prior_sp=(2.0, 2.0), prior_prev=(29.0, 98.0),
              delta_se=0.16, delta_sp=0.10),
         assurance_for="se"),
    case("C10", "assurance_for='se', P1 design, Jeffreys Se analysis", P1,
         assurance_for="se", analysis_prior_se=(0.5, 0.5)),
    case("C11", "assurance_for='sp', P1 (analysis = design)", P1,
         assurance_for="sp"),
    case("C12", "assurance_for='sp' + lower_half_width + uniform Sp analysis",
         dict(prior_se=(25.0, 7.0), prior_sp=(30.0, 2.0),
              prior_prev=(13.56, 122.06), delta_se=0.10, delta_sp=0.05),
         assurance_for="sp", width_type="lower_half_width",
         analysis_prior_sp=(1.0, 1.0)),
    case("C13", "pilot VAP a0=1 cohort use_in='both' (priors=)",
         dict(delta_se=0.16, delta_sp=0.20),
         pilot=dict(VAP_PILOT, a0=1.0, design="cohort", use_in="both")),
    case("C14", "pilot VAP a0=0.5 cohort use_in='design' (base as analysis)",
         dict(delta_se=0.16, delta_sp=0.25),
         pilot=dict(VAP_PILOT, a0=0.5, design="cohort", use_in="design")),
    case("C15", "pilot case_control a0=1 use_in='both' (prevalence not updated)",
         dict(delta_se=0.20, delta_sp=0.20),
         pilot=dict(tp=30, fn=5, tn=40, fp=25, base_prior_se=(1.0, 1.0),
                    base_prior_sp=(1.0, 1.0), base_prior_prev=(12.0, 43.0),
                    a0=1.0, design="case_control", use_in="both")),
    case("C16", "alpha_ci=0.90, target 0.90, P1 design, Jeffreys analysis", P1,
         alpha_ci=0.90, target=0.9,
         analysis_prior_se=(0.5, 0.5), analysis_prior_sp=(0.5, 0.5)),
    case("C17", "infeasible within N_max=1000 (deltas 0.02): no-crossing path",
         dict(P1, delta_se=0.02, delta_sp=0.02)),
]


# --------------------------------------------------------------------------
# R driver (generated)
# --------------------------------------------------------------------------
def rnum(v):
    return repr(float(v)) if not isinstance(v, int) else repr(v) + "L"


def rvec(v):
    return "c(" + ", ".join(repr(float(t)) for t in v) + ")"


def r_args(c):
    parts = []
    if c["pilot"] is None:
        parts += [f"prior_se={rvec(c['prior_se'])}", f"prior_sp={rvec(c['prior_sp'])}",
                  f"prior_prev={rvec(c['prior_prev'])}"]
        if c["analysis_prior_se"] is not None:
            parts.append(f"analysis_prior_se={rvec(c['analysis_prior_se'])}")
        if c["analysis_prior_sp"] is not None:
            parts.append(f"analysis_prior_sp={rvec(c['analysis_prior_sp'])}")
    parts += [f"delta_se={repr(float(c['delta_se']))}",
              f"delta_sp={repr(float(c['delta_sp']))}",
              f"target_assurance={repr(float(c['target']))}",
              f"alpha_ci={repr(float(c['alpha_ci']))}",
              f"width_type={c['width_type']!r}".replace("'", '"'),
              f"assurance_for={c['assurance_for']!r}".replace("'", '"')]
    return "list(" + ", ".join(parts) + ")"


def r_pilot(c):
    p = c["pilot"]
    if p is None:
        return "NULL"
    return ("list(tp=%dL, fn=%dL, tn=%dL, fp=%dL, base_prior_se=%s, base_prior_sp=%s, "
            "base_prior_prev=%s, a0=%s, design=\"%s\", use_in=\"%s\")" % (
                p["tp"], p["fn"], p["tn"], p["fp"], rvec(p["base_prior_se"]),
                rvec(p["base_prior_sp"]), rvec(p["base_prior_prev"]),
                repr(float(p["a0"])), p["design"], p["use_in"]))


R_HEAD = r'''
suppressMessages(library(dtasamplesize))
emit <- function(id, key, val) cat(sprintf("%s|%s|%s\n", id, key, val))
num <- function(x) if (is.na(x)) "NA" else sprintf("%.17g", x)
emit("META", "version", as.character(packageVersion("dtasamplesize")))
emit("META", "libpath", find.package("dtasamplesize"))
run_case <- function(id, args, pilot, grid) {
  if (!is.null(pilot)) {
    pp <- do.call(bam_prior_from_pilot, pilot)
    for (nm in c("prior_se","prior_sp","prior_prev","analysis_prior_se","analysis_prior_sp")) {
      v <- pp[[nm]]
      emit(id, paste0("pilot_", nm), paste(sprintf("%.17g", v), collapse = ","))
    }
    args$priors <- pp
  }
  base <- c(list(n_range = 20:60, B = 1000, N_max = %NMAX%L), args)
  r <- suppressWarnings(do.call(bam_sample_size, base))
  emit(id, "target_reached", as.integer(r$target_reached))
  emit(id, "N_total", num(r$N_total))
  emit(id, "joint_assurance", num(r$joint_assurance))
  emit(id, "max_assurance_evaluated", num(r$max_assurance_evaluated))
  emit(id, "N_at_max_assurance", num(r$N_at_max_assurance))
  emit(id, "search_type", r$search_type)
  emit(id, "eff_analysis_prior_se", paste(sprintf("%.17g", r$analysis_prior_se), collapse = ","))
  emit(id, "eff_analysis_prior_sp", paste(sprintf("%.17g", r$analysis_prior_sp), collapse = ","))
  Ns <- grid
  if (isTRUE(r$target_reached)) Ns <- c(Ns, r$N_total + (-1:1))
  Ns <- sort(unique(Ns[Ns >= 2]))
  for (N in Ns) {
    rr <- suppressWarnings(do.call(bam_sample_size, c(base, list(N_range = N))))
    emit(id, paste0("A@", N), num(rr$max_assurance_evaluated))
  }
}
'''


def build_driver():
    lines = [R_HEAD.replace("%NMAX%", str(NMAX))]
    for c in CASES:
        grid = "c(" + ", ".join(f"{g}L" for g in c["grid"]) + ")"
        lines.append(f'run_case("{c["id"]}", {r_args(c)}, {r_pilot(c)}, {grid})')
    lines.append('emit("META", "done", "1")')
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------
def parse(text):
    out = {}
    for ln in text.splitlines():
        if ln.count("|") >= 2:
            cid, key, val = ln.split("|", 2)
            out.setdefault(cid, {})[key] = val
    return out


def fnum(s):
    return float("nan") if s == "NA" else float(s)


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else os.path.normpath(
        os.path.join(HERE, "..", "output"))
    os.makedirs(outdir, exist_ok=True)
    logf = open(os.path.join(outdir, "sep_priors_check.log"), "w", encoding="utf-8")

    def log(s=""):
        print(s)
        logf.write(s + "\n")
        logf.flush()

    t0 = time.time()
    log("sep_priors_check.py -- dtasamplesize 0.7.0 vs ref_bam_sep_priors.py")
    log(f"started {time.strftime('%Y-%m-%d %H:%M:%S')}  tolerance |dA| <= {TOL:g}, N identical")
    log(f"Python {sys.version.split()[0]}")

    # ---- package side (black box) ----
    drv = os.path.join(outdir, "sep_priors_pkg_driver.R")
    with open(drv, "w", encoding="utf-8") as f:
        f.write(build_driver())
    env = dict(os.environ)
    libs = [p for p in (LIB_PKG, LIB_CRAN) if p]
    if libs:
        env["R_LIBS"] = os.pathsep.join(libs)
        env.pop("R_LIBS_USER", None)
    log(f"Rscript driver: {drv}")
    log(f"Rscript = {RSCRIPT}")
    log(f"R_LIBS = {env.get('R_LIBS', '(R default library paths)')}")
    t1 = time.time()
    p = subprocess.run([RSCRIPT, "--vanilla", drv], env=env, capture_output=True,
                       text=True, encoding="utf-8")
    log(f"Rscript exit code {p.returncode}, {time.time() - t1:.1f} s")
    with open(os.path.join(outdir, "sep_priors_pkg_raw.txt"), "w", encoding="utf-8") as f:
        f.write(p.stdout)
        f.write("\n--- STDERR ---\n" + p.stderr)
    if p.returncode != 0 or "META|done|1" not in p.stdout:
        log("R FAILED:\n" + p.stderr[-3000:])
        sys.exit(2)
    pkg = parse(p.stdout)
    log(f"package version: {pkg['META']['version']}   installed at: {pkg['META']['libpath']}")

    # ---- Python side + comparison ----
    rows = []
    n_fail = 0

    def add(cid, quantity, N, a_pkg, a_py, kind="A"):
        nonlocal n_fail
        if kind == "N":
            ok = (a_pkg == a_py) or (math.isnan(a_pkg) and math.isnan(a_py))
            diff = float("nan") if (math.isnan(a_pkg) or math.isnan(a_py)) else a_pkg - a_py
        else:
            diff = abs(a_pkg - a_py)
            ok = diff <= TOL
        if not ok:
            n_fail += 1
        rows.append(dict(case=cid, quantity=quantity, N=N, package=a_pkg,
                         python=a_py, abs_diff=diff, ok=int(ok)))
        return ok

    summary = []
    for c in CASES:
        cid = c["id"]
        r = pkg[cid]
        tc = time.time()
        # priors
        if c["pilot"] is None:
            dse, dsp, dprev = c["prior_se"], c["prior_sp"], c["prior_prev"]
            anse, ansp = c["analysis_prior_se"], c["analysis_prior_sp"]
        else:
            pp = c["pilot"]
            u = rsp.prior_from_pilot(pp["tp"], pp["fn"], pp["tn"], pp["fp"],
                                     pp["base_prior_se"], pp["base_prior_sp"],
                                     pp["base_prior_prev"], pp["a0"], pp["design"],
                                     pp["use_in"])
            dse, dsp, dprev = u["prior_se"], u["prior_sp"], u["prior_prev"]
            anse, ansp = u["analysis_prior_se"], u["analysis_prior_sp"]
            for nm in ("prior_se", "prior_sp", "prior_prev",
                       "analysis_prior_se", "analysis_prior_sp"):
                pk = [float(t) for t in r["pilot_" + nm].split(",")]
                py = list(u[nm])
                for j in (0, 1):
                    add(cid, f"pilot_update:{nm}[{j}]", "", pk[j], py[j])
        eff_se = [float(t) for t in r["eff_analysis_prior_se"].split(",")]
        eff_sp = [float(t) for t in r["eff_analysis_prior_sp"].split(",")]
        exp_se = anse if anse is not None else dse
        exp_sp = ansp if ansp is not None else dsp
        for j in (0, 1):
            add(cid, f"effective_analysis_prior_se[{j}]", "", eff_se[j], float(exp_se[j]))
            add(cid, f"effective_analysis_prior_sp[{j}]", "", eff_sp[j], float(exp_sp[j]))

        m = rsp.BAMExactSep(dse, dsp, dprev, c["delta_se"], c["delta_sp"],
                            level=c["alpha_ci"], nmax=NMAX, analysis_prior_se=anse,
                            analysis_prior_sp=ansp, width_type=c["width_type"],
                            assurance_for=c["assurance_for"])
        s = m.search(c["target"], 2, NMAX)
        reached_pkg = int(r["target_reached"])
        add(cid, "target_reached", "", float(reached_pkg), float(s["target_reached"]), "N")
        add(cid, "N_total", "", fnum(r["N_total"]),
            float(s["N_total"]) if s["N_total"] is not None else float("nan"), "N")
        if s["target_reached"]:
            add(cid, "joint_assurance=A(N_total)", s["N_total"],
                fnum(r["joint_assurance"]), s["assurance"])
        else:
            add(cid, "max_assurance_evaluated", "", fnum(r["max_assurance_evaluated"]),
                s["max_assurance_evaluated"])
            add(cid, "N_at_max_assurance", "", fnum(r["N_at_max_assurance"]),
                float(s["N_at_max"]), "N")
        # grid of fixed N (N_range of a single value)
        for key in sorted((k for k in r if k.startswith("A@")),
                          key=lambda k: int(k[2:])):
            N = int(key[2:])
            add(cid, "A(N) fixed-N (max_assurance_evaluated)", N, fnum(r[key]),
                m.assurance(N))
        # C01: also against ref_bam (unextended)
        if cid == "C01":
            rb = ref_bam.BAMExact(dse, dsp, dprev, c["delta_se"], c["delta_sp"],
                                  level=c["alpha_ci"], nmax=NMAX)
            Nrb, Arb = rb.search(c["target"], 2, NMAX)
            add(cid, "ref_bam.BAMExact N_total vs package", "", fnum(r["N_total"]),
                float(Nrb), "N")
            add(cid, "ref_bam.BAMExact A(N_total) vs package", Nrb,
                fnum(r["joint_assurance"]), Arb)
            for key in sorted((k for k in r if k.startswith("A@")),
                              key=lambda k: int(k[2:])):
                N = int(key[2:])
                add(cid, "ref_bam.BAMExact A(N) vs extension (python-python)", N,
                    rb.assurance(N), m.assurance(N))
        sub = [x for x in rows if x["case"] == cid]
        mx = max((x["abs_diff"] for x in sub if x["quantity"].find("N_total") < 0
                  and not math.isnan(x["abs_diff"])), default=float("nan"))
        summary.append((cid, c["desc"], fnum(r["N_total"]),
                        s["N_total"], fnum(r["joint_assurance"]), s["assurance"],
                        mx, all(x["ok"] for x in sub), time.time() - tc))
        log(f"{cid} done ({time.time() - tc:.1f} s)  N_pkg={r['N_total']}  N_py={s['N_total']}")

    # ---- CSV ----
    csvp = os.path.join(outdir, "sep_priors_check.csv")
    with open(csvp, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for x in rows:
            xx = dict(x)
            for k in ("package", "python", "abs_diff"):
                xx[k] = repr(x[k]) if isinstance(x[k], float) else x[k]
            w.writerow(xx)

    # ---- report ----
    log("")
    log("=" * 118)
    log(f"{'case':<5}{'N_pkg':>7}{'N_py':>7}{'A_pkg(N)':>20}{'A_py(N)':>20}"
        f"{'max|dA| all A':>16}  ok   description")
    for cid, desc, npk, npy, apk, apy, mx, ok, _ in summary:
        npy_s = "NA" if npy is None else str(npy)
        apy_s = "NA" if apy is None else f"{apy:.12f}"
        apk_s = "NA" if math.isnan(apk) else f"{apk:.12f}"
        npk_s = "NA" if math.isnan(npk) else str(int(npk))
        log(f"{cid:<5}{npk_s:>7}{npy_s:>7}{apk_s:>20}{apy_s:>20}{mx:>16.3e}  "
            f"{'OK ' if ok else 'FAIL'} {desc}")
    log("=" * 118)
    log(f"{len(rows)} comparisons, {n_fail} failures (tolerance {TOL:g}, N identical)")
    big = max((x["abs_diff"] for x in rows if not math.isnan(x["abs_diff"])
               and not x["quantity"].startswith("ref_bam.BAMExact A(N) vs extension")),
              default=0.0)
    log(f"largest |difference| over all package-vs-Python comparisons: {big:.3e}")
    pr = [x for x in rows if x["quantity"].startswith("ref_bam.BAMExact A(N) vs extension")]
    if pr:
        log("max |ref_bam - extension| (python-python, C01): "
            f"{max(x['abs_diff'] for x in pr):.3e}")
    log(f"CSV: {csvp}")
    log(f"total time {time.time() - t0:.1f} s")
    log("RESULT: " + ("ALL CASES AGREE" if n_fail == 0 else f"{n_fail} DISAGREEMENTS"))
    logf.close()
    sys.exit(0 if n_fail == 0 else 1)


if __name__ == "__main__":
    main()
