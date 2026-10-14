#!/usr/bin/env python3
"""
Project 03, Stage 6: regenerate the solver files of selected layouts at rib t = 1.5 mm (4 mm solid / 8 mm skin,
skin 1.0 mm, Stage 4 placeholder pressure) and KEEP them, for rendering. Same scripts and arguments as
run_layout_study.py; no cleanup. Each result is compared with the driver's stored record; a mismatch stops the run.
Cases: L000 (baseline), L006 (best buckling), L050 (lowest tip deflection at 1.5 mm); buckling solve for L006 only.
Usage: /usr/bin/python3 scripts/run_render_cases.py [L000,L006,L050]
"""
import csv
import json
import os
import re
import subprocess
import sys
import time

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(PROJECT)
PY = "/usr/bin/python3"
CONDA_SH = "/home/harsh/miniconda3/etc/profile.d/conda.sh"
BUCKLE_FOR = {"L006", "L100", "L126"}
LOGDIR = "results/stage6/render_cases"
os.makedirs(LOGDIR, exist_ok=True)


def mem_avail_gb():
    for ln in open("/proc/meminfo"):
        if ln.startswith("MemAvailable:"):
            return int(ln.split()[1]) / 1048576.0


def sh(cmd, env, label, case):
    t0 = time.time()
    print("[%s] %s ..." % (case, label), flush=True)
    r = subprocess.run(cmd, env=env, capture_output=True, text=True, errors="replace")
    out = r.stdout + r.stderr
    open("%s/%s_%s.log" % (LOGDIR, case, label), "w").write(out)
    print("[%s] %s: rc=%d  %.0f s" % (case, label, r.returncode, time.time() - t0), flush=True)
    if r.returncode != 0:
        sys.exit("STOP: %s failed for %s; see %s/%s_%s.log" % (label, case, LOGDIR, case, label))
    return out


def solve(job, case, label):
    while mem_avail_gb() < 7.0:
        print("waiting for memory (%.1f GB free)" % mem_avail_gb(), flush=True)
        time.sleep(30)
    out = sh(["bash", "-c", "source %s && conda activate ccx223 && exec %s -u scripts/run_ccx_watchdog.py %s 8"
              % (CONDA_SH, PY, job)], os.environ.copy(), label, case)
    if "SOLVE OK" not in out:
        sys.exit("STOP: no 'SOLVE OK' for %s %s" % (case, label))
    print("   " + [l for l in out.splitlines() if "peak RSS" in l][-1:][0].strip() if "peak RSS" in out else "")


def main():
    ids = sys.argv[1].split(",") if len(sys.argv) > 1 else ["L000", "L006", "L050"]
    layouts = {r["id"]: r for r in csv.DictReader(open("results/stage6_layouts.csv"))}
    for lid in ids:
        L = layouts[lid]
        env = os.environ.copy()
        for k in ("RIB_SPACING", "NETGEN_KW", "RIB_STATIONS", "LAYOUT_ID", "RIB_T"):
            env.pop(k, None)
        env.update(RIB_STATIONS=L["stations"], LAYOUT_ID=lid, RIB_T="1.5")
        tag = "_" + lid
        raw = "mesh/stage3_skin/raw_s4_k8%s_closeedgefacnone.npz" % tag
        sjob = "simulation/stage3_skin_s4_k8%s_closeedgefacnone_t1_n3p8" % tag
        job = "simulation/stage4_s4_k8%s_closeedgefacnone_t1_n3p8" % tag
        print("\n===== %s  stations %s =====" % (lid, L["stations"]), flush=True)
        sh([PY, "-u", "scripts/build_stage3_solid.py", "ribs"], env, "geom", lid)
        sh([PY, "-u", "scripts/build_stage3_skin_geom.py"], env, "imprint", lid)
        env_m = dict(env, NETGEN_KW="closeedgefac=None")
        for attempt in (1, 2):
            if os.path.exists(raw):
                os.remove(raw)
            try:
                sh(["timeout", "150", PY, "-u", "scripts/mesh_stage3_skin_raw.py", "4.0", "8.0"], env_m, "mesh%d" % attempt, lid)
                if os.path.exists(raw):
                    break
            except SystemExit:
                if attempt == 2:
                    raise
                print("   mesh attempt 1 failed (known intermittent Netgen failure); retrying once", flush=True)
        sh([PY, "-u", "scripts/build_stage3_skin_inp.py", raw, "1.0", "372.78", "3.8"], env, "inp", lid)
        sh([PY, "-u", "scripts/apply_stage4_pressure_load.py", sjob + ".inp", job + ".inp", "full"], env, "load", lid)
        solve(job, lid, "solve-static")
        sh([PY, "-u", "scripts/extract_stage3_skin.py", job], env, "extract-static", lid)
        rec = json.load(open("results/stage6/records/%s.json" % lid))["metrics"]
        got = json.load(open("results/%s_results.json" % os.path.basename(job)))["uy_spar_axis_mm"]
        rel = abs(got - rec["tip_mm"]) / rec["tip_mm"]
        print("[%s] TIP CHECK: now %.6f, record %.6f, rel diff %.2e -> %s" % (lid, got, rec["tip_mm"], rel, "PASS" if rel < 1e-6 else "FAIL"))
        if rel >= 1e-6:
            sys.exit("STOP: tip does not reproduce the driver record for %s" % lid)
        if lid in BUCKLE_FOR:
            bjob = job + "_buckle"
            sh([PY, "scripts/make_buckle_inp.py", job + ".inp", bjob + ".inp", "6"], env, "buckle-make", lid)
            solve(bjob, lid, "solve-buckle")
            out = sh([PY, "scripts/extract_stage3_buckle.py", bjob], env, "extract-buckle", lid)
            crit = float(re.search(r"Critical buckling factor: ([0-9.eE+-]+)", out).group(1))
            d = abs(crit - rec["buckling_crit"])
            print("[%s] BUCKLING CHECK: now %.4f, record %.4f, abs diff %.1e -> %s" % (lid, crit, rec["buckling_crit"], d, "PASS" if d < 1e-4 else "FAIL"))
            if d >= 1e-4:
                sys.exit("STOP: buckling does not reproduce the driver record for %s" % lid)
    print("\nALL CASES REPRODUCED. Files kept in simulation/ (stage4_s4_k8_L0xx_closeedgefacnone_t1_n3p8*).")


if __name__ == "__main__":
    main()
