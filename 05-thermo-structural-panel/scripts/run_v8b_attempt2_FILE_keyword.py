#!/usr/bin/python3
"""V8b: hand-over of the temperature history from a heat transfer run (.frd) to a structural run.

Usage: /usr/bin/python3 scripts/run_v8b.py

1. Thermal run (C3D20 mesh v8a_n20_z2, V8a-attempt-2 pulse, dt = 0.5 s) writing the .frd of every increment and the
   centre-column temperatures at every increment.
2. Structural run A (C3D20R mesh v6c_n20_z2, NLGEOM static with fixed increments, dt = 0.5 s) reading the .frd with *TEMPERATURE, FILE=.
   Plate held at three corner nodes only (rigid-body suppression), initial temperature 20 degC.
3. Structural run B (same, dt = 0.25 s) to characterise how CalculiX treats times between thermal frames.
Criteria (docs/specification.md, Amendment 7):
  1. run A, all nodes at t = 5, 10, ..., 30 s: max |T_struct - T_thermal| <= 0.005 K
  2. run B, centre column at multiples of 0.5 s: <= 0.005 K
  3. run B, other times from 0.75 s: deviation from linear / previous / next rule reported (no pass/fail)
  4. thermal run reproduces V8a attempt 2 at the centre column (<= 0.005 K) and all runs have zero warnings
Writes results/v8b_transfer.txt; exits with status 1 on failure.
"""
import math
import os
import re
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ccxio
import run_v8a as v8a
from material import MATERIAL

CCX, PY = v8a.CCX, v8a.PY
MESH_T = "v8a_n20_z2"
MESH_S = "v6c_n20_z2"
DT_T = 0.5
TOL = 0.005


def thermal_deck():
    m = MATERIAL
    return ("*HEADING\n"
            "Project 05 V8b: thermal run for the hand-over test, dt = %g s\n"
            "*INCLUDE, INPUT=../../mesh/%s.inp\n"
            "*MATERIAL, NAME=TI64\n"
            "*CONDUCTIVITY\n%.6g\n"
            "*DENSITY\n%.6g\n"
            "*SPECIFIC HEAT\n%.6g\n"
            "*SOLID SECTION, ELSET=PANEL, MATERIAL=TI64\n"
            "*INITIAL CONDITIONS, TYPE=TEMPERATURE\n"
            "NALL, %g\n"
            "*AMPLITUDE, NAME=PULSE\n"
            "0., 0., %g, 1., %g, 1.\n"
            "*STEP, INC=2000\n"
            "*HEAT TRANSFER, DIRECT\n"
            "%g, %g\n"
            "*INCLUDE, INPUT=flux.inp\n"
            "*NODE FILE\n"
            "NT\n"
            "*NODE PRINT, NSET=COL_CEN\n"
            "NT\n"
            "*END STEP\n") % (DT_T, MESH_T, m["k"], m["rho"], m["cp"], v8a.T0, v8a.T_RAMP, v8a.T_END, DT_T, v8a.T_END)


def struct_deck(dt, freq, nset, rb):
    m = MATERIAL
    return ("*HEADING\n"
            "Project 05 V8b: structural run reading temperatures from thermal.frd, dt = %g s\n"
            "*INCLUDE, INPUT=../../mesh/%s.inp\n"
            "*NSET, NSET=RB_A\n%d,\n"
            "*NSET, NSET=RB_B\n%d,\n"
            "*NSET, NSET=RB_D\n%d,\n"
            "*MATERIAL, NAME=TI64\n"
            "*ELASTIC\n%.9g, %.4g\n"
            "*EXPANSION, ZERO=%g\n%.6g\n"
            "*SOLID SECTION, ELSET=PANEL, MATERIAL=TI64\n"
            "*INITIAL CONDITIONS, TYPE=TEMPERATURE\n"
            "NALL, %g\n"
            "*BOUNDARY\n"
            "RB_A, 1, 3\n"
            "RB_B, 2, 3\n"
            "RB_D, 3, 3\n"
            "*STEP, NLGEOM, INC=2000\n"
            "*STATIC\n"
            "%g, %g, 1.E-3, %g\n"
            "*TEMPERATURE, FILE=thermal.frd\n"
            "*NODE PRINT, NSET=%s, FREQUENCY=%d\n"
            "NT\n"
            "*END STEP\n") % (dt, MESH_S, rb[0], rb[1], rb[2], m["E"], m["nu"], v8a.T0, m["alpha"], v8a.T0, dt, v8a.T_END, dt, nset, freq)


def run_ccx(case, job):
    with open(os.path.join(case, "ccx.log"), "w") as log:
        rc = subprocess.run([CCX, "-i", job], cwd=case, env=dict(os.environ, OMP_NUM_THREADS="4"),
                            stdout=log, stderr=subprocess.STDOUT).returncode
    with open(os.path.join(case, "ccx.log")) as f:
        txt = f.read()
    if rc != 0 or "Job finished" not in txt:
        raise SystemExit("ccx failed for %s (exit %d); log %s; last lines:\n%s"
                         % (job, rc, os.path.join(case, "ccx.log"), "\n".join(txt.splitlines()[-18:])))
    return txt.count("*WARNING")


def temp_blocks(dat, setname):
    bl = [b for b in ccxio.parse_dat(dat) if b["title"].lower().startswith("temperatures") and b["set"] == setname]
    out = {}
    for b in sorted(bl, key=lambda blk: blk["time"]):
        rows = np.array(b["rows"])
        o = np.argsort(rows[:, 0])
        out[round(b["time"], 6)] = (rows[o, 0].astype(int), rows[o, 1])
    return out


def maxdiff(a, b):
    if not np.array_equal(a[0], b[0]):
        raise SystemExit("node id lists differ between compared blocks")
    return float(np.abs(a[1] - b[1]).max())


def corner(nodes, x, y):
    hit = [n for n, p in nodes.items() if abs(p[0] - x) < 1e-9 and abs(p[1] - y) < 1e-9 and abs(p[2]) < 1e-9]
    if len(hit) != 1:
        raise SystemExit("corner (%g, %g, 0) found %d times" % (x, y, len(hit)))
    return int(hit[0])


def main():
    root = os.path.dirname(HERE)
    out = ["V8b: temperature hand-over, thermal C3D20 (%s) -> structural C3D20R (%s)" % (MESH_T, MESH_S)]
    mp_t = os.path.join(root, "mesh", MESH_T + ".inp")
    mp_s = os.path.join(root, "mesh", MESH_S + ".inp")
    par = ccxio.read_params(mp_s)
    nodes = ccxio.read_nodes(mp_s)
    rb = (corner(nodes, 0.0, 0.0), corner(nodes, par["L"], 0.0), corner(nodes, 0.0, par["W"]))
    out.append("rigid-body supports: node %d (ux,uy,uz) at (0,0), node %d (uy,uz) at (L,0), node %d (uz) at (0,W)" % rb)
    warns = {}

    # 1. thermal run
    job_t = "v8b_thermal"
    case_t = os.path.join(root, "simulation", job_t)
    os.makedirs(case_t, exist_ok=True)
    r = subprocess.run([PY, os.path.join(HERE, "make_flux_top.py"), mp_t, str(v8a.Q0), os.path.join(case_t, "flux.inp"), "PULSE",
                        "--gauss", repr(v8a.SIGMA), str(v8a.XC), str(v8a.YC)], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit("make_flux_top.py failed: %s%s" % (r.stdout, r.stderr))
    out.append(r.stdout.strip())
    with open(os.path.join(case_t, job_t + ".inp"), "w") as f:
        f.write(thermal_deck())
    warns["thermal"] = run_ccx(case_t, job_t)
    frd = os.path.join(case_t, job_t + ".frd")
    if not os.path.exists(frd):
        raise SystemExit("thermal run wrote no .frd file")
    out.append("thermal run done: %d warning(s), .frd %.1f MB" % (warns["thermal"], os.path.getsize(frd) / 1e6))
    print("\n".join(out), flush=True)
    th = temp_blocks(os.path.join(case_t, job_t + ".dat"), "COL_CEN")
    ref_col = temp_blocks(os.path.join(root, "simulation", "v8a2_n20_dt0p5", "v8a2_n20_dt0p5.dat"), "COL_CEN")
    ref_all = temp_blocks(os.path.join(root, "simulation", "v8a2_n20_dt0p5", "v8a2_n20_dt0p5.dat"), "NALL")
    d4 = max(maxdiff(th[t], ref_col[t]) for t in ref_col)
    out.append("criterion 4a: thermal run versus V8a attempt 2, centre column at %s s: max |dT| = %.3e K" % (sorted(ref_col), d4))

    # 2. structural runs
    res = {}
    for label, dt, freq, nset in (("A", 0.5, 10, "NALL"), ("B", 0.25, 1, "COL_CEN")):
        job = "v8b_struct_%s" % label.lower()
        case = os.path.join(root, "simulation", job)
        os.makedirs(case, exist_ok=True)
        shutil.copy(frd, os.path.join(case, "thermal.frd"))
        with open(os.path.join(case, job + ".inp"), "w") as f:
            f.write(struct_deck(dt, freq, nset, rb))
        warns[label] = run_ccx(case, job)
        res[label] = temp_blocks(os.path.join(case, job + ".dat"), nset)
        out.append("structural run %s (dt = %g s): %d warning(s), %d temperature blocks on %s" % (label, dt, warns[label], len(res[label]), nset))
        print(out[-1], flush=True)

    # criterion 1: run A all nodes
    a_times = sorted(ref_all)
    if sorted(res["A"]) != a_times:
        raise SystemExit("run A output times %s differ from the thermal NALL times %s" % (sorted(res["A"]), a_times))
    c1v = max(maxdiff(res["A"][t], ref_all[t]) for t in a_times)
    out.append("criterion 1: run A, all %d nodes at t = %s s: max |T_struct - T_thermal| = %.3e K (limit %.3f K)" % (len(ref_all[a_times[0]][0]), a_times, c1v, TOL))

    # criteria 2 and 3: run B
    ids0 = th[sorted(th)[0]][0]
    TT = dict(th)
    TT[0.0] = (ids0, np.full(len(ids0), v8a.T0))
    aligned, dev = [], {"linear": [], "previous": [], "next": []}
    for t, blk in sorted(res["B"].items()):
        k = t / DT_T
        if abs(k - round(k)) < 1e-6:
            aligned.append(maxdiff(blk, TT[round(round(k) * DT_T, 6)]))
        elif t >= 0.75:
            k0 = int(math.floor(k))
            t0, t1 = round(k0 * DT_T, 6), round((k0 + 1) * DT_T, 6)
            f = (t - t0) / DT_T
            a, b = TT[t0], TT[t1]
            dev["linear"].append(maxdiff(blk, (a[0], (1 - f) * a[1] + f * b[1])))
            dev["previous"].append(maxdiff(blk, a))
            dev["next"].append(maxdiff(blk, b))
    c2v = max(aligned)
    out.append("criterion 2: run B, centre column at %d multiples of 0.5 s: max |dT| = %.3e K (limit %.3f K)" % (len(aligned), c2v, TOL))
    n_un = len(dev["linear"])
    out.append("criterion 3 (characterisation): run B, %d times that are not multiples of 0.5 s (from 0.75 s), max |dT| per rule:" % n_un)
    for rule in ("linear", "previous", "next"):
        out.append("    %-9s %.4e K" % (rule, max(dev[rule])))
    best = min(dev, key=lambda r: max(dev[r]))
    out.append("    best rule: %s" % best)

    c1, c2 = c1v <= TOL, c2v <= TOL
    c4 = d4 <= TOL and all(v == 0 for v in warns.values())
    out.append("criterion 1 -> %s; criterion 2 -> %s; criterion 4 (reproduction %.3e K, warnings %s) -> %s"
               % ("pass" if c1 else "FAIL", "pass" if c2 else "FAIL", d4, warns, "pass" if c4 else "FAIL"))
    ok = c1 and c2 and c4
    out.append("RESULT: %s" % ("PASS" if ok else "FAIL"))
    print("\n".join(out[-9:]))
    with open(os.path.join(root, "results", "v8b_transfer.txt"), "w") as f:
        f.write("\n".join(out) + "\n")
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
