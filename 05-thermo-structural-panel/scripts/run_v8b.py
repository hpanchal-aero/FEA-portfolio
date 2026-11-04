#!/usr/bin/python3
"""V8b: hand-over of a transient temperature history from a heat transfer run to a structural run.

Usage: /usr/bin/python3 scripts/run_v8b.py

Attempt 2 (scripts/run_v8b_attempt2_FILE_keyword.py, on record) showed that *TEMPERATURE, FILE= takes the last frame of the
thermal results and the step ramps linearly to it. Attempt 3 therefore uses one structural step per thermal frame, each with an
explicit *TEMPERATURE field (all nodes) taken from the thermal .dat; each step ramps from the previous frame to the next.

1. Thermal run (C3D20 mesh v8a_n20_z2, V8a-attempt-2 pulse, dt = 0.5 s): temperatures of all nodes at every increment.
2. Structural run A (C3D20R mesh v6c_n20_z2, NLGEOM, 60 steps of 0.5 s, one increment each), all nodes printed every step.
   Plate held at three corner nodes (rigid-body suppression), initial temperature 20 degC.
3. Structural run B (same, two increments of 0.25 s per step), centre column printed every increment.
Criteria (docs/specification.md, Amendments 7 and 9):
  1. run A, all nodes at t = 5, 10, ..., 30 s: max |T_struct - T_thermal| <= 0.005 K (all 60 frames also reported)
  2. run B, centre column at multiples of 0.5 s: <= 0.005 K
  3. run B, other times: deviation from linear / previous / next rule reported (no pass/fail)
  4. thermal run reproduces V8a attempt 2 (all nodes at 5, 10, ..., 30 s, <= 0.005 K) and all runs have zero warnings
Writes results/v8b_transfer.txt; exits with status 1 on failure.
"""
import math
import os
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
NFRAMES = int(round(v8a.T_END / DT_T))
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
            "*NODE PRINT, NSET=NALL\n"
            "NT\n"
            "*NODE PRINT, NSET=COL_CEN\n"
            "NT\n"
            "*END STEP\n") % (DT_T, MESH_T, m["k"], m["rho"], m["cp"], v8a.T0, v8a.T_RAMP, v8a.T_END, DT_T, v8a.T_END)


def struct_deck(sub_dt, nset, rb):
    m = MATERIAL
    head = ("*HEADING\n"
            "Project 05 V8b: structural run, one step per thermal frame, increment %g s\n"
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
            "RB_D, 3, 3\n") % (sub_dt, MESH_S, rb[0], rb[1], rb[2], m["E"], m["nu"], v8a.T0, m["alpha"], v8a.T0)
    steps = []
    for k in range(1, NFRAMES + 1):
        steps.append("*STEP, NLGEOM, INC=2000\n"
                     "*STATIC\n"
                     "%g, %g, 1.E-3, %g\n"
                     "*INCLUDE, INPUT=temps/frame_%03d.inp\n"
                     "*NODE PRINT, NSET=%s\n"
                     "NT\n"
                     "*END STEP\n" % (sub_dt, DT_T, sub_dt, k, nset))
    return head + "".join(steps)


def write_frames(case, thn):
    d = os.path.join(case, "temps")
    os.makedirs(d, exist_ok=True)
    for k in range(1, NFRAMES + 1):
        ids, T = thn[round(k * DT_T, 6)]
        with open(os.path.join(d, "frame_%03d.inp" % k), "w") as f:
            f.write("*TEMPERATURE\n")
            f.writelines("%d, %.7e\n" % (i, t) for i, t in zip(ids, T))


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
    sim = os.path.join(root, "simulation")
    out = ["V8b attempt 3: temperature hand-over, thermal C3D20 (%s) -> structural C3D20R (%s), one structural step per thermal frame"
           % (MESH_T, MESH_S)]
    mp_t = os.path.join(root, "mesh", MESH_T + ".inp")
    mp_s = os.path.join(root, "mesh", MESH_S + ".inp")
    par = ccxio.read_params(mp_s)
    nodes = ccxio.read_nodes(mp_s)
    rb = (corner(nodes, 0.0, 0.0), corner(nodes, par["L"], 0.0), corner(nodes, 0.0, par["W"]))
    out.append("rigid-body supports: node %d (ux,uy,uz) at (0,0), node %d (uy,uz) at (L,0), node %d (uz) at (0,W)" % rb)
    warns = {}

    job_t = "v8b_thermal"
    case_t = os.path.join(sim, job_t)
    os.makedirs(case_t, exist_ok=True)
    r = subprocess.run([PY, os.path.join(HERE, "make_flux_top.py"), mp_t, str(v8a.Q0), os.path.join(case_t, "flux.inp"), "PULSE",
                        "--gauss", repr(v8a.SIGMA), str(v8a.XC), str(v8a.YC)], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit("make_flux_top.py failed: %s%s" % (r.stdout, r.stderr))
    out.append(r.stdout.strip())
    with open(os.path.join(case_t, job_t + ".inp"), "w") as f:
        f.write(thermal_deck())
    warns["thermal"] = run_ccx(case_t, job_t)
    dat_t = os.path.join(case_t, job_t + ".dat")
    thn = temp_blocks(dat_t, "NALL")
    th = temp_blocks(dat_t, "COL_CEN")
    expect = [round(k * DT_T, 6) for k in range(1, NFRAMES + 1)]
    if sorted(thn) != expect or sorted(th) != expect:
        raise SystemExit("thermal output times unexpected: %s" % sorted(thn))
    out.append("thermal run done: %d warning(s), %d frames of %d nodes" % (warns["thermal"], len(thn), len(thn[expect[0]][0])))
    print("\n".join(out), flush=True)
    ref_dat = os.path.join(sim, "v8a2_n20_dt0p5", "v8a2_n20_dt0p5.dat")
    ref_all = temp_blocks(ref_dat, "NALL")
    d4 = max(maxdiff(thn[t], ref_all[t]) for t in ref_all)
    out.append("criterion 4a: thermal run versus V8a attempt 2, all nodes at %s s: max |dT| = %.3e K" % (sorted(ref_all), d4))

    res = {}
    for label, sub_dt, nset in (("A", 0.5, "NALL"), ("B", 0.25, "COL_CEN")):
        job = "v8b_struct_%s" % label.lower()
        case = os.path.join(sim, job)
        os.makedirs(case, exist_ok=True)
        write_frames(case, thn)
        with open(os.path.join(case, job + ".inp"), "w") as f:
            f.write(struct_deck(sub_dt, nset, rb))
        warns[label] = run_ccx(case, job)
        res[label] = temp_blocks(os.path.join(case, job + ".dat"), nset)
        out.append("structural run %s (increment %g s): %d warning(s), %d temperature blocks on %s" % (label, sub_dt, warns[label], len(res[label]), nset))
        print(out[-1], flush=True)
    missing_a = [t for t in expect if t not in res["A"]]
    extra_a = sorted(t for t in res["A"] if t not in expect)
    if missing_a or len(res["B"]) != 2 * NFRAMES:
        raise SystemExit("unexpected output: run A missing frames %s, run B has %d blocks" % (missing_a, len(res["B"])))
    out.append("run A: %d extra output blocks from subdivided steps (increment cut back by CalculiX) at t = %s s" % (len(extra_a), extra_a))

    reg = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0]
    c1v = max(maxdiff(res["A"][t], thn[t]) for t in reg)
    c1all = max(maxdiff(res["A"][t], thn[t]) for t in expect)
    out.append("criterion 1: run A, all %d nodes at t = %s s: max |T_struct - T_thermal| = %.3e K (limit %.3f K); all %d frames: %.3e K"
               % (len(thn[expect[0]][0]), reg, c1v, TOL, NFRAMES, c1all))

    ids0 = th[expect[0]][0]
    TT = dict(th)
    TT[0.0] = (ids0, np.full(len(ids0), v8a.T0))
    aligned, dev = [], {"linear": [], "previous": [], "next": []}
    for t, blk in sorted(res["B"].items()):
        k = t / DT_T
        if abs(k - round(k)) < 1e-6:
            aligned.append(maxdiff(blk, TT[round(round(k) * DT_T, 6)]))
        else:
            k0 = int(math.floor(k))
            t0, t1 = round(k0 * DT_T, 6), round((k0 + 1) * DT_T, 6)
            f = (t - t0) / DT_T
            a, b = TT[t0], TT[t1]
            dev["linear"].append(maxdiff(blk, (a[0], (1 - f) * a[1] + f * b[1])))
            dev["previous"].append(maxdiff(blk, a))
            dev["next"].append(maxdiff(blk, b))
    c2v = max(aligned)
    out.append("criterion 2: run B, centre column at %d multiples of 0.5 s: max |dT| = %.3e K (limit %.3f K)" % (len(aligned), c2v, TOL))
    out.append("criterion 3 (characterisation): run B, %d times that are not multiples of 0.5 s, max |dT| per rule:" % len(dev["linear"]))
    for rule in ("linear", "previous", "next"):
        out.append("    %-9s %.4e K" % (rule, max(dev[rule])))
    out.append("    best rule: %s" % min(dev, key=lambda q: max(dev[q])))

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
