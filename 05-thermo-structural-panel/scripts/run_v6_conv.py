#!/usr/bin/python3
"""V6 mesh convergence: linear buckling factor of the membrane-restrained plate versus mesh.

Usage: /usr/bin/python3 scripts/run_v6_conv.py [--with80]
For each case (in-plane elements n per side, through-thickness elements nz) the script generates the mesh
(make_panel_mesh.py), the membrane restraint (frame_constraints.py --faces xy), the mechanical equivalent of a 1 K rise
(frame_prestrain.py), writes a two-step deck (preload, then *BUCKLE 6), runs CalculiX and reads the buckling factors.
Preload check: the total x reaction on face x = 0 must equal E alpha / (1 - nu) * h * W (tolerance 0.1 %).
Reference (docs/specification.md, Amendment 1): dT_cr = 2 pi^2 h^2 / (12 (1 + nu) alpha b^2).
Grid convergence index (Roache) from n = 10, 20, 40 at constant ratio r = 2 for each nz; monotone convergence required.
Pass criteria (fixed before the run): all runs finish without warnings; preload reaction within 0.1 % on every mesh;
n = 40, nz = 2 within 5 % of the closed form (the V6 tolerance); monotone convergence for nz = 2 so that the GCI exists.
Writes results/v6_convergence.txt and results/v6_convergence.csv. Exit status 1 on FAIL.
"""
import argparse
import math
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ccxio
from check_v6 import read_factors
from material import MATERIAL

CCX = "/home/harsh/miniconda3/envs/ccx223/bin/ccx"
PY = "/usr/bin/python3"
L = 0.5
H = 0.002
SERIES = [(10, 2), (20, 2), (40, 2), (10, 4), (20, 4), (40, 4)]
TOL = 0.05
RF_TOL = 0.001


def closed_form():
    m = MATERIAL
    return 2 * math.pi ** 2 * H ** 2 / (12 * (1 + m["nu"]) * m["alpha"] * L ** 2)


def call(args):
    r = subprocess.run(args, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit("command failed: %s\n%s%s" % (" ".join(args), r.stdout, r.stderr))
    return r.stdout


def deck(tag):
    m = MATERIAL
    return ("*HEADING\n"
            "Project 05 V6c convergence case %s: linear buckling, preload = mechanical equivalent of 1 K\n"
            "*INCLUDE, INPUT=../../mesh/%s.inp\n"
            "*INCLUDE, INPUT=../../mesh/%s_framepre.inp\n"
            "*NSET, NSET=FRAME_ALL\n"
            "FRAME_X0, FRAME_X1, FRAME_Y0, FRAME_Y1\n"
            "*MATERIAL, NAME=TI64\n"
            "*ELASTIC\n"
            "%.6g, %.6g\n"
            "*SOLID SECTION, ELSET=PANEL, MATERIAL=TI64\n"
            "*BOUNDARY\n"
            "FRAME_ALL, 3, 3\n"
            "*STEP\n"
            "*STATIC\n"
            "1., 1., 1.E-5, 1.\n"
            "*INCLUDE, INPUT=../../mesh/%s_framepre_bc.inp\n"
            "*NODE PRINT, NSET=FACE_X0, TOTALS=ONLY\n"
            "RF\n"
            "*END STEP\n"
            "*STEP\n"
            "*BUCKLE\n"
            "6\n"
            "*NODE PRINT, NSET=P_CEN\n"
            "U\n"
            "*END STEP\n") % (tag, tag, tag, m["E"], m["nu"], tag)


def run_case(root, n, nz):
    tag = "v6c_n%d_z%d" % (n, nz)
    call([PY, os.path.join(HERE, "make_panel_mesh.py"), tag, "--nx", str(n), "--ny", str(n), "--nz", str(nz)])
    call([PY, os.path.join(HERE, "frame_constraints.py"), tag, "--faces", "xy"])
    call([PY, os.path.join(HERE, "frame_prestrain.py"), tag])
    nnodes = len(ccxio.read_nodes(os.path.join(root, "mesh", tag + ".inp")))
    info = os.path.join(root, "mesh", tag + "_info.txt")
    quality = []
    if os.path.exists(info):
        with open(info) as f:
            quality = [s.strip() for s in f if any(k in s.lower() for k in ("aspect", "sj", "jacobian"))]
    case = os.path.join(root, "simulation", tag)
    os.makedirs(case, exist_ok=True)
    with open(os.path.join(case, tag + ".inp"), "w") as f:
        f.write(deck(tag))
    t0 = time.time()
    with open(os.path.join(case, "ccx.log"), "w") as log:
        rc = subprocess.run([CCX, "-i", tag], cwd=case, env=dict(os.environ, OMP_NUM_THREADS="4"),
                            stdout=log, stderr=subprocess.STDOUT).returncode
    secs = time.time() - t0
    with open(os.path.join(case, "ccx.log")) as f:
        logtxt = f.read()
    if rc != 0 or "Job finished" not in logtxt:
        raise SystemExit("ccx failed for %s (exit %d); see %s" % (tag, rc, os.path.join(case, "ccx.log")))
    dat = os.path.join(case, tag + ".dat")
    factors = read_factors(dat)
    if len(factors) != 6:
        raise SystemExit("%s: expected 6 factors, found %d" % (tag, len(factors)))
    blocks = [b for b in ccxio.parse_dat(dat) if b["title"].lower().startswith("total force") and b["set"] == "FACE_X0"]
    if not blocks:
        raise SystemExit("%s: no reaction block for FACE_X0" % tag)
    rf = float(blocks[0]["rows"][0][0])
    m = MATERIAL
    rf_exp = m["E"] * m["alpha"] / (1 - m["nu"]) * H * L
    return {"n": n, "nz": nz, "tag": tag, "nodes": nnodes, "secs": secs, "warn": logtxt.count("*WARNING"),
            "rf": rf, "rf_dev": 100.0 * (abs(rf) - rf_exp) / rf_exp, "f": factors, "quality": quality}


def gci(fc, fm, ff, r=2.0):
    d1, d2 = fm - fc, ff - fm
    if d1 * d2 <= 0 or d2 == 0:
        return None
    p = math.log(d1 / d2) / math.log(r)
    fext = ff + (ff - fm) / (r ** p - 1.0)
    eps = abs((ff - fm) / ff)
    return p, fext, 100.0 * 1.25 * eps / (r ** p - 1.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--with80", action="store_true")
    args = ap.parse_args()
    root = os.path.dirname(HERE)
    cases = list(SERIES)
    if args.with80:
        cases.append((80, 2))
    d1 = closed_form()
    out = ["V6c mesh convergence of the linear buckling factor, restrained 0.5 x 0.5 x 0.002 m Ti-6Al-4V plate, C3D20R",
           "closed form dT_cr(1,1) = %.5f K; expected edge reaction E alpha/(1-nu) h W = %.3f N" %
           (d1, MATERIAL["E"] * MATERIAL["alpha"] / (1 - MATERIAL["nu"]) * H * L),
           "   n  nz    nodes  time[s] warn   RF_x[N]  RFdev[%]    mode1[K]    mode2[K]    mode3[K]   err1[%]"]
    print("\n".join(out), flush=True)
    res = {}
    for n, nz in cases:
        r = run_case(root, n, nz)
        res[(n, nz)] = r
        line = " %3d  %2d  %7d  %7.1f %4d  %9.3f  %8.4f  %10.5f  %10.5f  %10.5f  %+8.3f" % (
            n, nz, r["nodes"], r["secs"], r["warn"], r["rf"], r["rf_dev"], r["f"][0], r["f"][1], r["f"][2],
            100.0 * (r["f"][0] - d1) / d1)
        out.append(line)
        print(line, flush=True)
        for q in r["quality"]:
            out.append("      mesh: " + q)
            print("      mesh: " + q, flush=True)
    out.append("")
    ok = True
    for nz in (2, 4):
        if all((n, nz) in res for n in (10, 20, 40)):
            fc, fm, ff = (res[(n, nz)]["f"][0] for n in (10, 20, 40))
            g = gci(fc, fm, ff)
            if g is None:
                out.append("nz = %d: not monotone (%.5f, %.5f, %.5f): no GCI" % (nz, fc, fm, ff))
                if nz == 2:
                    ok = False
            else:
                p, fext, gc = g
                out.append("nz = %d: observed order p = %.2f, Richardson value %.5f K (%+.3f %% from closed form), GCI(40) = %.3f %%"
                           % (nz, p, fext, 100 * (fext - d1) / d1, gc))
    if (40, 2) in res and (40, 4) in res:
        a, b = res[(40, 2)]["f"][0], res[(40, 4)]["f"][0]
        out.append("through-thickness effect at n = 40: nz = 4 minus nz = 2 = %+.5f K (%+.3f %%)" % (b - a, 100 * (b - a) / a))
    if any(r["warn"] for r in res.values()):
        out.append("FAIL reason: warnings in at least one run")
        ok = False
    if any(abs(r["rf_dev"]) > 100 * RF_TOL for r in res.values()):
        out.append("FAIL reason: preload reaction outside %.1f %%" % (100 * RF_TOL))
        ok = False
    if abs(100.0 * (res[(40, 2)]["f"][0] - d1) / d1) > 100 * TOL:
        out.append("FAIL reason: n = 40, nz = 2 outside %.1f %%" % (100 * TOL))
        ok = False
    out.append("RESULT: %s" % ("PASS" if ok else "FAIL"))
    print("\n".join(out[len(cases) + 3 + sum(len(r["quality"]) for r in res.values()):]))
    with open(os.path.join(root, "results", "v6_convergence.txt"), "w") as f:
        f.write("\n".join(out) + "\n")
    with open(os.path.join(root, "results", "v6_convergence.csv"), "w") as f:
        f.write("n,nz,nodes,seconds,warnings,rf_x_N,rf_dev_percent,mode1,mode2,mode3,mode4,mode5,mode6\n")
        for (n, nz), r in res.items():
            f.write("%d,%d,%d,%.1f,%d,%.4f,%.5f,%s\n" % (n, nz, r["nodes"], r["secs"], r["warn"], r["rf"], r["rf_dev"],
                                                         ",".join("%.6f" % v for v in r["f"])))
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
