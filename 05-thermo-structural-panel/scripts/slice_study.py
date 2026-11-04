#!/usr/bin/python3
"""Edge-layer study for the uniform-heating case (plane-strain slice of the plate).

Usage: /usr/bin/python3 scripts/slice_study.py [--levels 2 4 8 16] [--dT 2]

A thin slice of the plate (L = 0.5 m, h = 2 mm, W = 2 mm, one element in y, uy = 0 on both y faces)
represents the plate far from its y edges: plane strain with eps_y = 0, which for uniform heating is the
biaxially restrained state, so the interior reference is sigma_x = -E alpha dT / (1 - nu).
Level n: n elements through the thickness, 100 n along x (element size halves with each level, dx/dz = 2.5).
Each level is run with three support conditions on the x edges:
  midline:  ux = 0 only on the mid-thickness line of the x edges (single pinned line)
  full:     ux = 0 on the whole x faces (ideal restraint, the closed form is exact, edge rotation also restrained)
  membrane: thickness-mean of ux = 0 on every through-thickness line (scripts/frame_constraints.py), rotation free
uz = 0 on the mid-thickness lines of all four edges in all cases.
Reported: sigma_x over the elements whose centroid is within 50 mm of the plate centre (mean, extremes,
signed deviation of the mean from the reference, largest deviation) and the largest |sigma_z|.
Writes results/slice_study.csv and results/slice_study.txt; the .dat files are deleted after extraction.
"""
import argparse
import csv
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ccxio
from material import MATERIAL

CCX = "/home/harsh/miniconda3/envs/ccx223/bin/ccx"
PY = "/usr/bin/python3"
BC = {"midline": "FRAME_X0, 1, 1\nFRAME_X1, 1, 1\n",
      "full": "FACE_X0, 1, 1\nFACE_X1, 1, 1\n",
      "membrane": ""}


def deck(tag, dT, bc):
    m = MATERIAL
    inc = "*INCLUDE, INPUT=../../mesh/%s_frame.inp\n" % tag if bc == "membrane" else ""
    return ("*HEADING\n"
            "Project 05 edge-layer study: plane-strain slice, uniform dT = %g K, linear\n"
            "*INCLUDE, INPUT=../../mesh/%s.inp\n"
            "%s"
            "*NSET, NSET=FRAME_ALL\n"
            "FRAME_X0, FRAME_X1, FRAME_Y0, FRAME_Y1\n"
            "*MATERIAL, NAME=TI64\n"
            "*ELASTIC\n"
            "%.6E, %.4f\n"
            "*EXPANSION, ZERO=20.\n"
            "%.6E\n"
            "*SOLID SECTION, ELSET=PANEL, MATERIAL=TI64\n"
            "*INITIAL CONDITIONS, TYPE=TEMPERATURE\n"
            "NALL, 20.\n"
            "*BOUNDARY\n"
            "%s"
            "FACE_Y0, 2, 2\n"
            "FACE_Y1, 2, 2\n"
            "FRAME_ALL, 3, 3\n"
            "*STEP\n"
            "*STATIC\n"
            "1., 1., 1.E-5, 1.\n"
            "*TEMPERATURE\n"
            "NALL, %.6f\n"
            "*EL PRINT, ELSET=PANEL\n"
            "S\n"
            "*END STEP\n") % (dT, tag, inc, m["E"], m["nu"], m["alpha"], BC[bc], 20.0 + dT)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--levels", type=int, nargs="+", default=[2, 4, 8, 16])
    ap.add_argument("--dT", type=float, default=2.0)
    a = ap.parse_args()
    root = os.path.dirname(HERE)
    m = MATERIAL
    ref = -m["E"] * m["alpha"] * a.dT / (1.0 - m["nu"])
    header = ("level  nx   nz  dx[mm] dz[mm]   nodes  support   sxx_mean[Pa]   sxx_min[Pa]   sxx_max[Pa]  "
              "mean dev[%]  max|dev|[%]  max|szz|[%]  time[s]")
    out = ["reference sigma_x = -E alpha dT/(1-nu) = %.6e Pa (dT = %g K)" % (ref, a.dT), header]
    print("\n".join(out), flush=True)
    rows = []
    for n in a.levels:
        nx, nz = 100 * n, n
        tag = "s_n%d" % n
        mesh = os.path.join(root, "mesh", tag + ".inp")
        subprocess.run([PY, os.path.join(HERE, "make_panel_mesh.py"), tag, "--L", "0.5", "--W", "0.002", "--h", "0.002",
                        "--nx", str(nx), "--ny", "1", "--nz", str(nz)], check=True, stdout=subprocess.DEVNULL)
        subprocess.run([PY, os.path.join(HERE, "frame_constraints.py"), tag, "--faces", "x"], check=True, stdout=subprocess.DEVNULL)
        par = ccxio.read_params(mesh)
        nodes = ccxio.read_nodes(mesh)
        eids, conn = ccxio.read_elements(mesh)
        cx = np.array([np.mean([nodes[int(k)][0] for k in c[:8]]) for c in conn])
        sel = np.sort(eids[np.abs(cx - par["L"] / 2) < 0.05])
        for bc in ("midline", "full", "membrane"):
            job = "slice_%s_%s" % (tag, bc)
            case = os.path.join(root, "simulation", job)
            os.makedirs(case, exist_ok=True)
            with open(os.path.join(case, job + ".inp"), "w") as f:
                f.write(deck(tag, a.dT, bc))
            t0 = time.time()
            with open(os.path.join(case, "ccx.log"), "w") as log:
                rc = subprocess.run([CCX, "-i", job], cwd=case, env=dict(os.environ, OMP_NUM_THREADS="4"),
                                    stdout=log, stderr=subprocess.STDOUT).returncode
            with open(os.path.join(case, "ccx.log")) as f:
                logtxt = f.read()
            if rc != 0 or "Job finished" not in logtxt:
                raise SystemExit("ccx failed for %s (exit %d); see %s" % (job, rc, os.path.join(case, "ccx.log")))
            dat = os.path.join(case, job + ".dat")
            blocks = ccxio.parse_dat(dat)
            st = [b for b in blocks if b["title"].lower().startswith("stresses") and b["set"] == "PANEL"]
            if not st:
                raise SystemExit("no stress block in %s" % dat)
            s = np.array(st[-1]["rows"])
            keep = np.isin(s[:, 0].astype(np.int64), sel)
            sxx, szz = s[keep, 2], s[keep, 4]
            dev_mean = 100.0 * (sxx.mean() - ref) / abs(ref)
            dev_max = 100.0 * np.abs(sxx - ref).max() / abs(ref)
            szz_pct = 100.0 * np.abs(szz).max() / abs(ref)
            dt = time.time() - t0
            os.remove(dat)
            nwarn = logtxt.count("*WARNING")
            rows.append([n, nx, nz, 0.5e3 / nx, 2.0 / nz, len(nodes), bc, sxx.mean(), sxx.min(), sxx.max(), dev_mean, dev_max, szz_pct, dt])
            line = "%5d %4d %4d %7.4f %6.4f %7d  %-8s %13.6e %13.6e %13.6e %11.4e %12.4e %12.4e %7.1f" % tuple(rows[-1])
            if nwarn:
                line += "  (%d warning(s) in ccx.log)" % nwarn
            out.append(line)
            print(line, flush=True)
            with open(os.path.join(root, "results", "slice_study.csv"), "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["level", "nx", "nz", "dx_mm", "dz_mm", "nodes", "support", "sxx_mean_Pa", "sxx_min_Pa", "sxx_max_Pa",
                            "mean_dev_pct", "max_dev_pct", "max_szz_pct", "time_s"])
                w.writerows(rows)
    with open(os.path.join(root, "results", "slice_study.txt"), "w") as f:
        f.write("\n".join(out) + "\n")
    print("written results/slice_study.csv and results/slice_study.txt")


if __name__ == "__main__":
    main()
