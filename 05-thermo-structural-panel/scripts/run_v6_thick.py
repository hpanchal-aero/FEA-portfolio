#!/usr/bin/python3
"""V6d thickness sweep of the linear buckling factor (same model as V6, mesh 20 x 20 x 2, converged in V6c).

Usage: /usr/bin/python3 scripts/run_v6_thick.py
Thicknesses 1, 2, 3 and 5 mm (500 x 500 mm plate). The preload is the mechanical equivalent of a reference rise DT_REF = 0.1 K
(frame_prestrain.py --dT), so every buckling factor of interest is well above 1; factor x DT_REF is the critical temperature
rise. Attempt 1 used 1 K and lost the first mode of the 1 mm plate (record: results/v6d_attempt1_unit_preload_FAIL.txt).
For each thickness the closed form dT_cr = 2 pi^2 h^2 / (12 (1 + nu) alpha b^2) is evaluated with the same h. Reported:
gap = (FE - closed)/closed, the gap divided by h/L, and the log-log slope of |gap| against h between 1 and 5 mm.
Classification fixed before the first attempt and unchanged: slope 0.7-1.3 -> first order in h/L (edge-zone type effect);
1.7-2.3 -> second order (shear-deformation type); otherwise unclassified.
Pass criteria (unchanged): all runs finish without warnings; preload reaction within 0.1 % of E alpha/(1-nu) h W DT_REF;
every thickness in 1-3 mm within 5 % of its closed form. Guard: the run stops if the lowest factor is below 2.
Writes results/v6_thickness.txt and results/v6_thickness.csv. Exit status 1 on FAIL.
"""
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
from run_v6_conv import call, deck, CCX, PY, L

N = 20
NZ = 2
DT_REF = 0.1
THICK_MM = [1, 2, 3, 5]
IN_RANGE_MM = [1, 2, 3]
TOL = 0.05
RF_TOL = 0.001


def closed(h):
    m = MATERIAL
    return 2 * math.pi ** 2 * h ** 2 / (12 * (1 + m["nu"]) * m["alpha"] * L ** 2)


def run_case(root, hmm):
    h = hmm * 1e-3
    tag = "v6d_h%dmm_n%d_z%d" % (hmm, N, NZ)
    call([PY, os.path.join(HERE, "make_panel_mesh.py"), tag, "--h", "%g" % h, "--nx", str(N), "--ny", str(N), "--nz", str(NZ)])
    call([PY, os.path.join(HERE, "frame_constraints.py"), tag, "--faces", "xy"])
    call([PY, os.path.join(HERE, "frame_prestrain.py"), tag, "--dT", "%g" % DT_REF])
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
    raw = read_factors(dat)
    if len(raw) != 6:
        raise SystemExit("%s: expected 6 factors, found %d" % (tag, len(raw)))
    if raw[0] < 2.0:
        raise SystemExit("%s: lowest factor %.4g is below 2; factors below the shift value may be missing, reduce DT_REF" % (tag, raw[0]))
    blocks = [b for b in ccxio.parse_dat(dat) if b["title"].lower().startswith("total force") and b["set"] == "FACE_X0"]
    if not blocks:
        raise SystemExit("%s: no reaction block for FACE_X0" % tag)
    rf = float(blocks[0]["rows"][0][0])
    m = MATERIAL
    rf_exp = m["E"] * m["alpha"] / (1 - m["nu"]) * h * L * DT_REF
    nodes = ccxio.read_nodes(os.path.join(root, "mesh", tag + ".inp"))
    zmax = max(nodes[k][2] for k in nodes)
    if abs(zmax - h) > 1e-9:
        raise SystemExit("%s: mesh thickness %g differs from requested %g" % (tag, zmax, h))
    return {"hmm": hmm, "h": h, "secs": secs, "warn": logtxt.count("*WARNING"), "rf_dev": 100.0 * (abs(rf) - rf_exp) / rf_exp,
            "raw": raw, "f": [v * DT_REF for v in raw], "ref": closed(h)}


def main():
    root = os.path.dirname(HERE)
    out = ["V6d thickness sweep, 0.5 x 0.5 m plate, mesh %d x %d x %d C3D20R, restrained, mechanical-equivalent preload of %g K" % (N, N, NZ, DT_REF),
           "critical rise = buckling factor x %g K" % DT_REF,
           "  h[mm]   h/L[%]  time[s] warn RFdev[%]  factor1   dTcr FE[K]  closed[K]   gap[%]   gap/(h/L)   mode2/mode1  (closed 2.5)"]
    print("\n".join(out), flush=True)
    res = {}
    for hmm in THICK_MM:
        r = run_case(root, hmm)
        res[hmm] = r
        r["gap"] = (r["f"][0] - r["ref"]) / r["ref"]
        r["gapn"] = r["gap"] / (r["h"] / L)
        line = "  %4d   %6.3f  %7.1f %4d %8.4f  %8.4f  %10.5f %10.5f  %+8.4f  %+9.4f   %10.5f" % (
            hmm, 100 * r["h"] / L, r["secs"], r["warn"], r["rf_dev"], r["raw"][0], r["f"][0], r["ref"], 100 * r["gap"], r["gapn"], r["f"][1] / r["f"][0])
        out.append(line)
        print(line, flush=True)
    out.append("")
    g1, g5 = res[1]["gap"], res[5]["gap"]
    if g1 * g5 > 0:
        s = math.log(abs(g5) / abs(g1)) / math.log(5.0)
        cls = "first order in h/L (edge-zone type effect)" if 0.7 <= s <= 1.3 else (
            "second order (shear-deformation type)" if 1.7 <= s <= 2.3 else "unclassified")
        out.append("log-log slope of |gap| against h between 1 and 5 mm: %.2f -> %s" % (s, cls))
    else:
        out.append("gap changes sign between 1 and 5 mm: slope undefined -> unclassified")
    ok = True
    if any(r["warn"] for r in res.values()):
        out.append("FAIL reason: warnings in at least one run")
        ok = False
    if any(abs(r["rf_dev"]) > 100 * RF_TOL for r in res.values()):
        out.append("FAIL reason: preload reaction outside %.1f %%" % (100 * RF_TOL))
        ok = False
    for hmm in IN_RANGE_MM:
        if abs(res[hmm]["gap"]) > TOL:
            out.append("FAIL reason: %d mm outside %.1f %%" % (hmm, 100 * TOL))
            ok = False
    out.append("RESULT: %s" % ("PASS" if ok else "FAIL"))
    print("\n".join(out[len(THICK_MM) + 3:]))
    with open(os.path.join(root, "results", "v6_thickness.txt"), "w") as f:
        f.write("\n".join(out) + "\n")
    with open(os.path.join(root, "results", "v6_thickness.csv"), "w") as f:
        f.write("h_mm,seconds,warnings,rf_dev_percent,dT_ref,factor1,factor2,factor3,factor4,factor5,factor6,closed_form\n")
        for hmm, r in res.items():
            f.write("%d,%.1f,%d,%.5f,%g,%s,%.6f\n" % (hmm, r["secs"], r["warn"], r["rf_dev"], DT_REF, ",".join("%.6f" % v for v in r["raw"]), r["ref"]))
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
