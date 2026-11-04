#!/usr/bin/python3
"""V5: transient heating of a slab with a constant surface flux, insulated back face.

Usage: /usr/bin/python3 scripts/run_v5.py

Case: initial temperature T0 = 20 degC, uniform flux q on the top face (z = h) from t = 0, bottom and side
faces insulated. The field is uniform in x and y, so a small in-plane mesh (8 x 8 elements) gives the same
result as the full plate; the through-thickness mesh is refined (nz = 2, 4, 8; nz = 2 is the plate model).
Backward-Euler time stepping with fixed steps (*HEAT TRANSFER, DIRECT), dt = 0.05, 0.025, 0.0125 s, up to
t = 4 s; the step line allows up to 2000 increments (the CalculiX default of 100 stops longer runs).
Reference (Carslaw and Jaeger; slab 0 < z < h, insulated at z = 0, flux q into z = h), xi = z/h, tau = alpha t/h^2:
  T - T0 = (q h / k) [ tau + (3 xi^2 - 1)/6 - (2/pi^2) sum_n ((-1)^n / n^2) cos(n pi xi) exp(-n^2 pi^2 tau) ]
Compared at t = 0.5, 1, 2, 4 s on all nodes of the centre column.
Tolerance (fixed in docs/specification.md before the run): 2 % of the surface temperature rise, at the finest
time step, for the plate's through-thickness mesh (nz = 2). The refinement in nz and dt is reported.
Writes results/v5_transient.txt; exits with status 1 on failure.
"""
import math
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ccxio
from material import MATERIAL

CCX = "/home/harsh/miniconda3/envs/ccx223/bin/ccx"
PY = "/usr/bin/python3"
Q = 5.0e4
T0 = 20.0
TOTAL = 4.0
CHECK_TIMES = [0.5, 1.0, 2.0, 4.0]
DTS = [0.05, 0.025, 0.0125]
NZS = [2, 4, 8]
TOL = 0.02


def exact(z, t, h):
    m = MATERIAL
    alpha = m["k"] / (m["rho"] * m["cp"])
    tau = alpha * t / h ** 2
    xi = z / h
    n = np.arange(1, 401)[:, None]
    series = (((-1.0) ** n) / n ** 2 * np.cos(n * math.pi * xi[None, :]) * np.exp(-(n ** 2) * math.pi ** 2 * tau)).sum(axis=0)
    return T0 + Q * h / m["k"] * (tau + (3 * xi ** 2 - 1) / 6.0 - 2.0 / math.pi ** 2 * series)


def deck(mesh, dt):
    m = MATERIAL
    return ("*HEADING\n"
            "Project 05 V5: transient slab, constant top-face flux %g W/m2, insulated back face, dt = %g s\n"
            "*INCLUDE, INPUT=../../mesh/%s.inp\n"
            "*MATERIAL, NAME=TI64\n"
            "*CONDUCTIVITY\n%.6g\n"
            "*DENSITY\n%.6g\n"
            "*SPECIFIC HEAT\n%.6g\n"
            "*SOLID SECTION, ELSET=PANEL, MATERIAL=TI64\n"
            "*INITIAL CONDITIONS, TYPE=TEMPERATURE\n"
            "NALL, %g\n"
            "*AMPLITUDE, NAME=FLUXAMP\n"
            "0., 1., 100., 1.\n"
            "*STEP, INC=2000\n"
            "*HEAT TRANSFER, DIRECT\n"
            "%g, %g\n"
            "*INCLUDE, INPUT=flux.inp\n"
            "*NODE PRINT, NSET=COL_CEN\n"
            "NT\n"
            "*END STEP\n") % (Q, dt, mesh, m["k"], m["rho"], m["cp"], T0, dt, TOTAL)


def main():
    root = os.path.dirname(HERE)
    results = {}
    out = ["V5 transient slab, constant top flux q = %g W/m2, T0 = %g degC, h = 0.002 m, in-plane mesh 8 x 8" % (Q, T0),
           "diffusion time h^2/alpha = %.4f s" % (0.002 ** 2 * MATERIAL["rho"] * MATERIAL["cp"] / MATERIAL["k"]),
           "  nz     dt[s]   t[s]   T_top FE   T_top exact   max|T-T_exact| [K]   error [% of surface rise]"]
    print("\n".join(out), flush=True)
    for nz in NZS:
        tagm = "v5_n8_z%d" % nz
        mesh_path = os.path.join(root, "mesh", tagm + ".inp")
        subprocess.run([PY, os.path.join(HERE, "make_panel_mesh.py"), tagm, "--nx", "8", "--ny", "8", "--nz", str(nz)],
                       check=True, stdout=subprocess.DEVNULL)
        h = ccxio.read_params(mesh_path)["h"]
        nodes = ccxio.read_nodes(mesh_path)
        for dt in DTS:
            tag = ("%g" % dt).replace(".", "p")
            job = "v5_z%d_dt%s" % (nz, tag)
            case = os.path.join(root, "simulation", job)
            os.makedirs(case, exist_ok=True)
            r = subprocess.run([PY, os.path.join(HERE, "make_flux_top.py"), mesh_path, str(Q), os.path.join(case, "flux.inp"), "FLUXAMP"],
                               capture_output=True, text=True)
            if r.returncode != 0:
                raise SystemExit("make_flux_top.py failed: %s%s" % (r.stdout, r.stderr))
            if "total power 12500.000000 W" not in r.stdout:
                raise SystemExit("unexpected heater power: %s" % r.stdout)
            with open(os.path.join(case, job + ".inp"), "w") as f:
                f.write(deck(tagm, dt))
            with open(os.path.join(case, "ccx.log"), "w") as log:
                rc = subprocess.run([CCX, "-i", job], cwd=case, env=dict(os.environ, OMP_NUM_THREADS="4"),
                                    stdout=log, stderr=subprocess.STDOUT).returncode
            with open(os.path.join(case, "ccx.log")) as f:
                logtxt = f.read()
            if rc != 0 or "Job finished" not in logtxt:
                raise SystemExit("ccx failed for %s (exit %d); see %s" % (job, rc, os.path.join(case, "ccx.log")))
            nwarn = logtxt.count("*WARNING")
            blocks = [b for b in ccxio.parse_dat(os.path.join(case, job + ".dat"))
                      if b["title"].lower().startswith("temperatures") and b["set"] == "COL_CEN"]
            times = np.array([b["time"] for b in blocks])
            for tc in CHECK_TIMES:
                i = int(np.argmin(np.abs(times - tc)))
                if abs(times[i] - tc) > 1e-6:
                    raise SystemExit("no output at t = %g s for %s (nearest %g)" % (tc, job, times[i]))
                rows = np.array(blocks[i]["rows"])
                z = np.array([nodes[int(k)][2] for k in rows[:, 0]])
                o = np.argsort(z)
                z, T = z[o], rows[o, 1]
                Te = exact(z, tc, h)
                err = 100.0 * np.abs(T - Te).max() / (Te[-1] - T0)
                results[(nz, dt, tc)] = err
                line = "  %2d  %8.4f %5.2f %10.4f %12.4f %20.4e %24.4e" % (nz, dt, tc, T[-1], Te[-1], np.abs(T - Te).max(), err)
                if nwarn:
                    line += "  (%d warning(s))" % nwarn
                out.append(line)
                print(line, flush=True)
    out.append("error [% of surface rise] at t = 0.5 s and t = 1 s (the later times are at round-off):")
    for tc in (0.5, 1.0):
        out.append("  t = %.1f s       dt = %s" % (tc, "  ".join("%8g" % d for d in DTS)))
        for nz in NZS:
            out.append("    nz = %d      %s" % (nz, "  ".join("%8.4f" % results[(nz, d, tc)] for d in DTS)))
    ok = all(results[(NZS[0], DTS[-1], tc)] <= 100 * TOL for tc in CHECK_TIMES)
    out.append("tolerance %.1f %% of the surface rise, nz = %d, dt = %g s" % (100 * TOL, NZS[0], DTS[-1]))
    out.append("RESULT: %s" % ("PASS" if ok else "FAIL"))
    text = "\n".join(out) + "\n"
    print("\n".join(out[-(2 * (len(NZS) + 1) + 3):]))
    with open(os.path.join(root, "results", "v5_transient.txt"), "w") as f:
        f.write(text)
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
