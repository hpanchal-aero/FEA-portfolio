#!/usr/bin/python3
"""Project 05 baseline: Gaussian heating pulse on the 2 mm Ti-6Al-4V panel, simply supported frame, imperfection 0.1 h.

Usage: /usr/bin/python3 scripts/run_baseline.py

1. Thermal run (C3D20 mesh v8a_n20_z2, verified in V8a): temperatures of all nodes at every 0.5 s increment.
2. Structural run (C3D20R imperfect mesh v7_n20_z2_a0p1 plus the frame equations of v6c_n20_z2, as in V7): 60 NLGEOM steps of
   0.5 s, each ramping from one thermal frame to the next (verified in V8b); uz = 0 on the mid-thickness frame lines.
   Printed every increment: displacements of all nodes, centre-column temperatures, total reaction of all nodes (set ALLN)
   and of face x = 0, element stresses (parsed in the next stage).
Acceptance (docs/specification.md, Amendment 10), code correctness only:
  B1 run reaches t = 30 s, every frame time present, zero warnings (cutbacks reported)
  B2 centre-column temperatures equal the thermal history at the 60 frames within 0.005 K
  B3 total reaction (x, y, z) <= 1e-3 of the largest |Rx| on face x = 0, at every printed increment
Writes results/p05_baseline.txt and results/p05_baseline_history.csv; exits with status 1 if B1 to B3 fail.
"""
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ccxio
import run_v8a as v8a
import run_v8b as v8b
from material import MATERIAL

CCX, PY = v8a.CCX, v8a.PY
MESH_T = v8b.MESH_T
BASE = "v6c_n20_z2"
MESH_S = "v7_n20_z2_a0p1"
A_OVER_H = 0.1
SIM_T, SIM_S = "p05_base_thermal", "p05_base_struct"
TOL_T = 0.005
TOL_EQ = 1.0e-3


def struct_deck():
    m = MATERIAL
    head = ("*HEADING\n"
            "Project 05 baseline: Gaussian heating pulse, 2 mm Ti-6Al-4V panel, simply supported frame, imperfection %g h\n"
            "*INCLUDE, INPUT=../../mesh/%s.inp\n"
            "*INCLUDE, INPUT=../../mesh/%s_frame.inp\n"
            "*NSET, NSET=FRAME_ALL\n"
            "FRAME_X0, FRAME_X1, FRAME_Y0, FRAME_Y1\n"
            "*NSET, NSET=ALLN\n"
            "NALL\n"
            "*MATERIAL, NAME=TI64\n"
            "*ELASTIC\n%.6g, %.6g\n"
            "*EXPANSION, ZERO=20.\n%.6g\n"
            "*SOLID SECTION, ELSET=PANEL, MATERIAL=TI64\n"
            "*INITIAL CONDITIONS, TYPE=TEMPERATURE\n"
            "NALL, %g\n"
            "*BOUNDARY\n"
            "FRAME_ALL, 3, 3\n") % (A_OVER_H, MESH_S, BASE, m["E"], m["nu"], m["alpha"], v8a.T0)
    steps = []
    for k in range(1, v8b.NFRAMES + 1):
        steps.append("*STEP, NLGEOM, INC=2000\n"
                     "*STATIC\n"
                     "%g, %g, 1.E-3, %g\n"
                     "*INCLUDE, INPUT=temps/frame_%03d.inp\n"
                     "*NODE PRINT, NSET=NALL\n"
                     "U\n"
                     "*NODE PRINT, NSET=COL_CEN\n"
                     "NT\n"
                     "*NODE PRINT, NSET=ALLN, TOTALS=ONLY\n"
                     "RF\n"
                     "*NODE PRINT, NSET=FACE_X0, TOTALS=ONLY\n"
                     "RF\n"
                     "*EL PRINT, ELSET=PANEL\n"
                     "S\n"
                     "*END STEP\n" % (v8b.DT_T, v8b.DT_T, v8b.DT_T, k))
    return head + "".join(steps)


def thermal_run(root, sim):
    mp_t = os.path.join(root, "mesh", MESH_T + ".inp")
    case = os.path.join(sim, SIM_T)
    os.makedirs(case, exist_ok=True)
    r = subprocess.run([PY, os.path.join(HERE, "make_flux_top.py"), mp_t, str(v8a.Q0), os.path.join(case, "flux.inp"), "PULSE",
                        "--gauss", repr(v8a.SIGMA), str(v8a.XC), str(v8a.YC)], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit("make_flux_top.py failed: %s%s" % (r.stdout, r.stderr))
    with open(os.path.join(case, SIM_T + ".inp"), "w") as f:
        f.write(v8b.thermal_deck())
    warn = v8b.run_ccx(case, SIM_T)
    dat = os.path.join(case, SIM_T + ".dat")
    thn = v8b.temp_blocks(dat, "NALL")
    th = v8b.temp_blocks(dat, "COL_CEN")
    expect = [round(k * v8b.DT_T, 6) for k in range(1, v8b.NFRAMES + 1)]
    if sorted(thn) != expect or sorted(th) != expect:
        raise SystemExit("thermal output times unexpected: %s" % sorted(thn))
    return r.stdout.strip(), warn, thn, th, expect


def blocks_of(blocks, prefix, setname):
    return sorted([b for b in blocks if b["title"].lower().startswith(prefix) and b["set"] == setname], key=lambda b: b["time"])


def main():
    root = os.path.dirname(HERE)
    sim = os.path.join(root, "simulation")
    m = MATERIAL
    out = ["Project 05 baseline: Gaussian pulse (q0 = %g W/m2, sigma = L/6, ramp %g s, soak to %g s), 2 mm Ti-6Al-4V, imperfection %g h toward +z"
           % (v8a.Q0, v8a.T_RAMP, v8a.T_END, A_OVER_H)]
    pw, warn_t, thn, th, expect = thermal_run(root, sim)
    out.append(pw)
    out.append("thermal run: %d warning(s), %d frames" % (warn_t, len(thn)))
    print("\n".join(out), flush=True)

    case = os.path.join(sim, SIM_S)
    os.makedirs(case, exist_ok=True)
    v8b.write_frames(case, thn)
    with open(os.path.join(case, SIM_S + ".inp"), "w") as f:
        f.write(struct_deck())
    warn_s = v8b.run_ccx(case, SIM_S)
    blocks = ccxio.parse_dat(os.path.join(case, SIM_S + ".dat"))
    ub = blocks_of(blocks, "displacements", "NALL")
    tb = blocks_of(blocks, "temperatures", "COL_CEN")
    ra = blocks_of(blocks, "total force", "ALLN")
    rx = blocks_of(blocks, "total force", "FACE_X0")
    if not (len(ub) == len(tb) == len(ra) == len(rx)) or not ub:
        raise SystemExit("unequal block counts: U %d, NT %d, RF all %d, RF x0 %d" % (len(ub), len(tb), len(ra), len(rx)))
    t = np.array([b["time"] for b in ub])
    for other in (tb, ra, rx):
        if np.abs(np.array([b["time"] for b in other]) - t).max() > 1e-6:
            raise SystemExit("output times of the print requests differ")
    tset = set(round(x, 6) for x in t)
    missing = [x for x in expect if x not in tset]
    extra = sorted(x for x in tset if x not in set(expect))
    out.append("structural run: %d warning(s), %d output increments, extra increments from cutbacks at t = %s s, missing frames %s"
               % (warn_s, len(t), extra, missing))

    mp_i = os.path.join(root, "mesh", MESH_S + ".inp")
    mp_f = os.path.join(root, "mesh", BASE + ".inp")
    par = ccxio.read_params(mp_f)
    L, W, h = par["L"], par["W"], par["h"]
    nodes_i, nodes_f = ccxio.read_nodes(mp_i), ccxio.read_nodes(mp_f)
    ids_ref = np.array(sorted(nodes_f))
    dz0 = np.array([nodes_i[i][2] - nodes_f[i][2] for i in ids_ref])
    xy = np.array([[nodes_f[i][0], nodes_f[i][1]] for i in ids_ref])
    cen = [i for i in ids_ref if abs(nodes_f[i][0] - L / 2) < 1e-9 and abs(nodes_f[i][1] - W / 2) < 1e-9 and abs(nodes_f[i][2] - h / 2) < 1e-9]
    if len(cen) != 1:
        raise SystemExit("centre mid-plane node found %d times" % len(cen))
    ic = int(np.where(ids_ref == cen[0])[0][0])
    col_ids = th[expect[0]][0]
    top_col = int(np.argmax([nodes_f[int(i)][2] for i in col_ids]))

    wc, wmax, loc, ttop, fx0, rsum = [], [], [], [], [], []
    for b_u, b_t, b_a, b_x in zip(ub, tb, ra, rx):
        rows = np.array(b_u["rows"])
        o = np.argsort(rows[:, 0])
        if not np.array_equal(rows[o, 0].astype(int), ids_ref):
            raise SystemExit("node ids of the displacement block differ from the mesh")
        w = dz0 + rows[o, 3]
        j = int(np.argmax(np.abs(w)))
        wc.append(w[ic] / h)
        wmax.append(w[j] / h)
        loc.append(xy[j])
        rt = np.array(b_t["rows"])
        ttop.append(rt[np.argsort(rt[:, 0])][top_col, 1])
        fx0.append(b_x["rows"][0][0])
        rsum.append([b_a["rows"][0][0], b_a["rows"][0][1], b_a["rows"][0][2]])
    wc, wmax, ttop, fx0, rsum = np.array(wc), np.array(wmax), np.array(ttop), np.array(fx0), np.array(rsum)

    b2 = max(v8b.maxdiff(b_t_ids, (th[round(tt, 6)][0], th[round(tt, 6)][1]))
             for tt, b_t_ids in ((tt, (np.array(bt["rows"])[np.argsort(np.array(bt["rows"])[:, 0]), 0].astype(int),
                                       np.array(bt["rows"])[np.argsort(np.array(bt["rows"])[:, 0]), 1]))
                                 for tt, bt in zip(t, tb)) if round(tt, 6) in th)
    fref = np.abs(fx0).max()
    b3 = np.abs(rsum).max() / fref
    ncr = m["E"] * h * m["alpha"] * 2.30101 / (1 - m["nu"])
    out.append("")
    out.append("   t[s]  T_top,centre[C]  W_centre/h  max|w|/h    at (x, y) [m]       Rx(face x=0) [N]  N/Ncr(1)")
    for k in range(len(t)):
        if abs(t[k] / 2.5 - round(t[k] / 2.5)) < 1e-6 or round(t[k], 6) in extra:
            out.append("  %6.3f %14.3f %11.4f %9.4f   (%.4f, %.4f) %18.1f %9.3f"
                       % (t[k], ttop[k], wc[k], wmax[k], loc[k][0], loc[k][1], fx0[k], abs(fx0[k]) / L / ncr))
    out.append("(N/Ncr(1) = |Rx|/L divided by E h alpha dT_cr/(1 - nu) with dT_cr = 2.30101 K, the uniform-heating critical value; for scale only)")
    b1 = warn_s == 0 and warn_t == 0 and not missing and abs(t[-1] - v8a.T_END) < 1e-6
    out.append("B1 (completes, all frames present, zero warnings): %s" % ("pass" if b1 else "FAIL"))
    out.append("B2 (structural centre-column temperatures vs thermal history): max |dT| = %.3e K (limit %.3f K) -> %s" % (b2, TOL_T, "pass" if b2 <= TOL_T else "FAIL"))
    out.append("B3 (global equilibrium): max |total reaction| / max |Rx(face x=0)| = %.3e (limit %.0e) -> %s; largest |Rx(face x=0)| = %.1f N"
               % (b3, TOL_EQ, "pass" if b3 <= TOL_EQ else "FAIL", fref))
    ok = b1 and b2 <= TOL_T and b3 <= TOL_EQ
    out.append("RESULT: %s" % ("PASS" if ok else "FAIL"))
    print("\n".join(out[3:]))
    with open(os.path.join(root, "results", "p05_baseline.txt"), "w") as f:
        f.write("\n".join(out) + "\n")
    with open(os.path.join(root, "results", "p05_baseline_history.csv"), "w") as f:
        f.write("time_s,T_top_centre_C,W_centre_over_h,wmax_over_h,wmax_x_m,wmax_y_m,Rx_face_x0_N,sumRx_N,sumRy_N,sumRz_N\n")
        for k in range(len(t)):
            f.write("%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.4f,%.6e,%.6e,%.6e\n" % (t[k], ttop[k], wc[k], wmax[k], loc[k][0], loc[k][1], fx0[k], rsum[k][0], rsum[k][1], rsum[k][2]))
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
