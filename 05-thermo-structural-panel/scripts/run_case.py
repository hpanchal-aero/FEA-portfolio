#!/usr/bin/python3
"""Project 05 parametric case runner (docs/specification.md, Amendments 12 and 13).

Usage: /usr/bin/python3 scripts/run_case.py --h MM --bc rest|free|clamp --imp A_OVER_H [--n 20|40] [--name NAME] [--keep-dat] [--regress]
  e.g. /usr/bin/python3 scripts/run_case.py --h 2 --bc rest --imp 0.1 --regress
  A negative imperfection is best written --imp=-0.1.

Steps: builds the meshes if missing (base C3D20R, thermal C3D20, frame equations, imperfect copy); runs the thermal pulse once per
(h, n) with the temperatures of all nodes at every 0.5 s increment (cached in simulation/p5_thermal_h<h>_n<n>); runs the structural
case as 60 NLGEOM steps of 0.5 s, one thermal frame per step (verified in V8b); extracts the histories; applies G1 to G5.
Edge conditions: rest = frame equations (thickness-mean normal in-plane displacement zero) plus uz = 0 on the mid-thickness edge lines;
free = uz = 0 on the mid-thickness edge lines only, rigid-body modes removed by ux, uy at the centre node and uy at the (L, W/2) node;
clamp = normal in-plane displacement and uz fixed at all nodes of the four edge faces, tangential free.
Writes results/cases/NAME_history.csv and results/cases/NAME.txt. Exit status 1: a check failed; 2: stop rule triggered.
"""
import argparse
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
NFR, DT = v8b.NFRAMES, v8b.DT_T
L_PLATE = 0.5
BOUND = 1.0 / 0.789
STOP = 0.75 * MATERIAL["sigma_y"]
TOL_E, TOL_T, TOL_BAL, TOL_SUP, TOL_REG = 5.0e-4, 0.005, 1.0e-3, 1.0e-4, 1.0e-6


def fmt_tag(x):
    return ("m" if x < 0 else "") + ("%g" % abs(x)).replace(".", "p")


def run(cmd, what):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit("%s failed (exit %d):\n%s%s" % (what, r.returncode, r.stdout[-1500:], r.stderr[-1500:]))
    return r.stdout


def ensure_meshes(root, hmm, n, imp):
    hs = ("%g" % hmm).replace(".", "p")
    base = "p5_h%s_n%d" % (hs, n)
    thermal = base + "_t"
    imp_tag = "%s_a%s" % (base, fmt_tag(imp))
    mdir = os.path.join(root, "mesh")
    common = ["--L", "0.5", "--W", "0.5", "--h", repr(hmm * 1e-3), "--nx", str(n), "--ny", str(n), "--nz", "2"]
    for tag, et in ((base, "C3D20R"), (thermal, "C3D20")):
        if not os.path.exists(os.path.join(mdir, tag + ".inp")):
            run([PY, os.path.join(HERE, "make_panel_mesh.py"), tag] + common + ["--etype", et], "make_panel_mesh " + tag)
    if not os.path.exists(os.path.join(mdir, base + "_frame.inp")):
        run([PY, os.path.join(HERE, "frame_constraints.py"), base, "--faces", "xy"], "frame_constraints " + base)
    if not os.path.exists(os.path.join(mdir, imp_tag + ".inp")):
        run([PY, os.path.join(HERE, "make_imperfect_mesh.py"), base, imp_tag, "--a_over_h=%g" % imp], "make_imperfect_mesh " + imp_tag)
    return hs, base, thermal, imp_tag


def thermal_deck(mesh_tag):
    m = MATERIAL
    return ("*HEADING\nProject 05 case runner: thermal run, Gaussian pulse, dt = %g s\n"
            "*INCLUDE, INPUT=../../mesh/%s.inp\n"
            "*MATERIAL, NAME=TI64\n*CONDUCTIVITY\n%.6g\n*DENSITY\n%.6g\n*SPECIFIC HEAT\n%.6g\n"
            "*SOLID SECTION, ELSET=PANEL, MATERIAL=TI64\n"
            "*INITIAL CONDITIONS, TYPE=TEMPERATURE\nNALL, %g\n"
            "*AMPLITUDE, NAME=PULSE\n0., 0., %g, 1., %g, 1.\n"
            "*STEP, INC=2000\n*HEAT TRANSFER, DIRECT\n%g, %g\n"
            "*INCLUDE, INPUT=flux.inp\n"
            "*NODE PRINT, NSET=NALL\nNT\n"
            "*NODE PRINT, NSET=COL_CEN\nNT\n"
            "*END STEP\n") % (DT, mesh_tag, m["k"], m["rho"], m["cp"], v8a.T0, v8a.T_RAMP, v8a.T_END, DT, v8a.T_END)


def ensure_thermal(root, hs, hmm, n, thermal_tag, keep):
    """Returns (case directory name, dict with the thermal summary)."""
    name = "p5_thermal_h%s_n%d" % (hs, n)
    case = os.path.join(root, "simulation", name)
    info_path = os.path.join(case, "thermal_info.txt")
    if os.path.exists(info_path) and os.path.exists(os.path.join(case, "temps", "frame_%03d.inp" % NFR)):
        info = dict(line.split("=", 1) for line in open(info_path).read().split("\n") if "=" in line)
        return name, dict((k, float(v)) for k, v in info.items())
    mp = os.path.join(root, "mesh", thermal_tag + ".inp")
    os.makedirs(case, exist_ok=True)
    out = run([PY, os.path.join(HERE, "make_flux_top.py"), mp, str(v8a.Q0), os.path.join(case, "flux.inp"), "PULSE",
               "--gauss", repr(v8a.SIGMA), str(v8a.XC), str(v8a.YC)], "make_flux_top")
    P = float(out.split("total power")[1].split("W")[0])
    with open(os.path.join(case, name + ".inp"), "w") as f:
        f.write(thermal_deck(thermal_tag))
    warn = v8b.run_ccx(case, name)
    dat = os.path.join(case, name + ".dat")
    thn = v8b.temp_blocks(dat, "NALL")
    th = v8b.temp_blocks(dat, "COL_CEN")
    expect = [round(k * DT, 6) for k in range(1, NFR + 1)]
    if sorted(thn) != expect or sorted(th) != expect:
        raise SystemExit("thermal output times unexpected: %s" % sorted(thn)[:5])
    par = ccxio.read_params(mp)
    L, W, h = par["L"], par["W"], par["h"]
    nodes = ccxio.read_nodes(mp)
    eids, conn = ccxio.read_elements(mp)
    wts = v8a.node_weights(nodes, eids, conn, L, W, h)
    ids = thn[expect[0]][0]
    wv = np.array([wts[int(i)] for i in ids])
    C = MATERIAL["rho"] * MATERIAL["cp"] * L * W * h
    dev = 0.0
    for t in (5.0, 10.0, 15.0, 20.0, 25.0, 30.0):
        mean = float((wv * (thn[t][1] - v8a.T0)).sum() / (L * W * h))
        ref = P * v8a.discrete_seconds(t, DT) / C
        dev = max(dev, abs(mean - ref) / ref)
    top = int(np.argmax([nodes[int(i)][2] for i in th[expect[0]][0]]))
    peak = [float(thn[t][1].max()) for t in expect]
    topc = [float(th[t][1][top]) for t in expect]
    v8b.write_frames(case, thn)
    with open(os.path.join(case, "thermal_summary.csv"), "w") as f:
        f.write("time_s,T_peak_C,T_top_centre_C\n")
        f.writelines("%.6f,%.6f,%.6f\n" % (t, p, c) for t, p, c in zip(expect, peak, topc))
    info = {"power_W": P, "warnings": float(warn), "energy_dev": dev, "T_peak_end_C": peak[-1], "T_top_centre_end_C": topc[-1]}
    with open(info_path, "w") as f:
        f.writelines("%s=%.10g\n" % kv for kv in info.items())
    if not keep:
        os.remove(dat)
    return name, info


def struct_deck(base, imp_tag, thermal_dir, bc, a, rb_e):
    m = MATERIAL
    s = ["*HEADING\nProject 05 case: edge condition %s, imperfection %g h, Gaussian pulse, frames from %s\n" % (bc, a, thermal_dir),
         "*INCLUDE, INPUT=../../mesh/%s.inp\n" % imp_tag]
    if bc == "rest":
        s.append("*INCLUDE, INPUT=../../mesh/%s_frame.inp\n" % base)
    s.append("*NSET, NSET=FRAME_ALL\nFRAME_X0, FRAME_X1, FRAME_Y0, FRAME_Y1\n")
    s.append("*NSET, NSET=FACE_ALL\nFACE_X0, FACE_X1, FACE_Y0, FACE_Y1\n")
    s.append("*NSET, NSET=ALLN\nNALL\n")
    if bc == "free":
        s.append("*NSET, NSET=RB_E\n%d,\n*NSET, NSET=SUPP\nP_CEN, RB_E\n" % rb_e)
    s.append("*MATERIAL, NAME=TI64\n*ELASTIC\n%.6g, %.6g\n*EXPANSION, ZERO=20.\n%.6g\n"
             "*SOLID SECTION, ELSET=PANEL, MATERIAL=TI64\n"
             "*INITIAL CONDITIONS, TYPE=TEMPERATURE\nNALL, %g\n*BOUNDARY\n" % (m["E"], m["nu"], m["alpha"], v8a.T0))
    if bc == "rest":
        s.append("FRAME_ALL, 3, 3\n")
    elif bc == "free":
        s.append("FRAME_ALL, 3, 3\nP_CEN, 1, 2\nRB_E, 2, 2\n")
    else:
        s.append("FACE_X0, 1, 1\nFACE_X1, 1, 1\nFACE_Y0, 2, 2\nFACE_Y1, 2, 2\nFACE_ALL, 3, 3\n")
    for k in range(1, NFR + 1):
        s.append("*STEP, NLGEOM, INC=2000\n*STATIC\n%g, %g, 1.E-3, %g\n*INCLUDE, INPUT=../%s/temps/frame_%03d.inp\n"
                 "*NODE PRINT, NSET=NALL\nU\n"
                 "*NODE PRINT, NSET=COL_CEN\nNT\n"
                 "*NODE PRINT, NSET=ALLN, TOTALS=ONLY\nRF\n"
                 "*NODE PRINT, NSET=FACE_X0, TOTALS=ONLY\nRF\n"
                 "*NODE PRINT, NSET=FACE_X1, TOTALS=ONLY\nRF\n" % (DT, DT, DT, thermal_dir, k))
        if bc == "free":
            s.append("*NODE PRINT, NSET=SUPP\nRF\n")
        s.append("*EL PRINT, ELSET=PANEL\nS\n*END STEP\n")
    return "".join(s)


def blocks_of(blocks, prefix, setname):
    return sorted([b for b in blocks if b["title"].lower().startswith(prefix) and b["set"] == setname], key=lambda b: b["time"])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h", type=float, required=True, help="thickness in mm")
    ap.add_argument("--bc", choices=["rest", "free", "clamp"], required=True)
    ap.add_argument("--imp", type=float, required=True, help="signed imperfection A/h")
    ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--name")
    ap.add_argument("--keep-dat", action="store_true")
    ap.add_argument("--regress", action="store_true")
    a = ap.parse_args()
    root = os.path.dirname(HERE)
    sim = os.path.join(root, "simulation")
    m = MATERIAL
    hs, base, thermal_tag, imp_tag = ensure_meshes(root, a.h, a.n, a.imp)
    name = a.name or "h%s_%s_a%s_n%d" % (hs, a.bc, fmt_tag(a.imp), a.n)
    out = ["Project 05 case %s: h = %g mm, edge condition %s, imperfection %g h (%s), in-plane mesh %d x %d" %
           (name, a.h, a.bc, a.imp, "toward +z" if a.imp > 0 else "toward -z", a.n, a.n)]
    tdir, tinfo = ensure_thermal(root, hs, a.h, a.n, thermal_tag, a.keep_dat)
    g1 = tinfo["warnings"] == 0 and tinfo["energy_dev"] <= TOL_E
    out.append("thermal: %d warning(s), power %.4f W, energy deviation %.2e (limit %.0e), peak temperature at 30 s %.2f C (rise %.1f K) -> G1 %s"
               % (tinfo["warnings"], tinfo["power_W"], tinfo["energy_dev"], TOL_E, tinfo["T_peak_end_C"], tinfo["T_peak_end_C"] - v8a.T0, "pass" if g1 else "FAIL"))
    print("\n".join(out), flush=True)

    mp_b = os.path.join(root, "mesh", base + ".inp")
    mp_i = os.path.join(root, "mesh", imp_tag + ".inp")
    par = ccxio.read_params(mp_b)
    L, W, h = par["L"], par["W"], par["h"]
    nodes_b, nodes_i = ccxio.read_nodes(mp_b), ccxio.read_nodes(mp_i)
    ids_ref = np.array(sorted(nodes_b))
    dz0 = np.array([nodes_i[i][2] - nodes_b[i][2] for i in ids_ref])
    xy = np.array([[nodes_b[i][0], nodes_b[i][1]] for i in ids_ref])

    def find(x, y, z):
        hit = [i for i in ids_ref if abs(nodes_b[i][0] - x) < 1e-9 and abs(nodes_b[i][1] - y) < 1e-9 and abs(nodes_b[i][2] - z) < 1e-9]
        if len(hit) != 1:
            raise SystemExit("node at (%g, %g, %g) found %d times" % (x, y, z, len(hit)))
        return int(hit[0])
    cen = find(L / 2, W / 2, h / 2)
    ic = int(np.where(ids_ref == cen)[0][0])
    rb_e = find(L, W / 2, h / 2)

    case = os.path.join(sim, name)
    os.makedirs(case, exist_ok=True)
    with open(os.path.join(case, name + ".inp"), "w") as f:
        f.write(struct_deck(base, imp_tag, tdir, a.bc, a.imp, rb_e))
    warn_s = v8b.run_ccx(case, name)
    dat = os.path.join(case, name + ".dat")
    blocks = ccxio.parse_dat(dat)
    ub = blocks_of(blocks, "displacements", "NALL")
    tb = blocks_of(blocks, "temperatures", "COL_CEN")
    ra = blocks_of(blocks, "total force", "ALLN")
    r0 = blocks_of(blocks, "total force", "FACE_X0")
    r1 = blocks_of(blocks, "total force", "FACE_X1")
    sb = blocks_of(blocks, "stresses", "PANEL")
    sup = blocks_of(blocks, "forces", "SUPP") if a.bc == "free" else []
    cnt = [len(ub), len(tb), len(ra), len(r0), len(r1), len(sb)] + ([len(sup)] if a.bc == "free" else [])
    if len(set(cnt)) != 1 or not ub:
        raise SystemExit("unequal block counts: %s (U, NT, RF all, RF x0, RF x1, S%s)" % (cnt, ", SUPP" if a.bc == "free" else ""))
    t = np.array([b["time"] for b in ub])
    for other in (tb, ra, r0, r1, sb) + ((sup,) if sup else ()):
        if np.abs(np.array([b["time"] for b in other]) - t).max() > 1e-6:
            raise SystemExit("output times of the print requests differ")
    expect = [round(k * DT, 6) for k in range(1, NFR + 1)]
    tset = set(round(x, 6) for x in t)
    missing = [x for x in expect if x not in tset]
    extra = sorted(x for x in tset if x not in set(expect))
    g2 = warn_s == 0 and not missing and abs(t[-1] - v8a.T_END) < 1e-6
    out.append("structural: %d warning(s), %d output increments, cutbacks at t = %s s, missing frames %s -> G2 %s" % (warn_s, len(t), extra, missing, "pass" if g2 else "FAIL"))

    thsum = np.genfromtxt(os.path.join(sim, tdir, "thermal_summary.csv"), delimiter=",", names=True)
    wc, wmax, loc, ttop, rx0, rx1, rsum, vmmax = [], [], [], [], [], [], [], []
    tcol = {}
    for b_u, b_t, b_a, b_0, b_1, b_s in zip(ub, tb, ra, r0, r1, sb):
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
        rt = rt[np.argsort(rt[:, 0])]
        tcol[round(b_t["time"], 6)] = (rt[:, 0].astype(int), rt[:, 1])
        rx0.append(b_0["rows"][0][0])
        rx1.append(b_1["rows"][0][0])
        rsum.append(b_a["rows"][0][:3])
        r = np.array(b_s["rows"])
        sxx, syy, szz, sxy, sxz, syz = [r[:, k] for k in range(2, 8)]
        vmmax.append(float(np.sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2) + 3.0 * (sxy ** 2 + sxz ** 2 + syz ** 2)).max()))
    wc, wmax, rx0, rx1, rsum, vmmax = [np.array(v) for v in (wc, wmax, rx0, rx1, rsum, vmmax)]
    ttop = np.array([thsum["T_top_centre_C"][int(round(x / DT)) - 1] if round(x, 6) in set(expect) else np.nan for x in t])

    th_col = v8b.temp_blocks(os.path.join(sim, tdir, tdir + ".dat"), "COL_CEN") if os.path.exists(os.path.join(sim, tdir, tdir + ".dat")) else None
    g3v = float("nan")
    if th_col is not None:
        g3v = max(v8b.maxdiff(tcol[x], th_col[x]) for x in expect if x in tcol)
    else:
        top_ids = None
        g3v = max(abs(float(tcol[x][1].max()) - float(thsum["T_top_centre_C"][int(round(x / DT)) - 1])) for x in expect if x in tcol)
    g3 = g3v <= TOL_T
    out.append("hand-over: structural centre-column temperatures vs thermal history, max |dT| = %.3e K (limit %.3f K) -> G3 %s%s"
               % (g3v, TOL_T, "pass" if g3 else "FAIL", "" if th_col is not None else " (compared with the stored top-centre history)"))

    if a.bc == "free":
        sf = np.array([[abs(v) for row in b["rows"] for v in row[1:4]] for b in sup])
        ref = m["E"] * m["alpha"] * h * L * (tinfo["T_peak_end_C"] - v8a.T0)
        g4v = sf.max() / ref
        g4 = g4v <= TOL_SUP
        out.append("restraint check (free edges): largest rigid-body support reaction %.3e N = %.2e of E alpha h L dT (limit %.0e) -> G4 %s"
                   % (sf.max(), g4v, TOL_SUP, "pass" if g4 else "FAIL"))
    else:
        g4v = np.abs(rx0 + rx1).max() / np.abs(rx0).max()
        g4 = g4v <= TOL_BAL
        out.append("restraint check (face balance): max |Rx(x=0) + Rx(x=L)| / max|Rx(x=0)| = %.3e (limit %.0e) -> G4 %s" % (g4v, TOL_BAL, "pass" if g4 else "FAIL"))
    bound = BOUND * vmmax.max()
    g5 = bound < STOP
    out.append("stress: max von Mises %.2f MPa, bound %.2f MPa (%.1f %% of yield, stop bound %.0f MPa) -> G5 %s"
               % (vmmax.max() / 1e6, bound / 1e6, 100 * bound / m["sigma_y"], STOP / 1e6, "pass" if g5 else "STOP RULE TRIGGERED"))
    out.append("")
    out.append("   t[s]  T_top,c[C]  W_c/h   max|w|/h  at (x, y) [m]      Rx(x=0)[N]   max vM[MPa]")
    for k in range(len(t)):
        if abs(t[k] / 5.0 - round(t[k] / 5.0)) < 1e-6 or round(t[k], 6) in extra:
            out.append("  %6.3f %10.2f %8.4f %9.4f  (%.4f, %.4f) %13.1f %12.3f" % (t[k], ttop[k], wc[k], wmax[k], loc[k][0], loc[k][1], rx0[k], vmmax[k] / 1e6))
    out.append("end of pulse: W_c/h = %.4f (%.3f mm), max|w|/h = %.4f at (%.4f, %.4f), |Rx(x=0)| = %.1f N, T_peak = %.2f C"
               % (wc[-1], wc[-1] * h * 1e3, wmax[-1], loc[-1][0], loc[-1][1], abs(rx0[-1]), tinfo["T_peak_end_C"]))

    gate_ok = True
    if a.regress:
        ref = np.genfromtxt(os.path.join(root, "results", "p05_baseline_history.csv"), delimiter=",", names=True)
        if len(ref["time_s"]) != len(t) or np.abs(ref["time_s"] - t).max() > 1e-6:
            raise SystemExit("regression: output times differ from the baseline history")
        dw = np.abs(wc - ref["W_centre_over_h"]).max() / np.abs(ref["W_centre_over_h"]).max()
        dr = np.abs(rx0 - ref["Rx_face_x0_N"]).max() / np.abs(ref["Rx_face_x0_N"]).max()
        gate_ok = dw <= TOL_REG and dr <= TOL_REG
        out.append("REGRESSION GATE against results/p05_baseline_history.csv: W_centre/h relative difference %.3e, Rx(face x=0) relative difference %.3e (limit %.0e) -> %s"
                   % (dw, dr, TOL_REG, "pass" if gate_ok else "FAIL"))
    ok = g1 and g2 and g3 and g4 and gate_ok
    out.append("RESULT: %s%s" % ("PASS" if ok else "FAIL", "" if g5 else " (stop rule triggered, ask before continuing)"))
    print("\n".join(out[2:]))
    os.makedirs(os.path.join(root, "results", "cases"), exist_ok=True)
    with open(os.path.join(root, "results", "cases", name + ".txt"), "w") as f:
        f.write("\n".join(out) + "\n")
    with open(os.path.join(root, "results", "cases", name + "_history.csv"), "w") as f:
        f.write("time_s,T_top_centre_C,W_centre_over_h,wmax_over_h,wmax_x_m,wmax_y_m,Rx_face_x0_N,Rx_face_x1_N,max_vM_MPa,bound_MPa\n")
        for k in range(len(t)):
            f.write("%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.4f,%.4f,%.6f,%.6f\n" % (t[k], ttop[k], wc[k], wmax[k], loc[k][0], loc[k][1], rx0[k], rx1[k], vmmax[k] / 1e6, BOUND * vmmax[k] / 1e6))
    if not a.keep_dat:
        os.remove(dat)
    if not g5:
        sys.exit(2)
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
