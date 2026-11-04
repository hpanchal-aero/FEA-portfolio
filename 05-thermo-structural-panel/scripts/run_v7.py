#!/usr/bin/python3
"""V7: geometrically nonlinear thermal post-buckling of the restrained, imperfect plate versus the single-mode references.

Usage: /usr/bin/python3 scripts/run_v7.py [--n 20|40] [A_OVER_H ...]      (default n = 20, A/h = 0.1)
Meshes: mesh/v7_n<n>_z2_a<tag>.inp (imperfect, scripts/make_imperfect_mesh.py) and mesh/v6c_n<n>_z2_frame.inp must exist.
One *STATIC, NLGEOM step ramps a uniform *TEMPERATURE rise from 0 to tau_max = 2.5 dT_cr (dT_cr = 2.30101 K, closed form),
increments of 0.02 step time (= 0.05 in tau), automatic cutback allowed. Supports as in V3c/V6 (membrane restraint equations of
the base mesh, uz = 0 on the mid-thickness edge lines). Printed every increment: displacement of the centre mid-plane node
(set P_CEN) and the total reaction on face x = 0.
Compared: W = A + uz(P_CEN) (total deflection) against the corrected reference of Amendment 3 (primary) and against the
Amendment 2 reference (printed for the record only). Acceptance (docs/specification.md, Amendment 3):
  (1) the run reaches tau = 2.5 without solver errors or warnings;
  (2) blind (A = 0.01 h and 0.5 h): |W_FE/W_ref - 1| <= 5 % at tau = 1.2, 1.5, 2.0;
  (3) post hoc (A = 0.1 h, derived after the result was known): the same criterion.
  (4) mesh check, see scripts/check_v7_mesh.py.
Diagnostic, no criterion: membrane force N = |RF_x|/L in units of N_cr = E h alpha dT_cr/(1 - nu) against the corrected relation.
Writes results/v7_[n40_]a<tag>.txt and results/v7_[n40_]a<tag>_history.csv (the n = 20 names carry no n; the failed first
evaluation of A = 0.1 h against the Amendment 2 reference is kept in results/v7_a0p1_amendment2_FAIL.txt).
Exit status 1 if any criterion fails.
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
from material import MATERIAL
from ref_postbuckling import reference_wh, reference_wh_corrected, reference_n_corrected, dT_cr, H, L, TAUS

CCX = "/home/harsh/miniconda3/envs/ccx223/bin/ccx"
TAU_MAX = 2.5
DINC = 0.02
CRIT_TAUS = [1.2, 1.5, 2.0]
TOL = 0.05


def tagof(a):
    return "a" + ("%g" % a).replace(".", "p")


def deck(mesh_tag, base, tmax_c):
    m = MATERIAL
    return ("*HEADING\n"
            "Project 05 V7: geometrically nonlinear thermal post-buckling, imperfect restrained plate, mesh %s\n"
            "*INCLUDE, INPUT=../../mesh/%s.inp\n"
            "*INCLUDE, INPUT=../../mesh/%s_frame.inp\n"
            "*NSET, NSET=FRAME_ALL\n"
            "FRAME_X0, FRAME_X1, FRAME_Y0, FRAME_Y1\n"
            "*MATERIAL, NAME=TI64\n"
            "*ELASTIC\n"
            "%.6g, %.6g\n"
            "*EXPANSION, ZERO=20.\n"
            "%.6g\n"
            "*SOLID SECTION, ELSET=PANEL, MATERIAL=TI64\n"
            "*INITIAL CONDITIONS, TYPE=TEMPERATURE\n"
            "NALL, 20.\n"
            "*BOUNDARY\n"
            "FRAME_ALL, 3, 3\n"
            "*STEP, NLGEOM, INC=2000\n"
            "*STATIC\n"
            "%g, 1., 1.E-6, %g\n"
            "*TEMPERATURE\n"
            "NALL, %.9g\n"
            "*NODE PRINT, NSET=P_CEN\n"
            "U\n"
            "*NODE PRINT, NSET=FACE_X0, TOTALS=ONLY\n"
            "RF\n"
            "*END STEP\n") % (mesh_tag, mesh_tag, base, m["E"], m["nu"], m["alpha"], DINC, DINC, tmax_c)


def at(t_arr, y_arr, t):
    i = int(np.searchsorted(t_arr, t - 1e-9))
    if i < len(t_arr) and abs(t_arr[i] - t) < 1e-9:
        return y_arr[i], 0.0
    if i == 0 or i >= len(t_arr):
        raise SystemExit("time %g outside the printed range" % t)
    t0, t1 = t_arr[i - 1], t_arr[i]
    w = (t - t0) / (t1 - t0)
    return y_arr[i - 1] * (1 - w) + y_arr[i] * w, t1 - t0


def run_one(root, a, n):
    base = "v6c_n%d_z2" % n
    tag = "v7_n%d_z2_%s" % (n, tagof(a))
    mesh = os.path.join(root, "mesh", tag + ".inp")
    for p in (mesh, os.path.join(root, "mesh", base + "_frame.inp")):
        if not os.path.exists(p):
            raise SystemExit("missing %s" % p)
    d = dT_cr()
    name = ("v7_" if n == 20 else "v7_n%d_" % n) + tagof(a)
    case = os.path.join(root, "simulation", name)
    os.makedirs(case, exist_ok=True)
    with open(os.path.join(case, name + ".inp"), "w") as f:
        f.write(deck(tag, base, 20.0 + TAU_MAX * d))
    t0 = time.time()
    with open(os.path.join(case, "ccx.log"), "w") as log:
        rc = subprocess.run([CCX, "-i", name], cwd=case, env=dict(os.environ, OMP_NUM_THREADS="4"),
                            stdout=log, stderr=subprocess.STDOUT).returncode
    secs = time.time() - t0
    with open(os.path.join(case, "ccx.log")) as f:
        logtxt = f.read()
    nwarn = logtxt.count("*WARNING")
    finished = rc == 0 and "Job finished" in logtxt
    dat = os.path.join(case, name + ".dat")
    blocks = ccxio.parse_dat(dat) if os.path.exists(dat) else []
    ub = [b for b in blocks if b["title"].lower().startswith("displacements") and b["set"] == "P_CEN"]
    fb = [b for b in blocks if b["title"].lower().startswith("total force") and b["set"] == "FACE_X0"]
    if len(ub) != len(fb) or not ub:
        raise SystemExit("%s: %d displacement and %d reaction blocks (exit %d, finished %s); see %s" %
                         (name, len(ub), len(fb), rc, finished, os.path.join(case, "ccx.log")))
    t = np.array([b["time"] for b in ub])
    uz = np.array([b["rows"][0][3] for b in ub])
    rx = np.array([abs(b["rows"][0][0]) for b in fb])
    o = np.argsort(t)
    t, uz, rx = t[o], uz[o], rx[o]
    A = a * H
    W = (A + uz) / H
    m = MATERIAL
    ncr = m["E"] * H * m["alpha"] * d / (1 - m["nu"])
    N = rx / L / ncr
    reached = abs(t[-1] - 1.0) < 1e-9
    kind = "post hoc" if abs(a - 0.1) < 1e-12 else "blind"
    out = ["V7 nonlinear thermal post-buckling, imperfection A = %g h = %.4f mm (%s), mesh %s, dT_cr(closed form) = %.5f K" % (a, A * 1e3, kind, tag, d),
           "run: exit %d, finished %s, warnings %d, %d printed increments, last time %.6f (target 1), %.0f s" % (rc, finished, nwarn, len(t), t[-1], secs),
           "  tau    dT[K]   W_FE/h  W_ref3/h  err3[%]  W_ref2/h  err2[%]  N_FE/Ncr  N_ref3/Ncr   note",
           "  (ref3 = corrected reference, Amendment 3, primary; ref2 = Amendment 2 reference, for the record only)"]
    errs = {}
    for tau in TAUS:
        w, gap = at(t, W, tau / TAU_MAX)
        nn, _ = at(t, N, tau / TAU_MAX)
        r3 = reference_wh_corrected(a, tau)
        r2 = reference_wh(a, tau)
        n3 = reference_n_corrected(a, tau)
        e3 = 100.0 * (w / r3 - 1.0)
        e2 = 100.0 * (w / r2 - 1.0)
        errs[tau] = e3
        note = "" if gap == 0 else "interpolated over %.3f" % gap
        out.append("  %.2f  %7.4f  %7.4f  %8.4f  %+7.3f  %8.4f  %+7.3f  %8.4f  %9.4f   %s" % (tau, tau * d, w, r3, e3, r2, e2, nn, n3, note))
    for line in out:
        print(line)
    with open(os.path.join(root, "results", name + "_history.csv"), "w") as f:
        f.write("time,tau,dT_K,W_over_h,N_over_Ncr\n")
        for i in range(len(t)):
            f.write("%.6f,%.6f,%.6f,%.6f,%.6f\n" % (t[i], t[i] * TAU_MAX, t[i] * TAU_MAX * d, W[i], N[i]))
    ok_run = finished and reached and nwarn == 0
    ok_w = all(abs(errs[tau]) <= 100 * TOL for tau in CRIT_TAUS)
    tail = []
    if not ok_run:
        tail.append("  criterion 1 not met (finished %s, reached tau = 2.5: %s, warnings %d)" % (finished, reached, nwarn))
    for tau in CRIT_TAUS:
        if abs(errs[tau]) > 100 * TOL:
            tail.append("  criterion on W not met at tau = %.1f (%+.2f %%, tolerance %.0f %%)" % (tau, errs[tau], 100 * TOL))
    passed = ok_run and ok_w
    out.append("RESULT (%s): %s" % (kind, "PASS" if passed else "FAIL"))
    out.extend(tail)
    print("\n".join(out[len(TAUS) + 4:]))
    with open(os.path.join(root, "results", name + ".txt"), "w") as f:
        f.write("\n".join(out) + "\n")
    return kind, passed, errs


def main():
    root = os.path.dirname(HERE)
    args = sys.argv[1:]
    n = 20
    if args and args[0] == "--n":
        n = int(args[1])
        args = args[2:]
    if n not in (20, 40):
        raise SystemExit("--n must be 20 or 40")
    amps = [float(x) for x in args] or [0.1]
    summary = []
    for a in amps:
        kind, passed, errs = run_one(root, a, n)
        summary.append((a, kind, passed, errs))
    print("\nSUMMARY, in-plane mesh %d x %d (error of W_FE against the corrected reference at tau = 1.2, 1.5, 2.0; tolerance %.0f %%)" % (n, n, 100 * TOL))
    for a, kind, passed, errs in summary:
        print("  A = %-5g h  %-8s  %s   %s" % (a, kind, "PASS" if passed else "FAIL",
                                              "  ".join("%+.2f %%" % errs[tau] for tau in CRIT_TAUS)))
    if not all(s[2] for s in summary):
        sys.exit(1)


if __name__ == "__main__":
    main()
