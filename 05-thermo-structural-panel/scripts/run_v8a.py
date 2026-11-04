#!/usr/bin/python3
"""V8a: transient heating of the 500 x 500 x 2 mm panel by the baseline Gaussian pulse (Amendment 4).

Usage: /usr/bin/python3 scripts/run_v8a.py

Cases (mesh, dt): (20 x 20 x 2, 0.5 s), (20 x 20 x 2, 0.25 s), (40 x 40 x 2, 0.5 s); C3D20 (full integration, attempt 2) meshes v8a_n{20,40}_z2.
Flux: top face, Gaussian sigma = L/6 about the centre, q0 = 30 kW/m2, amplitude 0 -> 1 over 5 s, then 1 until 30 s.
Back face and edges insulated. Heat transfer, DIRECT (backward Euler), constant increments, output every 5 s.
Criteria (docs/specification.md, Amendment 5):
  1. volume-averaged temperature rise at t = 5, 10, ..., 30 s within 0.05 % of sum_n P f(t_n) dt / (rho c V)
  2. peak rise at 30 s (20 x 20) changes by < 1 % of the rise between dt = 0.5 s and dt = 0.25 s
  3. symmetry x<->y and x -> L-x within 0.01 K; centre top-bottom difference within 10 % of q0 h/(2k)
The volume average uses exact node weights of the 20-node hexahedron (3 x 3 x 3 Gauss), self-checked.
Writes results/v8a_thermal.txt and results/v8a_thermal_history.csv; exits with status 1 on failure.
"""
import math
import os
import re
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ccxio
from material import MATERIAL

CCX = "/home/harsh/miniconda3/envs/ccx223/bin/ccx"
PY = "/usr/bin/python3"
Q0 = 30000.0
SIGMA = 0.5 / 6.0
XC = YC = 0.25
T0 = 20.0
T_RAMP = 5.0
T_END = 30.0
CASES = [(20, 0.5), (20, 0.25), (40, 0.5)]
TOL_ENERGY = 5.0e-4
TOL_DT = 0.01
TOL_SYM = 0.01
TOL_TT = 0.10

XI = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1], [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1],
               [0, -1, -1], [1, 0, -1], [0, 1, -1], [-1, 0, -1], [0, -1, 1], [1, 0, 1], [0, 1, 1], [-1, 0, 1],
               [-1, -1, 0], [1, -1, 0], [1, 1, 0], [-1, 1, 0]], dtype=float)


def shape_and_derivs(xi, eta, zeta):
    """Shape functions and natural derivatives of the 20-node serendipity hexahedron (Abaqus/CalculiX order)."""
    N = np.zeros(20)
    dN = np.zeros((20, 3))
    for n, (a, b, c) in enumerate(XI):
        if a != 0 and b != 0 and c != 0:
            fa, fb, fc = 1 + a * xi, 1 + b * eta, 1 + c * zeta
            s = a * xi + b * eta + c * zeta - 2
            N[n] = fa * fb * fc * s / 8
            dN[n] = [a * fb * fc * (s + fa) / 8, b * fa * fc * (s + fb) / 8, c * fa * fb * (s + fc) / 8]
        elif a == 0:
            fb, fc = 1 + b * eta, 1 + c * zeta
            N[n] = (1 - xi ** 2) * fb * fc / 4
            dN[n] = [-xi * fb * fc / 2, b * (1 - xi ** 2) * fc / 4, c * (1 - xi ** 2) * fb / 4]
        elif b == 0:
            fa, fc = 1 + a * xi, 1 + c * zeta
            N[n] = (1 - eta ** 2) * fa * fc / 4
            dN[n] = [a * (1 - eta ** 2) * fc / 4, -eta * fa * fc / 2, c * (1 - eta ** 2) * fa / 4]
        else:
            fa, fb = 1 + a * xi, 1 + b * eta
            N[n] = (1 - zeta ** 2) * fa * fb / 4
            dN[n] = [a * (1 - zeta ** 2) * fb / 4, b * (1 - zeta ** 2) * fa / 4, -zeta * fa * fb / 2]
    return N, dN


def node_weights(nodes, eids, conn, L, W, h):
    """w[node] = integral of the node's shape function over the mesh; self-checked against exact integrals."""
    g = np.array([-math.sqrt(0.6), 0.0, math.sqrt(0.6)])
    gw = np.array([5.0 / 9, 8.0 / 9, 5.0 / 9])
    pts = [(g[i], g[j], g[k]) for i in range(3) for j in range(3) for k in range(3)]
    wg = np.array([gw[i] * gw[j] * gw[k] for i in range(3) for j in range(3) for k in range(3)])
    NN = np.zeros((27, 20))
    DN = np.zeros((27, 20, 3))
    for q, (x, y, z) in enumerate(pts):
        NN[q], DN[q] = shape_and_derivs(x, y, z)
    if abs(NN[5].sum() - 1.0) > 1e-12:
        raise SystemExit("shape functions do not sum to 1")
    w = {n: 0.0 for n in nodes}
    for e, c in zip(eids, conn):
        if len(c) < 20:
            raise SystemExit("element %d has %d nodes; V8a needs 20-node elements" % (e, len(c)))
        ids = [int(n) for n in c[:20]]
        X = np.array([nodes[n] for n in ids])
        det = np.linalg.det(np.einsum("gnk,nj->gkj", DN, X))
        if det.min() <= 0.0:
            raise SystemExit("element %d: non-positive Jacobian (node order not standard?)" % e)
        we = np.einsum("g,gi,g->i", wg, NN, det)
        for n, v in zip(ids, we):
            w[n] += v
    ids = list(w)
    wv = np.array([w[n] for n in ids])
    X = np.array([nodes[n] for n in ids])
    V = L * W * h
    checks = [("sum w", wv.sum(), V), ("sum w x", (wv * X[:, 0]).sum(), V * L / 2),
              ("sum w x^2", (wv * X[:, 0] ** 2).sum(), V * L * L / 3), ("sum w z^2", (wv * X[:, 2] ** 2).sum(), V * h * h / 3)]
    for name, got, ref in checks:
        if abs(got - ref) > 1e-10 * abs(ref):
            raise SystemExit("weight self-check failed: %s = %.12e, expected %.12e" % (name, got, ref))
    return w


def deck(mesh, dt, freq):
    m = MATERIAL
    return ("*HEADING\n"
            "Project 05 V8a: Gaussian pulse, q0 = %g W/m2, sigma = L/6, ramp %g s, soak to %g s, dt = %g s\n"
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
            "*NODE PRINT, NSET=COL_CEN, FREQUENCY=%d\n"
            "NT\n"
            "*NODE PRINT, NSET=NALL, FREQUENCY=%d\n"
            "NT\n"
            "*END STEP\n") % (Q0, T_RAMP, T_END, dt, mesh, m["k"], m["rho"], m["cp"], T0, T_RAMP, T_END, dt, T_END, freq, freq)


def discrete_seconds(t_end, dt):
    """sum_n f(t_n) dt with f the ramp evaluated at the end of each increment."""
    s = 0.0
    for i in range(1, int(round(t_end / dt)) + 1):
        s += min(i * dt / T_RAMP, 1.0) * dt
    return s


def key(x, y, z):
    return (int(round(x * 1e7)), int(round(y * 1e7)), int(round(z * 1e7)))


def run_case(root, n, dt, wcache, out, hist):
    m = MATERIAL
    mesh_tag = "v8a_n%d_z2" % n
    mesh_path = os.path.join(root, "mesh", mesh_tag + ".inp")
    if not os.path.exists(mesh_path):
        raise SystemExit("missing mesh %s" % mesh_path)
    par = ccxio.read_params(mesh_path)
    L, W, h = par["L"], par["W"], par["h"]
    nodes = ccxio.read_nodes(mesh_path)
    if n not in wcache:
        eids, conn = ccxio.read_elements(mesh_path)
        wcache[n] = node_weights(nodes, eids, conn, L, W, h)
        out.append("n = %d: node weights self-checks passed (sum, centroid, x^2, z^2 integrals to 1e-10)" % n)
    w = wcache[n]
    job = "v8a2_n%d_dt%s" % (n, ("%g" % dt).replace(".", "p"))
    case = os.path.join(root, "simulation", job)
    os.makedirs(case, exist_ok=True)
    r = subprocess.run([PY, os.path.join(HERE, "make_flux_top.py"), mesh_path, str(Q0), os.path.join(case, "flux.inp"), "PULSE",
                        "--gauss", repr(SIGMA), str(XC), str(YC)], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit("make_flux_top.py failed: %s%s" % (r.stdout, r.stderr))
    mt = re.search(r"(\d+) top faces, total power ([0-9.eE+-]+) W", r.stdout)
    if not mt or int(mt.group(1)) != n * n:
        raise SystemExit("unexpected flux generator output: %s" % r.stdout)
    P = float(mt.group(2))
    freq = int(round(T_RAMP / dt))
    with open(os.path.join(case, job + ".inp"), "w") as f:
        f.write(deck(mesh_tag, dt, freq))
    with open(os.path.join(case, "ccx.log"), "w") as log:
        rc = subprocess.run([CCX, "-i", job], cwd=case, env=dict(os.environ, OMP_NUM_THREADS="4"),
                            stdout=log, stderr=subprocess.STDOUT).returncode
    with open(os.path.join(case, "ccx.log")) as f:
        logtxt = f.read()
    if rc != 0 or "Job finished" not in logtxt:
        raise SystemExit("ccx failed for %s (exit %d); see %s" % (job, rc, os.path.join(case, "ccx.log")))
    nwarn = logtxt.count("*WARNING")
    blocks = [b for b in ccxio.parse_dat(os.path.join(case, job + ".dat")) if b["title"].lower().startswith("temperatures")]
    nall = sorted([b for b in blocks if b["set"] == "NALL"], key=lambda b: b["time"])
    col = sorted([b for b in blocks if b["set"] == "COL_CEN"], key=lambda b: b["time"])
    expect = [T_RAMP * (j + 1) for j in range(int(round(T_END / T_RAMP)))]
    if len(nall) != len(expect) or any(abs(b["time"] - t) > 1e-6 for b, t in zip(nall, expect)):
        raise SystemExit("%s: unexpected NALL output times %s" % (job, [b["time"] for b in nall]))
    C = m["rho"] * m["cp"] * L * W * h
    res = {"P": P, "nwarn": nwarn, "dev": []}
    out.append("case %s: P = %.6f W, rho c V = %.4f J/K, %d warning(s)" % (job, P, C, nwarn))
    out.append("    t[s]   mean rise FE [K]   discrete reference [K]   deviation [%]   continuous-time value [K]")
    for b in nall:
        rows = np.array(b["rows"])
        ids = rows[:, 0].astype(int)
        if len(ids) != len(nodes):
            raise SystemExit("%s: NALL block has %d nodes, mesh has %d" % (job, len(ids), len(nodes)))
        mean = sum(w[i] * (t - T0) for i, t in zip(ids, rows[:, 1])) / (L * W * h)
        ref = P * discrete_seconds(b["time"], dt) / C
        cont = P * (min(b["time"], T_RAMP) ** 2 / (2 * T_RAMP) + max(b["time"] - T_RAMP, 0.0)) / C
        dev = 100.0 * (mean - ref) / ref
        res["dev"].append(abs(dev))
        out.append("  %6.2f %18.5f %24.5f %15.5f %27.5f" % (b["time"], mean, ref, dev, cont))
    res["mean_end"] = mean
    rows = np.array(nall[-1]["rows"])
    ids = rows[:, 0].astype(int)
    T = dict((key(*nodes[i]), t) for i, t in zip(ids, rows[:, 1]))
    sym = 0.0
    for i, t in zip(ids, rows[:, 1]):
        x, y, z = nodes[i]
        for k2 in (key(y, x, z), key(L - x, y, z)):
            if k2 not in T:
                raise SystemExit("%s: symmetric partner node missing" % job)
            sym = max(sym, abs(t - T[k2]))
    imax = int(np.argmax(rows[:, 1]))
    res["sym"] = sym
    res["Tmax"] = rows[imax, 1]
    res["loc"] = nodes[int(ids[imax])]
    for b in col:
        r2 = np.array(b["rows"])
        zz = np.array([nodes[int(k)][2] for k in r2[:, 0]])
        o = np.argsort(zz)
        hist.append("%s,%g,%.6f,%.6f" % (job, b["time"], r2[o[-1], 1], r2[o[0], 1]))
    r2 = np.array(col[-1]["rows"])
    zz = np.array([nodes[int(k)][2] for k in r2[:, 0]])
    o = np.argsort(zz)
    res["tt"] = r2[o[-1], 1] - r2[o[0], 1]
    res["Ttop"] = r2[o[-1], 1]
    out.append("    peak T (all nodes) %.4f degC at (%.4f, %.4f, %.5f) m; centre top %.4f degC; rise %.4f K"
               % (res["Tmax"], res["loc"][0], res["loc"][1], res["loc"][2], res["Ttop"], res["Tmax"] - T0))
    out.append("    symmetry max |dT| %.3e K; centre top - bottom %.4f K (quasi-steady %.4f K)"
               % (sym, res["tt"], Q0 * h / (2 * m["k"])))
    return res


def main():
    root = os.path.dirname(HERE)
    m = MATERIAL
    out = ["V8a: Gaussian pulse thermal verification, q0 = %g W/m2, sigma = %.5f m, ramp %g s, soak to %g s" % (Q0, SIGMA, T_RAMP, T_END),
           "continuous-time total energy per unit power: 27.5 s; scheme adds dt/2 s for the ramp (right-endpoint rule)"]
    hist = ["case,t_s,T_top_centre_C,T_bottom_centre_C"]
    wcache, res = {}, {}
    for n, dt in CASES:
        res[(n, dt)] = run_case(root, n, dt, wcache, out, hist)
        print("\n".join(out[-8:]), flush=True)
    q_tt = Q0 * 0.002 / (2 * m["k"])
    c1 = all(max(r["dev"]) <= 100 * TOL_ENERGY for r in res.values())
    a, b = res[(20, 0.5)], res[(20, 0.25)]
    dtchg = abs(a["Tmax"] - b["Tmax"]) / (b["Tmax"] - T0)
    c2 = dtchg <= TOL_DT
    c3 = all(r["sym"] <= TOL_SYM and abs(r["tt"] - q_tt) <= TOL_TT * q_tt for r in res.values())
    meshchg = abs(res[(40, 0.5)]["Tmax"] - a["Tmax"]) / (a["Tmax"] - T0)
    out.append("")
    out.append("criterion 1 (energy, max deviation over all cases and times): %.5f %% (limit %.3f %%) -> %s"
               % (max(max(r["dev"]) for r in res.values()), 100 * TOL_ENERGY, "pass" if c1 else "FAIL"))
    out.append("criterion 2 (peak rise change dt 0.5 -> 0.25 s): %.4f %% (limit %.1f %%) -> %s" % (100 * dtchg, 100 * TOL_DT, "pass" if c2 else "FAIL"))
    out.append("criterion 3 (symmetry <= %.2f K; top-bottom within %.0f %% of %.4f K) -> %s" % (TOL_SYM, 100 * TOL_TT, q_tt, "pass" if c3 else "FAIL"))
    out.append("reported: peak rise %.3f K (20, dt 0.5), %.3f K (20, dt 0.25), %.3f K (40, dt 0.5); mesh 20 -> 40 changes the peak rise by %.4f %%"
               % (a["Tmax"] - T0, b["Tmax"] - T0, res[(40, 0.5)]["Tmax"] - T0, 100 * meshchg))
    ok = c1 and c2 and c3
    out.append("RESULT: %s" % ("PASS" if ok else "FAIL"))
    print("\n".join(out[-6:]))
    with open(os.path.join(root, "results", "v8a_thermal.txt"), "w") as f:
        f.write("\n".join(out) + "\n")
    with open(os.path.join(root, "results", "v8a_thermal_history.csv"), "w") as f:
        f.write("\n".join(hist) + "\n")
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
