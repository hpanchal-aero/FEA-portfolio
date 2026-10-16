#!/usr/bin/env python3
"""Constant-strain patch test of a CalculiX S4 shell with a rotated in-plane material axis.

Flat 4x4 mesh of 0.25 m squares in the XY plane. Orientation as in the uCRM skins:
material 1 along a=(cos55,sin55,0), b=(-sin55,cos55,0). Six strain states are imposed on all
nodes; the work of the reaction forces, W = sum(RF . u) = t*A * e^T Qbar e, gives the six
independent entries of the rotated plane-stress stiffness matrix.
"""
import glob
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

PROJ = Path(__file__).resolve().parent.parent
WORK = PROJ / "simulation" / "orientation_test"

E1, E2, NU12, G12 = 70e9, 45.5e9, 0.3, 16.9e9
E3, NU13, NU23, G13, G23 = 45.5e9, 0.3, 0.3, 16.08e9, 9.72e9
T, SIDE, N = 0.01, 1.0, 4
THETA = np.radians(55.0)
S = 1.0e-6


def find_ccx():
    env = os.environ.get("CCX")
    if env:
        return env
    cands = [c for c in sorted(glob.glob(str(Path.home() / "miniconda3/envs/ccx223/bin/ccx*")))
             if Path(c).is_file() and os.access(c, os.X_OK) and not c.endswith(".so")]
    if not cands:
        sys.exit("ABORT: no ccx binary found in the ccx223 env; set CCX=/path/to/ccx")
    return cands[0]


def analytic():
    nu21 = NU12 * E2 / E1
    d = 1.0 - NU12 * nu21
    Q11, Q22, Q12, Q66 = E1 / d, E2 / d, NU12 * E2 / d, G12
    m, n = np.cos(THETA), np.sin(THETA)
    Qb11 = Q11 * m**4 + 2 * (Q12 + 2 * Q66) * m**2 * n**2 + Q22 * n**4
    Qb22 = Q11 * n**4 + 2 * (Q12 + 2 * Q66) * m**2 * n**2 + Q22 * m**4
    Qb12 = (Q11 + Q22 - 4 * Q66) * m**2 * n**2 + Q12 * (m**4 + n**4)
    Qb66 = (Q11 + Q22 - 2 * Q12 - 2 * Q66) * m**2 * n**2 + Q66 * (m**4 + n**4)
    Qb16 = (Q11 - Q12 - 2 * Q66) * m**3 * n - (Q22 - Q12 - 2 * Q66) * m * n**3
    Qb26 = (Q11 - Q12 - 2 * Q66) * m * n**3 - (Q22 - Q12 - 2 * Q66) * m**3 * n
    return dict(Q11=Qb11, Q22=Qb22, Q66=Qb66, Q12=Qb12, Q16=Qb16, Q26=Qb26)


def write_case(name, e):
    e1, e2, e3 = e
    h = SIDE / N
    ids, coords, nid = {}, [], 1
    for j in range(N + 1):
        for i in range(N + 1):
            ids[(i, j)] = nid
            coords.append((nid, i * h, j * h))
            nid += 1
    lines = ["*NODE"] + [f"{k}, {x:.10f}, {y:.10f}, 0." for k, x, y in coords]
    lines += ["*ELEMENT, TYPE=S4, ELSET=EALL"]
    k = 0
    for j in range(N):
        for i in range(N):
            k += 1
            lines.append(f"{k}, {ids[(i, j)]}, {ids[(i + 1, j)]}, {ids[(i + 1, j + 1)]}, {ids[(i, j + 1)]}")
    lines += ["*NSET, NSET=ALLN, GENERATE", f"1, {nid - 1}, 1"]
    ca, sa = np.cos(THETA), np.sin(THETA)
    lines += ["*ORIENTATION, NAME=OR1, SYSTEM=RECTANGULAR", f"{ca:.10f}, {sa:.10f}, 0., {-sa:.10f}, {ca:.10f}, 0.",
              "*MATERIAL, NAME=MAT", "*ELASTIC, TYPE=ENGINEERING CONSTANTS",
              f"{E1:.6e}, {E2:.6e}, {E3:.6e}, {NU12}, {NU13}, {NU23}, {G12:.6e}, {G13:.6e}",
              f"{G23:.6e}, 0.",
              "*SHELL SECTION, ELSET=EALL, MATERIAL=MAT, ORIENTATION=OR1", f"{T}",
              "*STEP", "*STATIC", "*BOUNDARY"]
    disp = {}
    for k2, x, y in coords:
        ux = S * (e1 * x + 0.5 * e3 * y)
        uy = S * (e2 * y + 0.5 * e3 * x)
        disp[k2] = (ux, uy)
        lines += [f"{k2}, 1, 1, {ux:.10e}", f"{k2}, 2, 2, {uy:.10e}", f"{k2}, 3, 6, 0."]
    lines += ["*NODE PRINT, NSET=ALLN", "RF", "*END STEP"]
    (WORK / f"{name}.inp").write_text("\n".join(lines) + "\n")
    return disp


def read_rf(dat):
    rf, grab = {}, False
    for ln in dat.read_text().splitlines():
        if "forces" in ln:
            grab = True
            continue
        if grab:
            tok = ln.split()
            if len(tok) == 4:
                try:
                    rf[int(tok[0])] = (float(tok[1]), float(tok[2]))
                    continue
                except ValueError:
                    pass
            if rf and not ln.strip():
                break
    return rf


def main():
    ccx = find_ccx()
    WORK.mkdir(parents=True, exist_ok=True)
    print("ccx:", ccx, "|", subprocess.run([ccx, "-v"], capture_output=True, text=True).stdout.strip()[:40])
    cases = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1)]
    A_mat, rhs = [], []
    for k, e in enumerate(cases):
        name = f"patch{k}"
        disp = write_case(name, e)
        r = subprocess.run([ccx, name], cwd=WORK, capture_output=True, text=True, timeout=600)
        dat = WORK / f"{name}.dat"
        if r.returncode != 0 or not dat.exists():
            print(r.stdout[-1500:], r.stderr[-1500:])
            sys.exit(f"ABORT: ccx failed for {name}")
        rf = read_rf(dat)
        if len(rf) != len(disp):
            sys.exit(f"ABORT: parsed {len(rf)} reaction nodes, expected {len(disp)}")
        work = sum(rf[n][0] * disp[n][0] + rf[n][1] * disp[n][1] for n in disp)
        e1, e2, e3 = e
        A_mat.append([e1 * e1, e2 * e2, e3 * e3, 2 * e1 * e2, 2 * e1 * e3, 2 * e2 * e3])
        rhs.append(work / (S * S * T * SIDE * SIDE))
    sol = np.linalg.solve(np.array(A_mat, float), np.array(rhs))
    fe = dict(zip(["Q11", "Q22", "Q66", "Q12", "Q16", "Q26"], sol))
    an = analytic()
    print(f"\nrotation of material axis from global x: {np.degrees(THETA):.1f} deg")
    print("entry   analytic [GPa]    CalculiX [GPa]   diff / Q11")
    ok = True
    for key in ["Q11", "Q22", "Q66", "Q12", "Q16", "Q26"]:
        d = (fe[key] - an[key]) / an["Q11"]
        ok &= abs(d) <= 1e-3
        print(f"{key}   {an[key]/1e9:14.4f}   {fe[key]/1e9:14.4f}   {d:+.2e}")
    print("OVERALL:", "PASS" if ok else "FAIL")


if __name__ == "__main__":
    main()
