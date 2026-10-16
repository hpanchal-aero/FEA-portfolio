#!/usr/bin/env python3
"""Does CalculiX apply *ORIENTATION to S4 shells like Nastran's MCID?

Nastran rule: material 1-axis = projection of the MCID system's x-axis onto the element plane.
Test: orthotropic strip in a plane spanned by u=(cos phi,0,sin phi) and v=(0,1,0), orientation
a=(0,0,1) (x-axis), b=(0,1,0) (point in the xy-plane). Projection of z onto the plane is
parallel to u for any phi != 0, so material-1 is expected along u and material-2 along v.
"""
import glob
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

PROJ = Path(__file__).resolve().parent.parent
WORK = PROJ / "simulation" / "orientation_test"

E1, E2, E3 = 70e9, 30e9, 30e9
NU12 = NU13 = NU23 = 0.3
G12, G13, G23 = 11.67e9, 16.08e9, 9.72e9
T, L, W, NL, NW, F = 0.01, 2.0, 0.5, 16, 4, 500.0


def find_ccx():
    env = os.environ.get("CCX")
    if env:
        return env
    cands = [c for c in sorted(glob.glob(str(Path.home() / "miniconda3/envs/ccx223/bin/ccx*")))
             if Path(c).is_file() and os.access(c, os.X_OK) and not c.endswith(".so")]
    if not cands:
        sys.exit("ABORT: no ccx binary found in the ccx223 env; set CCX=/path/to/ccx")
    return cands[0]


def write_case(name, phi_deg, long_axis):
    phi = np.radians(phi_deg)
    u = np.array([np.cos(phi), 0.0, np.sin(phi)])
    v = np.array([0.0, 1.0, 0.0])
    ell, m = (u, v) if long_axis == "u" else (v, u)

    ids, nid = {}, 1
    coords = []
    for j in range(NW + 1):
        for i in range(NL + 1):
            ids[(i, j)] = nid
            coords.append((nid, ell * (L * i / NL) + m * (W * j / NW)))
            nid += 1
    elems = []
    for j in range(NW):
        for i in range(NL):
            elems.append([ids[(i, j)], ids[(i + 1, j)], ids[(i + 1, j + 1)], ids[(i, j + 1)]])
    clamp = [ids[(0, j)] for j in range(NW + 1)]
    load = [ids[(NL, j)] for j in range(NW + 1)]

    def nset(label, lst):
        out = [f"*NSET, NSET={label}"]
        for k in range(0, len(lst), 8):
            out.append(", ".join(str(x) for x in lst[k:k + 8]))
        return out

    lines = ["*NODE"]
    lines += [f"{n}, {p[0]:.10f}, {p[1]:.10f}, {p[2]:.10f}" for n, p in coords]
    lines += ["*ELEMENT, TYPE=S4, ELSET=EALL"]
    lines += [f"{k + 1}, " + ", ".join(str(n) for n in e) for k, e in enumerate(elems)]
    lines += nset("CLAMP", clamp) + nset("LOAD", load)
    lines += ["*ORIENTATION, NAME=OR1, SYSTEM=RECTANGULAR", "0.,0.,1.,0.,1.,0.",
              "*MATERIAL, NAME=MAT", "*ELASTIC, TYPE=ENGINEERING CONSTANTS",
              f"{E1:.6e}, {E2:.6e}, {E3:.6e}, {NU12}, {NU13}, {NU23}, {G12:.6e}, {G13:.6e}",
              f"{G23:.6e}, 0.",
              "*SHELL SECTION, ELSET=EALL, MATERIAL=MAT, ORIENTATION=OR1", f"{T}",
              "*STEP", "*STATIC", "*BOUNDARY", "CLAMP, 1, 6", "*CLOAD"]
    for n in load:
        for k in range(3):
            val = F / len(load) * ell[k]
            if abs(val) > 1e-12:
                lines.append(f"{n}, {k + 1}, {val:.8e}")
    lines += ["*NODE PRINT, NSET=LOAD", "U", "*END STEP"]
    (WORK / f"{name}.inp").write_text("\n".join(lines) + "\n")
    return ell, load


def read_disp(dat, nodes):
    disp, grab = {}, False
    for ln in dat.read_text().splitlines():
        if "displacements" in ln:
            grab = True
            continue
        if grab:
            tok = ln.split()
            if len(tok) == 4:
                try:
                    disp[int(tok[0])] = np.array([float(t) for t in tok[1:]])
                    continue
                except ValueError:
                    pass
            if disp and not ln.strip():
                break
    return np.mean([disp[n] for n in nodes], axis=0)


def main():
    ccx = find_ccx()
    WORK.mkdir(parents=True, exist_ok=True)
    print("ccx binary:", ccx)
    print(subprocess.run([ccx, "-v"], capture_output=True, text=True).stdout.strip()[:200])
    print(f"strip L={L} m, W={W} m, t={T} m, F={F} N, E1={E1/1e9:.0f} GPa, E2={E2/1e9:.0f} GPa\n")

    res = {}
    for phi in (90.0, 35.0):
        for ax in ("u", "v"):
            name = f"phi{int(phi)}_{ax}"
            ell, load = write_case(name, phi, ax)
            r = subprocess.run([ccx, name], cwd=WORK, capture_output=True, text=True, timeout=600)
            dat = WORK / f"{name}.dat"
            if r.returncode != 0 or not dat.exists():
                print(r.stdout[-1500:], r.stderr[-1500:])
                sys.exit(f"ABORT: ccx failed for {name} (rc={r.returncode})")
            d = read_disp(dat, load)
            delta = float(np.dot(d, ell))
            pred = F * L / ((E1 if ax == "u" else E2) * T * W)
            res[(phi, ax)] = (delta, pred)
            print(f"{name:10s} load along {ax}: delta = {delta:.4e} m   predicted {pred:.4e} m   "
                  f"ratio {delta / pred:.3f}   |off-axis| = {np.linalg.norm(d - delta * ell):.2e} m")

    print()
    ok_all = True
    for phi in (90.0, 35.0):
        du, dv = res[(phi, "u")][0], res[(phi, "v")][0]
        r = dv / du
        ok = (2.1 <= r <= 2.5) and all(abs(res[(phi, a)][0] / res[(phi, a)][1] - 1) <= 0.08 for a in ("u", "v"))
        ok_all &= ok
        print(f"phi={phi:4.0f}: delta_v/delta_u = {r:.3f}  -> {'PASS' if ok else 'FAIL'}")
    print("OVERALL:", "PASS" if ok_all else "FAIL")


if __name__ == "__main__":
    main()
