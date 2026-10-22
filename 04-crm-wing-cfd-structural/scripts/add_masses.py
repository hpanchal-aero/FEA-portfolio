#!/usr/bin/env python3
"""Add the attached masses of the uCRM-9 (engine, LE/TE, fuel) to the CalculiX deck as real point masses.

NASTRAN: CONM2 point masses (zero rotational inertia) tied to the wing by RBE3 (translations, equal weights,
reference DOFs 123456); the engine masses are in addition tied rigidly (RBE2) to the master grid 20000056.
CalculiX: one 1-node MASS element per CONM2. The translations of every mass node are written with *EQUATION
in terms of the wing-node translations using the exact RBE3 least-squares map
    x = (u_ref, theta) = N^-1 sum_i w_i G_i^T u_i,   G_i = [I, -skew(d_i)],   N = sum_i w_i G_i^T G_i
    RBE3 mass node (reference grid):   u_ref = sum_i Cu_i u_i
    engine grid e (RBE2 to master m):  u_e = u_m + theta x d_e = sum_i (Cu_i - skew(d_e) Ct_i) u_i
Modes:
    static : base deck = static_v1 (aero nodal forces); adds gravity (*DLOAD GRAV) on shells and masses
             -> simulation/static_v3
    modal  : *FREQUENCY step with the same constraints and masses, no loads  -> simulation/modal_v1
Usage: /usr/bin/python3 scripts/add_masses.py coarse static|modal [--modes 10] [--fuel 20|full]
"""
import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from ucrm_bdf import ival, num, read_cards  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("level", choices=["coarse", "medium", "fine"])
ap.add_argument("mode", choices=["static", "modal"])
ap.add_argument("--modes", type=int, default=10)
ap.add_argument("--fuel", choices=["20", "full"], default="20")
args = ap.parse_args()
level, mode = args.level, args.mode

root = HERE.parent
fem = root / "data" / "ucrm9" / "uCRM9" / "Mesh" / level.capitalize() / "FEM"
base = root / "simulation" / "static_v1" / level / f"uCRM9_{level}_v1.inp"
if not base.exists():
    sys.exit(f"ABORT: {base} not found (run ucrm_to_ccx.py {level} first)")
if mode == "static":
    out = root / "simulation" / "static_v3" / level / f"uCRM9_{level}_v3.inp"
else:
    out = root / "simulation" / "modal_v1" / level / f"uCRM9_{level}_modal.inp"
out.parent.mkdir(parents=True, exist_ok=True)

# ---- base deck: split at the (single) *STEP ---------------------------------------------------------
lines = base.read_text().splitlines()
steps = [k for k, ln in enumerate(lines) if ln.strip().upper() == "*STEP"]
if len(steps) != 1:
    sys.exit("ABORT: expected exactly one *STEP in the base deck")
head, step = lines[:steps[0]], lines[steps[0]:]

wing, elsets, kw, max_el = {}, [], None, 0
for ln in head:
    s = ln.strip()
    if s.startswith("*"):
        kw = s.split(",")[0].upper()
        if kw == "*ELEMENT":
            m = re.search(r"ELSET=([A-Za-z0-9_]+)", s, re.I)
            if m:
                elsets.append(m.group(1))
        continue
    if kw == "*NODE" and s:
        t = [x.strip() for x in s.split(",")]
        wing[int(t[0])] = np.array([float(t[1]), float(t[2]), float(t[3])])
    elif kw == "*ELEMENT" and s:
        max_el = max(max_el, int(s.split(",")[0]))
if not wing or not elsets:
    sys.exit("ABORT: could not read nodes or element sets from the base deck")

bnd, kw = [], None
for ln in step:
    s = ln.strip()
    if s.startswith("*"):
        kw = s.split(",")[0].upper()
        continue
    if kw == "*BOUNDARY" and s:
        bnd.append(s)

# ---- mass files ----------------------------------------------------------------------------------------
fuel_file = f"uCRM-9_fuel_mass_{'20_percent' if args.fuel == '20' else 'full'}_{level}.bdf"
files = (("LE/TE", f"uCRM-9_LETE_masses_{level}.bdf"),
         ("engine", f"uCRM-9_engine_mass_{level}.bdf"),
         ("fuel", fuel_file))
mass_grid, conm, grp, rbe2m, rbe3 = {}, {}, {}, {}, []
for label, fn in files:
    for name, f in read_cards(fem / fn):
        if name == "GRID":
            if ival(f[1]) != 0 or ival(f[5]) != 0:
                sys.exit(f"ABORT: mass GRID with non-zero CP/CD in {fn}")
            g = ival(f[0])
            if g in mass_grid:
                sys.exit(f"ABORT: duplicate mass grid {g}")
            mass_grid[g] = np.array([num(f[2]), num(f[3]), num(f[4])])
            grp[g] = label
        elif name == "CONM2":
            if ival(f[2]) != 0 or any(num(x) != 0 for x in f[4:7]) or any(num(x) != 0 for x in f[8:14]):
                sys.exit(f"ABORT: CONM2 with CID, offset or inertia in {fn}")
            g = ival(f[1])
            conm[g] = conm.get(g, 0.0) + num(f[3])
        elif name == "RBE2":
            if f[2].strip() != "123456":
                sys.exit(f"ABORT: RBE2 with CM != 123456 in {fn}")
            for x in f[3:]:
                if x.strip():
                    rbe2m[ival(x)] = ival(f[1])
        elif name == "RBE3":
            if f[3].strip() != "123456" or f[5].strip() != "123" or any("." in x for x in f[6:]):
                sys.exit(f"ABORT: unsupported RBE3 layout in {fn}")
            rbe3.append((ival(f[2]), num(f[4]), [ival(x) for x in f[6:] if x.strip()]))

if set(mass_grid) != set(conm):
    sys.exit("ABORT: mass grids and CONM2 cards disagree")
if set(mass_grid) & set(wing):
    sys.exit("ABORT: a mass grid id clashes with a wing node id")
refs = {r for r, _, _ in rbe3}
if refs != {g for g in conm if g not in rbe2m}:
    sys.exit("ABORT: RBE3 reference grids and the masses that are not RBE2 slaves disagree")
if any(m not in refs for m in rbe2m.values()):
    sys.exit("ABORT: an RBE2 master is not an RBE3 reference grid")


def skew(d):
    return np.array([[0.0, -d[2], d[1]], [d[2], 0.0, -d[0]], [-d[1], d[0], 0.0]])


# ---- exact RBE3 maps -----------------------------------------------------------------------------------
rb = {}
rng = np.random.default_rng(1)
for ref, w, deps in rbe3:
    if any(n not in wing for n in deps):
        sys.exit(f"ABORT: RBE3 {ref}: dependent node missing from the wingbox")
    d = np.array([wing[n] for n in deps]) - mass_grid[ref]
    Gs, N = [], np.zeros((6, 6))
    for di in d:
        G = np.hstack([np.eye(3), -skew(di)])
        Gs.append(G)
        N += w * G.T @ G
    if np.linalg.cond(N) > 1e10:
        sys.exit(f"ABORT: RBE3 {ref}: ill-conditioned (cond {np.linalg.cond(N):.2e})")
    Ninv = np.linalg.inv(N)
    Cu, Ct = [], []
    for G in Gs:
        Ci = Ninv @ (w * G.T)
        Cu.append(Ci[:3])
        Ct.append(Ci[3:])
    a, b = rng.normal(size=3), rng.normal(size=3)
    ui = [a + np.cross(b, di) for di in d]
    ur = sum(C @ u for C, u in zip(Cu, ui))
    th = sum(C @ u for C, u in zip(Ct, ui))
    if np.linalg.norm(ur - a) > 1e-9 * max(1.0, np.linalg.norm(a)) or np.linalg.norm(th - b) > 1e-9 * max(1.0, np.linalg.norm(b)):
        sys.exit(f"ABORT: RBE3 {ref}: the map does not reproduce rigid-body motion")
    rb[ref] = (deps, Cu, Ct)

# ---- equations -----------------------------------------------------------------------------------------
eqs, neq = [], 0
new_id = {g: max(wing) + 1 + k for k, g in enumerate(sorted(mass_grid))}   # contiguous labels after the wing nodes


def add_eq(node, k, deps, C):
    global neq
    cmax = max(1.0, max(abs(C[i][k, l]) for i in range(len(deps)) for l in range(3)))
    terms = [(new_id[node], k + 1, 1.0)]
    for i, dep in enumerate(deps):
        for l in range(3):
            c = C[i][k, l]
            if abs(c) > 1e-12 * cmax:
                terms.append((dep, l + 1, -float(c)))
    eqs.append(str(len(terms)))
    for q in range(0, len(terms), 4):
        eqs.append(", ".join(f"{n}, {d_}, {c:.12e}" for n, d_, c in terms[q:q + 4]))
    neq += 1


for g in sorted(conm):
    if g in rbe2m:
        m = rbe2m[g]
        deps, Cu, Ct = rb[m]
        de = mass_grid[g] - mass_grid[m]
        C = [Cu[i] - skew(de) @ Ct[i] for i in range(len(deps))]
    else:
        deps, Cu, Ct = rb[g]
        C = Cu
    for k in range(3):
        add_eq(g, k, deps, C)

# ---- deck ----------------------------------------------------------------------------------------------
mnodes = sorted(mass_grid)
txt = list(head)
txt += ["** attached masses: MASS elements; RBE3 / RBE2 expressed with *EQUATION (exact least-squares map)",
        "*NODE, NSET=NMASS"]
txt += [f"{new_id[g]}, {mass_grid[g][0]:.9e}, {mass_grid[g][1]:.9e}, {mass_grid[g][2]:.9e}" for g in mnodes]
mass_sets = []
for k, g in enumerate(mnodes):
    nm = f"MS{k}"
    mass_sets.append(nm)
    txt += [f"*ELEMENT, TYPE=MASS, ELSET={nm}", f"{max_el + 1 + k}, {new_id[g]}", f"*MASS, ELSET={nm}", f"{conm[g]:.9e}"]
txt += ["*EQUATION"] + eqs

if mode == "static":
    inode = [k for k, ln in enumerate(step) if ln.strip().upper() == "*NODE FILE"]
    if len(inode) != 1:
        sys.exit("ABORT: expected exactly one *NODE FILE in the base step")
    dl = ["*DLOAD"] + [f"{nm}, GRAV, 9.81, 0., 0., -1." for nm in elsets + mass_sets]
    txt += step[:inode[0]] + dl + step[inode[0]:]
else:
    txt += ["*STEP", "*FREQUENCY", str(args.modes), "*BOUNDARY"] + bnd + ["*NODE FILE", "U", "*END STEP"]

out.write_text("\n".join(txt) + "\n")

print(f"wrote {out}")
print(f"mass nodes renumbered {min(new_id.values())}..{max(new_id.values())}; mass elements {max_el + 1}..{max_el + len(conm)} "
      f"(CalculiX sizes its arrays by the highest label)")
print(f"wing nodes {len(wing)}  mass nodes {len(mass_grid)}  MASS elements {len(conm)}  "
      f"RBE3 {len(rbe3)}  RBE2 slaves {len(rbe2m)}  equations {neq}")
print(f"RBE3 least-squares map reproduces rigid-body motion for all {len(rb)} RBE3 cards (asserted)")
tot = 0.0
for lab in ("LE/TE", "engine", "fuel"):
    m = sum(conm[g] for g in conm if grp[g] == lab)
    tot += m
    print(f"{lab:7s}: mass {m:10.1f} kg   weight {m * 9.81:12.1f} N")
print(f"total attached mass {tot:.1f} kg   weight {tot * 9.81:.1f} N  (shell weight is applied by *DLOAD in the static deck)")
