#!/usr/bin/env python3
"""Read-only consistency checks of the uCRM-9 FEM input files.

Usage: check_ucrm_inputs.py [coarse] [medium] [fine]
"""
import math
import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ucrm_bdf import ival, num, read_cards  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent / "data" / "ucrm9" / "uCRM9" / "Mesh"
LEVELS = sys.argv[1:] or ["coarse"]

# Values quoted from uCRM-9_Design_Specifications.xlsx / Results.xlsx
REF_DOF = {"coarse": 142428, "medium": 356328, "fine": 571716}
WING_STRUCT_MASS_BOTH = 34482.0   # kg, both wings, includes x1.5 factor per footnote
CL, MACH, S_REF, ALT_FT = 0.5, 0.85, 383.74, 37000.0


def isa_pressure(h_m):
    """ISA static pressure, 11-20 km layer (isothermal)."""
    return 22632.06 * math.exp(-(h_m - 11000.0) / 6341.62)


for level in LEVELS:
    fem = ROOT / level.capitalize() / "FEM"
    print("=" * 70)
    print(f"LEVEL: {level}   ({fem})")
    cards = read_cards(fem / f"uCRM-9_wingbox_{level}.bdf")

    grids, quads, pshell, mat2, cord, spc = {}, [], {}, {}, {}, []
    nonzero_cs = 0
    for name, f in cards:
        if name == "GRID":
            grids[ival(f[0])] = np.array([num(f[2]), num(f[3]), num(f[4])])
            if ival(f[1]) != 0 or ival(f[5]) != 0:
                nonzero_cs += 1
        elif name == "CQUADR":
            quads.append((ival(f[0]), ival(f[1]), [ival(f[i]) for i in range(2, 6)], ival(f[6])))
        elif name == "PSHELL":
            pshell[ival(f[0])] = dict(mid1=ival(f[1]), t=num(f[2]), mid2=ival(f[3]),
                                      r12=num(f[4]), mid3=ival(f[5]), tst=num(f[6]),
                                      nsm=num(f[7]), mid4=ival(f[10]))
        elif name == "MAT2":
            mat2[ival(f[0])] = dict(g11=num(f[1]), g12=num(f[2]), g13=num(f[3]),
                                    g22=num(f[4]), g23=num(f[5]), g33=num(f[6]),
                                    rho=num(f[7]))
        elif name == "CORD2R":
            cord[ival(f[0])] = [num(f[i]) for i in range(2, 11)]
        elif name == "SPC":
            i = 1
            while i + 1 < len(f):
                g = ival(f[i])
                if g > 0:
                    spc.append((g, f[i + 1].strip()))
                i += 3

    print(f"nodes {len(grids)}   quads {len(quads)}   PSHELL {len(pshell)}   "
          f"MAT2 {len(mat2)}   CORD2R {len(cord)}   SPC entries {len(spc)}")
    print(f"nodes with non-zero CP/CD: {nonzero_cs}")
    print(f"DOF = 6 x nodes = {6 * len(grids)}   reference = {REF_DOF[level]}   "
          f"match: {6 * len(grids) == REF_DOF[level]}")

    used = {n for _, _, ns, _ in quads for n in ns}
    print(f"nodes used by quads: {len(used)}   unused nodes: {len(grids) - len(used)}   "
          f"quad nodes missing from GRID: {len(used - set(grids))}")
    print(f"elements whose MCID equals PID: {sum(1 for _, p, _, m in quads if m == p)} / {len(quads)}")
    missing_p = {p for _, p, _, _ in quads} - set(pshell)
    missing_m = {m for p in pshell.values() for m in (p['mid1'], p['mid2'], p['mid3'], p['mid4'])} - set(mat2)
    print(f"PIDs used but undefined: {sorted(missing_p)}   MIDs referenced but undefined: {sorted(missing_m)}")

    xyz = np.array(list(grids.values()))
    print("bounding box [m]  x: %.3f..%.3f   y: %.3f..%.3f   z: %.3f..%.3f" % (
        xyz[:, 0].min(), xyz[:, 0].max(), xyz[:, 1].min(), xyz[:, 1].max(),
        xyz[:, 2].min(), xyz[:, 2].max()))

    # --- membrane matrix -> engineering constants -------------------------------------
    rows, bad, g13max = [], 0, 0.0
    for pid, p in sorted(pshell.items()):
        m = mat2[p["mid1"]]
        g13max = max(g13max, abs(m["g13"]), abs(m["g23"]))
        G = np.array([[m["g11"], m["g12"]], [m["g12"], m["g22"]]])
        if m["g11"] <= 0 or np.linalg.det(G) <= 0:
            bad += 1
            continue
        Sm = np.linalg.inv(G)
        rows.append((p["t"], 1 / Sm[0, 0], 1 / Sm[1, 1], -Sm[0, 1] / Sm[0, 0], m["g33"], m["rho"]))
    a = np.array(rows)
    print(f"\nmembrane matrices: {len(rows)} positive-definite, {bad} not;  "
          f"max |G13|,|G23| = {g13max:.3e} Pa")
    print("property      t [mm]    E1 [GPa]   E2 [GPa]   nu12     G12 [GPa]  rho [kg/m3]")
    for lab, fn in (("min", np.min), ("median", np.median), ("max", np.max)):
        r = fn(a, axis=0)
        print(f"{lab:8s} {r[0]*1e3:10.3f} {r[1]/1e9:11.3f} {r[2]/1e9:10.3f} {r[3]:8.4f} {r[4]/1e9:10.3f} {r[5]:12.1f}")
    for pid in (1, 2):
        p = pshell[pid]
        print(f"PID {pid}: t={p['t']*1e3:.3f} mm  12I/t3={p['r12']:.4f}  TS/T={p['tst']:.3f}  "
              f"NSM={p['nsm']:.3g}  MIDs=({p['mid1']},{p['mid2']},{p['mid3']},{p['mid4']})")

    # --- shell mass ----------------------------------------------------------------------
    mass = 0.0
    for _, pid, ns, _ in quads:
        p1, p2, p3, p4 = (grids[n] for n in ns)
        area = 0.5 * np.linalg.norm(np.cross(p3 - p1, p4 - p2))
        pr = pshell[pid]
        mass += area * (pr["t"] * mat2[pr["mid1"]]["rho"] + pr["nsm"])
    print(f"\nshell mass (area x t x rho), one wing : {mass:10.1f} kg")
    print(f"spreadsheet wing mass / 2             : {WING_STRUCT_MASS_BOTH / 2:10.1f} kg  (includes x1.5 factor)")
    print(f"spreadsheet wing mass / 2 / 1.5       : {WING_STRUCT_MASS_BOTH / 2 / 1.5:10.1f} kg")
    print(f"model mass x 1.5                      : {mass * 1.5:10.1f} kg")

    # --- boundary conditions -----------------------------------------------------------------
    sn = sorted({g for g, _ in spc})
    print(f"\nSPC: {len(sn)} distinct nodes, DOF strings {dict(Counter(c for _, c in spc))}")
    if sn:
        pts = np.array([grids[g] for g in sn if g in grids])
        print("SPC node bounding box [m]  x: %.3f..%.3f   y: %.3f..%.3f   z: %.3f..%.3f" % (
            pts[:, 0].min(), pts[:, 0].max(), pts[:, 1].min(), pts[:, 1].max(),
            pts[:, 2].min(), pts[:, 2].max()))
        print(f"SPC nodes missing from GRID: {sum(1 for g in sn if g not in grids)}")
        by = {}
        for g, c in spc:
            by.setdefault(c, []).append(grids[g][1])
        for c, ys in sorted(by.items()):
            print(f"   DOF {c:>7s}: {len(ys):5d} entries, y = {min(ys):.3f}..{max(ys):.3f} m")

    # --- loads and attached masses --------------------------------------------------------------
    print()
    for fname in (f"uCRM-9_aero_forces_1.0g_{level}.bdf",):
        fc = read_cards(fem / fname)
        F = np.zeros(3)
        nf, sids, offgrid = 0, set(), 0
        for name, f in fc:
            if name == "FORCE":
                mag = num(f[3])
                F += mag * np.array([num(f[4]), num(f[5]), num(f[6])])
                nf += 1
                sids.add(ival(f[0]))
                offgrid += ival(f[1]) not in grids
        print(f"{fname}: {nf} FORCE cards, SIDs {sorted(sids)}, nodes not in GRID: {offgrid}")
        print(f"   sum Fx = {F[0]:.4e} N   Fy = {F[1]:.4e} N   Fz = {F[2]:.4e} N")
    p_inf = isa_pressure(ALT_FT * 0.3048)
    q = 0.5 * 1.4 * p_inf * MACH ** 2
    print(f"   cross-check: ISA p({ALT_FT:.0f} ft) = {p_inf:.1f} Pa, q = {q:.1f} Pa, "
          f"CL*q*S/2 = {CL * q * S_REF / 2:.4e} N (approximate, not from the dataset)")

    for fname in (f"uCRM-9_engine_mass_{level}.bdf", f"uCRM-9_LETE_masses_{level}.bdf",
                  f"uCRM-9_fuel_mass_20_percent_{level}.bdf", f"uCRM-9_fuel_mass_full_{level}.bdf"):
        fc = read_cards(fem / fname)
        cnt = Counter(n for n, _ in fc)
        m_tot = sum(num(f[3]) for n, f in fc if n == "CONM2")
        print(f"{fname}: cards {dict(cnt)}   total CONM2 mass = {m_tot:.1f} kg")
