#!/usr/bin/env python3
"""Independent force balance of a CalculiX static deck (gravity by *DLOAD GRAV, point masses as T3D2 truss pairs tied by
*EQUATION). Read-only: parses <job>.inp and <job>.dat.

Usage: check_equilibrium.py /path/to/jobname      (without extension)

Accounting
  external load: CLOAD sum A, gravity W = sum over elements of mass * g (shells: area * t * rho, area = |d1 x d2| / 2 from the
  diagonals, 1/4 of the weight per node; T3D2: rho * A * L, 1/2 per node). Loads on dependent DOFs of *EQUATION blocks (first
  term dependent) are transferred to the independent DOFs: f_k += -(a_k / a_1) f_dep, recursively.
  printed RF at the print set NSPC is K*u: reaction + external load at those nodes. So
      reaction at NSPC = printed total - (CLOAD + transferred gravity) on all DOFs of NSPC nodes
  Closure: reaction + (A + W) = 0.
"""
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

base = Path(sys.argv[1]).resolve()
inp, dat = base.with_suffix(".inp"), base.with_suffix(".dat")
SPREADSHEET_SHELL_MASS = 11958.5

coords, elems, esets = {}, {}, defaultdict(list)
nspc, sect, dens, dload, cloads, eqs = set(), {}, {}, [], [], []
kw, opts, cur_mat, cur_sect, eq_cur = None, {}, None, None, None
for ln in inp.read_text().splitlines():
    s = ln.strip()
    if not s or s.startswith("**"):
        continue
    if s.startswith("*"):
        parts = [p.strip() for p in s.split(",")]
        kw = parts[0].upper()
        opts = {}
        for p in parts[1:]:
            if "=" in p:
                a, b = p.split("=", 1)
                opts[a.strip().upper()] = b.strip()
        if kw == "*MATERIAL":
            cur_mat = opts.get("NAME")
        if kw in ("*SHELL SECTION", "*SOLID SECTION"):
            cur_sect = opts["ELSET"]
            sect[cur_sect] = dict(kind=kw, mat=opts.get("MATERIAL"), value=None)
        eq_cur = None
        continue
    t = [x.strip() for x in s.split(",") if x.strip()]
    if kw == "*NODE":
        coords[int(t[0])] = np.array([float(t[1]), float(t[2]), float(t[3])])
    elif kw == "*ELEMENT":
        eid = int(t[0])
        elems[eid] = (opts.get("TYPE", "").upper(), [int(x) for x in t[1:]])
        esets[opts.get("ELSET")].append(eid)
    elif kw == "*NSET" and opts.get("NSET", "").upper() == "NSPC":
        nspc.update(int(x) for x in t)
    elif kw in ("*SHELL SECTION", "*SOLID SECTION"):
        sect[cur_sect]["value"] = float(t[0])
    elif kw == "*DENSITY":
        dens[cur_mat] = float(t[0])
    elif kw == "*DLOAD":
        if t[1].upper() == "GRAV":
            dload.append((t[0], float(t[2]), np.array([float(x) for x in t[3:6]])))
    elif kw == "*CLOAD":
        cloads.append((int(t[0]), int(t[1]), float(t[2])))
    elif kw == "*EQUATION":
        if eq_cur is None:
            eq_cur = dict(n=int(t[0]), terms=[])
        else:
            for i in range(0, len(t), 3):
                eq_cur["terms"].append((int(t[i]), int(t[i + 1]), float(t[i + 2])))
            if len(eq_cur["terms"]) == eq_cur["n"]:
                eqs.append(eq_cur["terms"])
                eq_cur = None

dep = {}
for terms in eqs:
    n0, d0, a0 = terms[0]
    dep[(n0, d0)] = [(n, d, -a / a0) for n, d, a in terms[1:]]


def transfer(key, f, out, depth=0):
    if depth > 20:
        sys.exit(f"ABORT: equation chain too deep at {key}")
    if key in dep:
        for n, d, r in dep[key]:
            transfer((n, d), r * f, out, depth + 1)
    else:
        out[key] += f


def elem_mass(elset):
    info = sect[elset]
    rho = dens[info["mat"]]
    tot, per = 0.0, []
    for eid in esets[elset]:
        etype, ns = elems[eid]
        P = [coords[n] for n in ns]
        if info["kind"] == "*SHELL SECTION":
            area = 0.5 * np.linalg.norm(np.cross(P[2] - P[0], P[3] - P[1]))
            m = area * info["value"] * rho
        else:
            m = rho * info["value"] * np.linalg.norm(P[1] - P[0])
        per.append((ns, m))
        tot += m
    return tot, per


m_shell = m_truss = 0.0
grav = defaultdict(float)
for elset, mag, vec in dload:
    tot, per = elem_mass(elset)
    if sect[elset]["kind"] == "*SHELL SECTION":
        m_shell += tot
    else:
        m_truss += tot
    g = mag * vec
    for ns, m in per:
        for n in ns:
            for d in range(3):
                if g[d] != 0.0:
                    transfer((n, d + 1), m * g[d] / len(ns), grav)
W = np.zeros(3)
W_nspc = np.zeros(3)
for (n, d), f in grav.items():
    W[d - 1] += f
    if n in nspc:
        W_nspc[d - 1] += f

cl = defaultdict(float)
for n, d, v in cloads:
    transfer((n, d), v, cl)
A, A_nspc = np.zeros(3), np.zeros(3)
for (n, d), f in cl.items():
    A[d - 1] += f
    if n in nspc:
        A_nspc[d - 1] += f

rf_tot, mode = None, None
for ln in dat.read_text().splitlines():
    s = ln.strip()
    if s.startswith("total force"):
        mode = "RFT"
        continue
    tok = ln.split()
    if mode == "RFT" and len(tok) == 3:
        rf_tot = np.array([float(x) for x in tok])
        mode = None
if rf_tot is None:
    sys.exit("ABORT: no 'total force' block in the .dat file")

R = rf_tot - A_nspc - W_nspc
res = R + A + W
print(f"deck: {len(coords)} nodes, {len(elems)} elements, {len(eqs)} equations, {len(nspc)} NSPC nodes")
print(f"shell mass in deck     : {m_shell:12.2f} kg   (spreadsheet {SPREADSHEET_SHELL_MASS:.1f} kg, ratio {m_shell / SPREADSHEET_SHELL_MASS:.4f})")
print(f"truss (point) masses   : {m_truss:12.2f} kg")
print(f"CLOAD sum A            : Fx={A[0]:+.4e} Fy={A[1]:+.4e} Fz={A[2]:+.4e} N")
print(f"gravity sum W (after equation transfer): Fx={W[0]:+.4e} Fy={W[1]:+.4e} Fz={W[2]:+.4e} N")
print(f"  of which on NSPC nodes (all DOFs): Fz={W_nspc[2]:+.4e} N   CLOAD on NSPC nodes: Fz={A_nspc[2]:+.4e} N")
print(f"printed RF total       : Fx={rf_tot[0]:+.4e} Fy={rf_tot[1]:+.4e} Fz={rf_tot[2]:+.4e} N")
print(f"reaction at NSPC       : Fx={R[0]:+.4e} Fy={R[1]:+.4e} Fz={R[2]:+.4e} N")
print(f"closure R + A + W      : Fx={res[0]:+.3e} Fy={res[1]:+.3e} Fz={res[2]:+.3e} N")
net = abs(A[2] + W[2])
print(f"net vertical load |A + W| = {net:.4e} N;  |Fz closure| = {abs(res[2]):.3e} N = {100 * abs(res[2]) / net:.3f} % of it")
