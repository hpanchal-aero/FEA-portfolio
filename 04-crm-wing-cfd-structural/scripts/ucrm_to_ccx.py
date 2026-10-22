#!/usr/bin/env python3
"""Convert the uCRM-9 wingbox (Nastran BDF, smeared-stiffener CQUADR shells) to a CalculiX S4 deck.

  default      : v1 = shells + orientations + SPC + aerodynamic nodal forces        -> simulation/static_v1
  --all-loads  : v2 = v1 + weight of the shells + engine, LE/TE and fuel masses      -> simulation/static_v2
                 as statically equivalent nodal forces (RBE2/RBE3 load transfer in closed form).

SI units (m, N, kg, Pa). Modelling choices (documented assumptions):
  * Membrane stiffness: MAT2 of MID1 inverted to E1, E2, nu12; G12 = MID1 G33.
  * Transverse shear: G13 = MID3 G11, G23 = MID3 G22.
  * E3, nu13, nu23 are placeholders required for CalculiX's 3D constants (E3 = E2, nu = 0.3).
  * Bending matrix (MID2), coupling matrix (MID4) and 12I/t^3 are NOT represented.
  * RBE3 (translations only, equal weights) is a load-sharing relation, not a stiffness: a load (F, M)
    at the reference node becomes f_i = w_i (a + b x d_i) on the dependent nodes, with sum f_i = F and
    sum d_i x f_i = M. The script asserts both.
  * Shell weight is applied as lumped nodal forces (rho t A / 4 per node).
"""
import argparse
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from ucrm_bdf import ival, num, read_cards  # noqa: E402

G_ACC = np.array([0.0, 0.0, -9.81])   # m/s^2, as GRAV 2 0 -9.81 0. 0. 1.0 in the NASTRAN deck

ap = argparse.ArgumentParser()
ap.add_argument("level", choices=["coarse", "medium", "fine"])
ap.add_argument("--all-loads", action="store_true")
ap.add_argument("--fuel", choices=["20", "full"], default="20")
ap.add_argument("--e3-factor", type=float, default=1.0, help="sensitivity: E3 = factor * E2 (default 1)")
ap.add_argument("--nu3", type=float, default=0.3, help="sensitivity: nu13 = nu23 (default 0.3)")
ap.add_argument("--gshear-scale", type=float, default=1.0, help="sensitivity: scale of G13 and G23 (default 1)")
ap.add_argument("--tag", default="", help="suffix for the output directory and job name")
args = ap.parse_args()
level, all_loads = args.level, args.all_loads
tag = "v2" if all_loads else "v1"
if args.tag:
    tag = f"{tag}_{args.tag}"

fem = HERE.parent / "data" / "ucrm9" / "uCRM9" / "Mesh" / level.capitalize() / "FEM"
out_dir = HERE.parent / "simulation" / f"static_{tag}" / level
out_dir.mkdir(parents=True, exist_ok=True)
out_path = out_dir / f"uCRM9_{level}_{tag}.inp"

cards = read_cards(fem / f"uCRM-9_wingbox_{level}.bdf")
grids, quads, pshell, mat2, cord, spc = {}, [], {}, {}, {}, []
for name, f in cards:
    if name == "GRID":
        if ival(f[1]) != 0 or ival(f[5]) != 0:
            sys.exit("ABORT: GRID with non-zero CP/CD found")
        grids[ival(f[0])] = (num(f[2]), num(f[3]), num(f[4]))
    elif name == "CQUADR":
        quads.append((ival(f[0]), ival(f[1]), [ival(f[i]) for i in range(2, 6)], ival(f[6])))
    elif name == "PSHELL":
        pshell[ival(f[0])] = dict(mid1=ival(f[1]), t=num(f[2]), mid2=ival(f[3]),
                                  mid3=ival(f[5]), nsm=num(f[7]), mid4=ival(f[10]))
    elif name == "MAT2":
        mat2[ival(f[0])] = dict(g11=num(f[1]), g12=num(f[2]), g22=num(f[4]), g33=num(f[6]), rho=num(f[7]))
    elif name == "CORD2R":
        cord[ival(f[0])] = [np.array([num(f[i]) for i in range(j, j + 3)]) for j in (2, 5, 8)]
    elif name == "SPC":
        i = 1
        while i + 1 < len(f):
            g = ival(f[i])
            if g > 0:
                spc.append((g, f[i + 1].strip()))
            i += 3

aero = defaultdict(lambda: np.zeros(3))
for name, f in read_cards(fem / f"uCRM-9_aero_forces_1.0g_{level}.bdf"):
    if name == "FORCE":
        if ival(f[2]) != 0:
            sys.exit("ABORT: FORCE with non-zero CID found")
        aero[ival(f[1])] += num(f[3]) * np.array([num(f[4]), num(f[5]), num(f[6])])


def axes(cid):
    A, B, C = cord[cid]
    z = (B - A) / np.linalg.norm(B - A)
    x = (C - A) - np.dot(C - A, z) * z
    x /= np.linalg.norm(x)
    return x, np.cross(z, x)


pid_mcid = defaultdict(set)
for eid, pid, ns, mcid in quads:
    pid_mcid[pid].add(mcid)
    if pid not in pshell or any(n not in grids for n in ns):
        sys.exit(f"ABORT: element {eid} references a missing property or node")
if any(len(v) != 1 for v in pid_mcid.values()):
    sys.exit("ABORT: a property uses more than one MCID")
for g, _ in spc:
    if g not in grids:
        sys.exit(f"ABORT: SPC node {g} missing")
for g in aero:
    if g not in grids:
        sys.exit(f"ABORT: loaded node {g} missing")

# --- extra loads for v2 -----------------------------------------------------------------------
shell_w = defaultdict(lambda: np.zeros(3))
mass_f = defaultdict(lambda: np.zeros(3))
groups = {}
m_shell = 0.0


def y_moment(fd):
    return sum(grids[n][1] * F[2] for n, F in fd.items())


if all_loads:
    for eid, pid, ns, _ in quads:
        p = [np.array(grids[n]) for n in ns]
        area = 0.5 * np.linalg.norm(np.cross(p[2] - p[0], p[3] - p[1]))
        pr = pshell[pid]
        me = area * (pr["t"] * mat2[pr["mid1"]]["rho"] + pr["nsm"])
        m_shell += me
        for n in ns:
            shell_w[n] += me * G_ACC / 4.0

    fuel_file = f"uCRM-9_fuel_mass_{'20_percent' if args.fuel == '20' else 'full'}_{level}.bdf"
    for label, fn in (("engine", f"uCRM-9_engine_mass_{level}.bdf"),
                      ("LE/TE", f"uCRM-9_LETE_masses_{level}.bdf"),
                      ("fuel", fuel_file)):
        mgrid, conm2, rbe2m, rbe3 = {}, {}, {}, []
        for name, f in read_cards(fem / fn):
            if name == "GRID":
                if ival(f[1]) != 0 or ival(f[5]) != 0:
                    sys.exit(f"ABORT: mass GRID with non-zero CP/CD in {fn}")
                mgrid[ival(f[0])] = np.array([num(f[2]), num(f[3]), num(f[4])])
            elif name == "CONM2":
                if ival(f[2]) != 0 or any(num(x) != 0 for x in f[4:7]) or any(num(x) != 0 for x in f[8:14]):
                    sys.exit(f"ABORT: CONM2 with CID, offset or inertia in {fn}")
                conm2[ival(f[1])] = conm2.get(ival(f[1]), 0.0) + num(f[3])
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

        refload = defaultdict(lambda: [np.zeros(3), np.zeros(3)])
        for g, m in conm2.items():
            F = m * G_ACC
            if g in rbe2m:
                mst = rbe2m[g]
                refload[mst][0] += F
                refload[mst][1] += np.cross(mgrid[g] - mgrid[mst], F)
            else:
                refload[g][0] += F
        if set(refload) != {r for r, _, _ in rbe3}:
            sys.exit(f"ABORT: reference nodes of masses and RBE3 cards disagree in {fn}")

        gf = defaultdict(lambda: np.zeros(3))
        for ref, w, deps in rbe3:
            if any(n not in grids for n in deps):
                sys.exit(f"ABORT: RBE3 dependent node missing from the wingbox in {fn}")
            F, M = refload[ref]
            d = np.array([grids[n] for n in deps]) - mgrid[ref]
            W, c = w * len(d), w * d.sum(axis=0)
            J = w * (np.sum(d * d) * np.eye(3) - d.T @ d)
            A = J - (np.dot(c, c) * np.eye(3) - np.outer(c, c)) / W
            if np.linalg.cond(A) > 1e10:
                sys.exit(f"ABORT: RBE3 {ref} dependents are degenerate (cond {np.linalg.cond(A):.2e})")
            b = np.linalg.solve(A, M - np.cross(c, F) / W)
            a = (F - np.cross(b, c)) / W
            f_i = w * (a + np.cross(b, d))
            tolF = 1e-9 * max(1.0, np.linalg.norm(F))
            tolM = 1e-9 * max(1.0, np.linalg.norm(M), np.linalg.norm(F) * np.linalg.norm(d, axis=1).max())
            if np.linalg.norm(f_i.sum(axis=0) - F) > tolF or np.linalg.norm(np.cross(d, f_i).sum(axis=0) - M) > tolM:
                sys.exit(f"ABORT: RBE3 {ref}: force/moment equivalence check failed")
            for n, fi in zip(deps, f_i):
                gf[n] += fi
                mass_f[n] += fi
        groups[label] = (sum(conm2.values()), gf)

L = []
L += [f"** uCRM-9 {level} wingbox, CalculiX deck {tag} ({'all loads' if all_loads else 'aero nodal forces only'})",
      "** Generated by scripts/ucrm_to_ccx.py. SI units: m, N, kg, Pa."]
L += ["*NODE, NSET=NALL"]
L += [f"{n}, {x:.9e}, {y:.9e}, {z:.9e}" for n, (x, y, z) in sorted(grids.items())]

orient_name, pid_orient, orient_axes = {}, {}, {}
for pid in sorted(pshell):
    x, y = axes(next(iter(pid_mcid[pid])))
    key = tuple(np.round(np.concatenate([x, y]), 6))
    if key not in orient_name:
        orient_name[key] = f"OR{len(orient_name) + 1}"
        orient_axes[orient_name[key]] = (x, y)
    pid_orient[pid] = orient_name[key]

by_pid = defaultdict(list)
for eid, pid, ns, _ in quads:
    by_pid[pid].append((eid, ns))
for pid in sorted(by_pid):
    L += [f"*ELEMENT, TYPE=S4, ELSET=P{pid}"]
    L += [f"{eid}, {ns[0]}, {ns[1]}, {ns[2]}, {ns[3]}" for eid, ns in by_pid[pid]]

for nm, (x, y) in orient_axes.items():
    L += [f"*ORIENTATION, NAME={nm}, SYSTEM=RECTANGULAR",
          f"{x[0]:.10f}, {x[1]:.10f}, {x[2]:.10f}, {y[0]:.10f}, {y[1]:.10f}, {y[2]:.10f}"]

n_bad = 0
for pid in sorted(pshell):
    p = pshell[pid]
    m1, m3 = mat2[p["mid1"]], mat2[p["mid3"]]
    G = np.array([[m1["g11"], m1["g12"]], [m1["g12"], m1["g22"]]])
    if m1["g11"] <= 0 or np.linalg.det(G) <= 0 or m3["g11"] <= 0 or m3["g22"] <= 0:
        n_bad += 1
        continue
    S = np.linalg.inv(G)
    E1, E2, nu12 = 1.0 / S[0, 0], 1.0 / S[1, 1], -S[0, 1] / S[0, 0]
    L += [f"*MATERIAL, NAME=M{pid}", "*ELASTIC, TYPE=ENGINEERING CONSTANTS",
          f"{E1:.9e}, {E2:.9e}, {args.e3_factor * E2:.9e}, {nu12:.6f}, {args.nu3:.6f}, {args.nu3:.6f}, "
          f"{m1['g33']:.9e}, {args.gshear_scale * m3['g11']:.9e}",
          f"{args.gshear_scale * m3['g22']:.9e}, 0.", "*DENSITY", f"{m1['rho']:.6f}"]
if n_bad:
    sys.exit(f"ABORT: {n_bad} properties with invalid material matrices")

for pid in sorted(pshell):
    L += [f"*SHELL SECTION, ELSET=P{pid}, MATERIAL=M{pid}, ORIENTATION={pid_orient[pid]}",
          f"{pshell[pid]['t']:.9e}"]

dofs = defaultdict(set)
for g, c in spc:
    for ch in c:
        if ch in "123456":
            dofs[g].add(int(ch))
spc_nodes = sorted(dofs)
L += ["*NSET, NSET=NSPC"]
for k in range(0, len(spc_nodes), 10):
    L.append(", ".join(str(n) for n in spc_nodes[k:k + 10]))

forces = defaultdict(lambda: np.zeros(3))
for src in (aero, shell_w, mass_f):
    for g, F in src.items():
        forces[g] += F

L += ["*STEP", "*STATIC", "*BOUNDARY"]
for g in spc_nodes:
    d = sorted(dofs[g])
    L += [f"{g}, 1, 6"] if d == [1, 2, 3, 4, 5, 6] else [f"{g}, {k}, {k}" for k in d]

L += ["*CLOAD"]
tot = np.zeros(3)
for g in sorted(forces):
    for k in range(3):
        if forces[g][k] != 0.0:
            L.append(f"{g}, {k + 1}, {forces[g][k]:.9e}")
    tot += forces[g]
L += ["*NODE FILE", "U", "*NODE PRINT, NSET=NALL", "U",
      "*NODE PRINT, NSET=NSPC, TOTALS=ONLY", "RF", "*END STEP"]

out_path.write_text("\n".join(L) + "\n")
print(f"wrote {out_path}")
print(f"settings: E3 factor {args.e3_factor}, nu13 = nu23 = {args.nu3}, shear scale {args.gshear_scale}")
print(f"nodes {len(grids)}  S4 elements {len(quads)}  properties {len(pshell)}  orientations {len(orient_axes)}")
print(f"constrained nodes {len(spc_nodes)}  fully clamped {sum(1 for g in spc_nodes if len(dofs[g]) == 6)}")
a_tot = sum(aero.values(), np.zeros(3))
print(f"\naero    : {len(aero):6d} nodes  F = ({a_tot[0]:+.4e}, {a_tot[1]:+.4e}, {a_tot[2]:+.4e}) N   "
      f"root bending moment sum(y Fz) = {y_moment(aero):+.4e} N m")
if all_loads:
    s_tot = sum(shell_w.values(), np.zeros(3))
    print(f"shells  : mass {m_shell:10.1f} kg  F = ({s_tot[0]:+.4e}, {s_tot[1]:+.4e}, {s_tot[2]:+.4e}) N   "
          f"sum(y Fz) = {y_moment(shell_w):+.4e} N m")
    inertial_m = y_moment(shell_w)
    for label, (m, gf) in groups.items():
        t = sum(gf.values(), np.zeros(3))
        print(f"{label:8s}: mass {m:10.1f} kg  F = ({t[0]:+.4e}, {t[1]:+.4e}, {t[2]:+.4e}) N   "
              f"sum(y Fz) = {y_moment(gf):+.4e} N m   (expected weight {m * 9.81:.1f} N)")
        inertial_m += y_moment(gf)
    ratio = (y_moment(aero) + inertial_m) / y_moment(aero)
    print(f"\nnet root bending moment / aero-only = {ratio:.3f}   "
          f"(crude scale factor for the v1 tip deflection; load shapes differ)")
print(f"TOTAL applied: Fx={tot[0]:+.4e}  Fy={tot[1]:+.4e}  Fz={tot[2]:+.4e} N   loaded nodes {len(forces)}")
