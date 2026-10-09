#!/usr/bin/env python3
"""
Project 03, Stage 3 -- extract a spar v2 + ribs (C3D10) + skin (S6) static solve.

Load case is read from the '** LOADCASE LIFT_N=.. NG=..' header of the .inp (default 24.5 N, n=1).
Weights (solid, skin, rib point weights) are scaled by NG.

Verification: tip deflection vs tapered-beam references (spar v2 only "B_lo" and spar v2 + fully
effective skin ring "B_hi", with the same loads incl. skin+rib weight, ribs as point weights).
B_hi is APPROXIMATE (vertical-offset ring, no Poisson/ribs) -- a plausibility bracket, not a bound.
Reaction-force equilibrium is NOT evaluated: FORC / *NODE PRINT totals are not interpretable for
CalculiX expanded shells (see development log).

Stress: von Mises on solid-only nodes (original IDs not used by any shell) and on CalculiX's expanded
shell nodes (IDs > number of *NODE entries) -- the expanded-node interpretation is ASSUMED, not verified.
Peaks near the fixed root are constraint-corner values, so peaks with z > 10 mm are reported separately.
Yield assumed 290 MPa (low end of typical 2024-T3 range; UNCITED).

Usage (system Python):
  /usr/bin/python3 scripts/extract_stage3_skin.py <job_path_without_ext>
"""
import json
import math
import os
import re
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_stage3_solid as b
from extract_stage3_solid import A_of_z, I_of_z, ZTAB, von_mises
from extract_stage3_ribs import RIBW

SPAN = 1000.0
E_MOD = 73100.0
RHO = 2780e-12
G = 9810.0
TOTAL_LIFT_N = 24.5            # default only; real value from the .inp header
YIELD_ASSUMED = 290.0          # MPa, UNCITED


def read_loadcase(path):
    lift, ng = 24.5, 1.0
    with open(path) as f:
        for _ in range(12):
            line = f.readline()
            m = re.search(r"LOADCASE\s+LIFT_N=([-+0-9.eE]+)\s+NG=([-+0-9.eE]+)", line)
            if m:
                lift, ng = float(m.group(1)), float(m.group(2))
                break
    return lift, ng


def parse_inp(path):
    xyz, solid, shell, root, cload = [], [], [], [], []
    mode, t_skin = None, None
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("*"):
                u = s.upper()
                if u.startswith("*NODE") and not u.startswith("*NODE FILE") and not u.startswith("*NODE PRINT"):
                    mode = "node"
                elif u.startswith("*ELEMENT"):
                    mode = "solid" if "C3D10" in u else ("shell" if "S6" in u else None)
                elif u.startswith("*NSET") and "NROOT" in u:
                    mode = "root"
                elif u.startswith("*CLOAD"):
                    mode = "cload"
                elif u.startswith("*SHELL SECTION"):
                    mode = "thick"
                else:
                    mode = None
                continue
            p = [t.strip() for t in s.split(",") if t.strip()]
            if mode == "node":
                xyz.append([float(v) for v in p[1:4]])
            elif mode == "solid":
                solid.append([int(v) for v in p[1:11]])
            elif mode == "shell":
                shell.append([int(v) for v in p[1:7]])
            elif mode == "root":
                root.extend(int(v) for v in p)
            elif mode == "cload":
                cload.append((int(p[0]), int(p[1]), float(p[2])))
            elif mode == "thick":
                t_skin = float(p[0]); mode = None
    return (np.array(xyz), np.array(solid), np.array(shell), np.array(root), cload, t_skin)


def read_frd(path, name, ncomp):
    ids, vals = [], []
    in_block = False
    with open(path) as f:
        for line in f:
            if line.startswith(" -4"):
                in_block = (line[5:13].strip() == name)
                continue
            if in_block:
                if line.startswith(" -3"):
                    in_block = False
                elif line.startswith(" -1"):
                    ids.append(int(line[3:13]))
                    vals.append([float(line[13 + 12*k:25 + 12*k]) for k in range(ncomp)])
    if not ids:
        print(f"ERROR: block {name} not found"); sys.exit(1)
    ids = np.array(ids)
    arr = np.zeros((ids.max() + 1, ncomp))
    arr[ids] = np.array(vals)
    return arr


def read_frd_nodes(path):
    ids, xyz = [], []
    on = False
    with open(path) as f:
        for line in f:
            if line.startswith("    2C"):
                on = True
                continue
            if on:
                if line.startswith(" -3"):
                    break
                if line.startswith(" -1"):
                    ids.append(int(line[3:13]))
                    xyz.append([float(line[13 + 12*k:25 + 12*k]) for k in range(3)])
    return np.array(ids), np.array(xyz)


def ring_props(z, t):
    out = b.outline_at(z)
    n = len(out)
    perim, I = 0.0, 0.0
    for i in range(n):
        x0, y0 = out[i]
        x1, y1 = out[(i + 1) % n]
        s = math.hypot(x1 - x0, y1 - y0)
        ym0 = y0 + math.copysign(t / 2.0, y0) if y0 != 0 else 0.0
        ym1 = y1 + math.copysign(t / 2.0, y1) if y1 != 0 else 0.0
        perim += s
        I += t * s / 3.0 * (ym0 ** 2 + ym0 * ym1 + ym1 ** 2)
    return perim, I


def beam_tip(t, include_skin_stiffness, lift, ng):
    ptab = np.array([ring_props(z, t) for z in ZTAB])
    perim_tab, Iring_tab = ptab[:, 0], ptab[:, 1]
    n = 20000
    z = np.linspace(0, SPAN, n)
    dz = z[1] - z[0]
    q0 = lift / (np.pi * SPAN / 4.0)
    A_skin = np.interp(z, ZTAB, perim_tab) * t
    q = q0 * np.sqrt(np.clip(1.0 - (z / SPAN) ** 2, 0, None)) - RHO * G * ng * (A_of_z(z) + A_skin)
    V = np.zeros(n)
    for i in range(n - 2, -1, -1):
        V[i] = V[i + 1] + 0.5 * (q[i] + q[i + 1]) * dz
    for zr, _, _, P in RIBW:
        V[z < zr - 1e-9] -= P * ng
    M = np.zeros(n)
    for i in range(n - 2, -1, -1):
        M[i] = M[i + 1] + 0.5 * (V[i] + V[i + 1]) * dz
    EI = E_MOD * (I_of_z(z) + (np.interp(z, ZTAB, Iring_tab) if include_skin_stiffness else 0.0))
    curv = M / EI
    th = np.zeros(n); v = np.zeros(n)
    for i in range(1, n):
        th[i] = th[i - 1] + 0.5 * (curv[i - 1] + curv[i]) * dz
    for i in range(1, n):
        v[i] = v[i - 1] + 0.5 * (th[i - 1] + th[i]) * dz
    return v[-1]


def main():
    if len(sys.argv) != 2:
        print("Usage: extract_stage3_skin.py <job>"); sys.exit(1)
    job = sys.argv[1]
    import re as _re_rs
    _m_rs = _re_rs.search(r"_rs(\d+)(?:_|$)", os.path.basename(job))
    _want_rs = float(_m_rs.group(1)) if _m_rs else 200.0
    if b.RIB_SPACING != _want_rs:
        print(f"ERROR: job name implies RIB_SPACING={_want_rs:g} but the environment gives {b.RIB_SPACING:g}; set RIB_SPACING={_want_rs:g} -- STOP."); sys.exit(1)
    inp, frd = job + ".inp", job + ".frd"
    lift, ng = read_loadcase(inp)
    xyz, solid, shell, root, cload, t_skin = parse_inp(inp)
    N0 = len(xyz)
    coords = np.vstack([np.zeros((1, 3)), xyz])
    print(f"Job: {job}\nLoad case: lift {lift:.4f} N per half-wing, load factor n = {ng}")
    print(f"Original nodes: {N0}; C3D10: {len(solid)}; S6: {len(shell)}; skin t = {t_skin} mm")

    U = read_frd(frd, "DISP", 3)
    S = read_frd(frd, "STRESS", 6)
    fid, fxyz = read_frd_nodes(frd)
    print(f"frd nodes: {len(fid)} (max id {fid.max()}); original IDs present {int((fid <= N0).sum())}, "
          f"expanded {int((fid > N0).sum())}")

    # ---- masses / loads ----
    P = coords[solid[:, :4]]
    vol = np.abs(np.einsum("ij,ij->i", np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]),
                           P[:, 3] - P[:, 0])) / 6.0
    Q = coords[shell[:, :3]]
    area = 0.5 * np.linalg.norm(np.cross(Q[:, 1] - Q[:, 0], Q[:, 2] - Q[:, 0]), axis=1)
    skin_area = area.sum()
    W_solid = RHO * G * ng * vol.sum()
    W_skin = RHO * G * ng * t_skin * skin_area
    cload_y = sum(f for nd, dof, f in cload if dof == 2)
    print(f"\nSolid volume {vol.sum():.4f} mm^3 -> weight at n={ng}: {W_solid:.4f} N "
          f"(mass {RHO * vol.sum() * 1e3:.4f} kg)")
    print(f"Skin area {skin_area:.4f} mm^2 (geometry 459622.6823), t={t_skin} -> weight {W_skin:.4f} N "
          f"(mass {RHO * t_skin * skin_area * 1e3:.4f} kg)")
    print(f"Applied lift: CLOAD sum {cload_y:.4f} N (header {lift:.4f}); net up {lift - W_solid - W_skin:.4f} N")
    print("Reaction-force equilibrium: NOT evaluated (FORC / NODE PRINT totals not interpretable for "
          "expanded shells; see development log).")

    # ---- tip deflection on solid nodes ----
    solid_ids = np.unique(solid)
    z = coords[:, 2]
    tipsel = solid_ids[(np.abs(z[solid_ids] - SPAN) < 2.0) & (np.abs(coords[solid_ids, 0]) < 4.5)]
    uy = U[tipsel, 1]
    if uy.min() == 0.0:
        print("WARNING: zero displacement row among tip nodes")
    uy_axis = uy.mean()
    v_lo = beam_tip(t_skin, False, lift, ng)
    v_hi = beam_tip(t_skin, True, lift, ng)
    frac = (v_lo - uy_axis) / (v_lo - v_hi) if v_lo != v_hi else float("nan")
    print("\n--- Tip deflection (Y), spar-axis solid nodes ---")
    print(f"FEA Uy = {uy_axis:.6f} mm (n={len(tipsel)}, spread {uy.min():.6f}..{uy.max():.6f})")
    print(f"Beam B_lo (spar only, skin+ribs as weight)      : {v_lo:.6f} mm")
    print(f"Beam B_hi (spar + fully effective skin ring)    : {v_hi:.6f} mm")
    print(f"FEA / B_hi - 1 = {100 * (uy_axis / v_hi - 1):+.2f}%;  skin effectiveness fraction "
          f"(B_lo -> B_hi) = {frac:.3f}")

    # ---- stress ----
    vm = von_mises(S)
    skin_ids = np.unique(shell)
    solid_only = np.setdiff1d(solid_ids, skin_ids)
    zs = coords[solid_only, 2]

    def report(label, ids, zarr, xyz_arr):
        for tag, m in (("all", np.ones(len(ids), bool)), ("z>10mm", zarr > 10.0)):
            if not m.any():
                continue
            v = vm[ids[m]]
            k = int(np.argmax(v))
            c = xyz_arr[m][k]
            print(f"  {label:12s} {tag:7s}: max VM {v.max():8.3f} MPa at ({c[0]:7.2f},{c[1]:7.2f},{c[2]:8.2f}); "
                  f"p99.9 {np.percentile(v, 99.9):8.3f}, p99 {np.percentile(v, 99):8.3f}; "
                  f"FoS(yield {YIELD_ASSUMED:.0f}) at max = {YIELD_ASSUMED / v.max():.2f}")

    print("\n--- Von Mises (local peaks; convergence NOT established; yield assumed, UNCITED) ---")
    report("solid-only", solid_only, zs, coords[solid_only])
    exp_mask = fid > N0
    exp_ids, exp_xyz = fid[exp_mask], fxyz[exp_mask]
    report("skin (exp.)", exp_ids, exp_xyz[:, 2], exp_xyz)
    print("(skin stress = CalculiX expanded-node values; interpretation assumed, not verified)")

    os.makedirs("results", exist_ok=True)
    out = f"results/{os.path.basename(job)}_results.json"
    with open(out, "w") as f:
        json.dump({
            "job": job, "lift_N": lift, "n": ng, "t_skin_mm": t_skin,
            "W_solid_N": float(W_solid), "W_skin_N": float(W_skin),
            "uy_spar_axis_mm": float(uy_axis), "bracket_lo_mm": float(v_lo),
            "bracket_hi_mm": float(v_hi), "skin_effectiveness_fraction": float(frac),
            "vm_skin_max_MPa": float(vm[exp_ids].max()),
            "vm_skin_max_z_gt_10_MPa": float(vm[exp_ids[exp_xyz[:, 2] > 10.0]].max()),
            "vm_solid_only_max_MPa": float(vm[solid_only].max()),
            "vm_solid_only_max_z_gt_10_MPa": float(vm[solid_only[zs > 10.0]].max()),
        }, f, indent=2)
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
