"""
Project 02, Config A -- extract CalculiX results and check against theory.

Usage:  /usr/bin/python3 scripts/extract_config_a.py <job_path_without_ext>

Reads coordinates from <job>.inp and DISP / FORC / STRESS blocks from <job>.frd
(fixed-width: node id at line[3:13], values in 12-char fields from col 13;
STRESS has 6 components on one -1 line: SXX SYY SZZ SXY SYZ SZX).

Equilibrium note: CalculiX RF is the internal nodal force at every node. The
reaction of the fixed base is the sum of RF over the constrained bottom nodes
(z=0). Loaded top nodes carry the opposite sign, so the sum over ALL nodes is
zero by construction and is NOT an equilibrium check.

Reference (uniaxial, prismatic assembly, uniform top pressure, fixed base):
    sigma = P,  delta = P*L/E.  Verifies axial stiffness and mean stress only.
"""

import json
import os
import sys
import numpy as np

E_MOD = 68900.0
PRESSURE = 0.1643
L = 100.0
TOP_AREA = 953.0
FORCE_APPLIED = PRESSURE * TOP_AREA
SIGMA_THEORY = PRESSURE
DELTA_THEORY = PRESSURE * L / E_MOD
RAIL_INNER = 41.5
TOL = 1e-3


def read_inp_nodes(path):
    xyz, n_elem, mode = [], 0, None
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("*"):
                u = s.upper()
                if u.startswith("*NODE") and not u.startswith("*NODE FILE"):
                    mode = "node"
                elif u.startswith("*ELEMENT"):
                    mode = "elem"
                else:
                    mode = None
                continue
            if mode == "node":
                p = s.split(",")
                xyz.append([float(v) for v in p[1:4]])
            elif mode == "elem":
                n_elem += 1
    return np.array(xyz), n_elem


def read_frd_block(path, name, ncomp, n_nodes):
    """Return array (n_nodes+1, ncomp); last block with this name wins."""
    data = None
    in_block = False
    with open(path) as f:
        for line in f:
            if line.startswith(" -4"):
                in_block = (line[5:13].strip() == name)
                if in_block:
                    data = np.zeros((n_nodes + 1, ncomp))
                continue
            if in_block:
                if line.startswith(" -3"):
                    in_block = False
                elif line.startswith(" -1"):
                    node = int(line[3:13])
                    for k in range(ncomp):
                        data[node, k] = float(line[13 + 12 * k: 25 + 12 * k])
    if data is None:
        print(f"ERROR: block {name} not found in {path}")
        sys.exit(1)
    return data


def von_mises(s):
    sxx, syy, szz, sxy, syz, szx = (s[:, i] for i in range(6))
    return np.sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2)
                   + 3.0 * (sxy ** 2 + syz ** 2 + szx ** 2))


def main():
    if len(sys.argv) != 2:
        print("Usage: extract_config_a.py <job_path_without_ext>")
        sys.exit(1)
    job = sys.argv[1]
    inp, frd = job + ".inp", job + ".frd"
    xyz, n_elem = read_inp_nodes(inp)
    n = len(xyz)
    coords = np.vstack([np.zeros((1, 3)), xyz])
    print(f"Job: {job}\nNodes: {n}, elements: {n_elem}")

    U = read_frd_block(frd, "DISP", 3, n)
    RF = read_frd_block(frd, "FORC", 3, n)
    S = read_frd_block(frd, "STRESS", 6, n)
    vm = von_mises(S)
    z = coords[:, 2]

    bottom = np.nonzero(np.abs(z[1:] - 0.0) < TOL)[0] + 1
    top = np.nonzero(np.abs(z[1:] - L) < TOL)[0] + 1

    print("\n--- Equilibrium ---")
    rb = RF[bottom].sum(axis=0)
    rt = RF[top].sum(axis=0)
    ra = RF[1:].sum(axis=0)
    eq_err = 100 * (rb[2] - FORCE_APPLIED) / FORCE_APPLIED
    print(f"Bottom nodes (z=0): {len(bottom)}   Top nodes (z={L:.0f}): {len(top)}")
    print(f"Sum RF, bottom (reaction): X={rb[0]:.6f}  Y={rb[1]:.6f}  Z={rb[2]:.6f} N")
    print(f"Sum RF, top:               X={rt[0]:.6f}  Y={rt[1]:.6f}  Z={rt[2]:.6f} N")
    print(f"Sum RF, all nodes:         X={ra[0]:.6f}  Y={ra[1]:.6f}  Z={ra[2]:.6f} N  (should be ~0)")
    print(f"Applied: {FORCE_APPLIED:.6f} N (compressive, -z); expected bottom reaction +{FORCE_APPLIED:.6f} N")
    print(f"Z equilibrium error (bottom): {eq_err:+.4f}%")

    print("\n--- Axial displacement at top (z=100) ---")
    uz = U[top, 2]
    uz_pct = 100 * (uz.mean() - (-DELTA_THEORY)) / (-DELTA_THEORY)
    print(f"Top nodes: {len(top)}")
    print(f"Uz mean={uz.mean():.8f}  min={uz.min():.8f}  max={uz.max():.8f} mm")
    print(f"Theory Uz = {-DELTA_THEORY:.8f} mm -> mean %diff = {uz_pct:+.4f}%")
    print(f"Max |Ux|={np.abs(U[1:,0]).max():.3e}, max |Uy|={np.abs(U[1:,1]).max():.3e} mm")

    print("\n--- Mean von Mises, slab 25 < z < 75 mm ---")
    slab = (z > 25.0) & (z < 75.0)
    slab[0] = False
    rail = (np.abs(coords[:, 0]) >= RAIL_INNER - TOL) & (np.abs(coords[:, 1]) >= RAIL_INNER - TOL)
    m_all = slab
    m_rail = slab & rail
    m_pan = slab & ~rail
    vm_all, vm_rail, vm_pan = vm[m_all].mean(), vm[m_rail].mean(), vm[m_pan].mean()
    for label, m, v in (("all", m_all, vm_all), ("rail region", m_rail, vm_rail),
                        ("panel region", m_pan, vm_pan)):
        print(f"{label:13s}: nodes={int(m.sum()):6d}  mean VM={v:.6f} MPa  "
              f"%diff vs {SIGMA_THEORY} = {100*(v-SIGMA_THEORY)/SIGMA_THEORY:+.3f}%")

    print("\n--- Peak von Mises (10 highest nodes) ---")
    order = np.argsort(vm[1:])[::-1][:10] + 1
    print(f"Peak VM = {vm[order[0]]:.6f} MPa  (mean-slab ratio {vm[order[0]]/vm_all:.2f}x)")
    for nd in order:
        c = coords[nd]
        print(f"  node {nd:6d}: VM={vm[nd]:.6f}  at ({c[0]:8.3f}, {c[1]:8.3f}, {c[2]:8.3f})")

    os.makedirs("results", exist_ok=True)
    out = f"results/{os.path.basename(job)}_results.json"
    res = {
        "job": job, "nodes": n, "elements": n_elem,
        "equations": 3 * n,
        "reaction_z_bottom_N": float(rb[2]), "applied_N": FORCE_APPLIED,
        "equilibrium_err_pct": float(eq_err),
        "rf_z_top_N": float(rt[2]), "rf_z_all_N": float(ra[2]),
        "uz_top_mean_mm": float(uz.mean()), "uz_theory_mm": -DELTA_THEORY,
        "uz_pct_diff": float(uz_pct),
        "mean_vm_slab_MPa": float(vm_all),
        "mean_vm_rail_MPa": float(vm_rail), "mean_vm_panel_MPa": float(vm_pan),
        "peak_vm_MPa": float(vm[order[0]]),
        "peak_vm_location": [float(v) for v in coords[order[0]]],
    }
    with open(out, "w") as f:
        json.dump(res, f, indent=2)
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
