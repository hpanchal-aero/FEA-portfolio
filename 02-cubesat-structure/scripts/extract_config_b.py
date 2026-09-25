"""
Project 02, Config B -- extract results and net-section stress.

Usage: /usr/bin/python3 scripts/extract_config_b.py <job_path_without_ext>

Net-section stress (pragmatic proxy, not an integrated force/area value):
mean SZZ (signed axial stress, not von Mises) over panel-region corner
nodes at z=50 +/- 1.5mm (the pocket throat for all current shapes), split
into "full section" (all such nodes) and "excluding top decile" (a corner-
singularity-exclusion proxy, since re-entrant corners also spike SZZ, not
just von Mises). The gap between the two is reported as a corner-
sensitivity indicator, not a converged stress-concentration factor.
"""
import json
import os
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extract_config_a import read_frd_block, von_mises  # noqa: E402

PRESSURE = 0.1643
FORCE_APPLIED = 156.5779
RAIL_INNER = 41.5
TOL = 1e-3
SECTION_Z = 50.0
SECTION_TOL = 1.5


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


def main():
    if len(sys.argv) != 2:
        print("Usage: extract_config_b.py <job_path_without_ext>")
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
    szz = S[:, 2]
    z = coords[:, 2]

    bottom = np.nonzero(np.abs(z[1:] - 0.0) < TOL)[0] + 1
    top = np.nonzero(np.abs(z[1:] - 100.0) < TOL)[0] + 1

    print("\n--- Equilibrium ---")
    rb = RF[bottom].sum(axis=0)
    eq_err = 100 * (rb[2] - FORCE_APPLIED) / FORCE_APPLIED
    print(f"Bottom reaction Z: {rb[2]:.6f} N (error {eq_err:+.4f}%)")

    print("\n--- Overall axial stiffness ---")
    uz_top_mean = U[top, 2].mean()
    stiffness = FORCE_APPLIED / abs(uz_top_mean)
    print(f"Uz top mean: {uz_top_mean:.8f} mm  ->  effective stiffness F/|d| = {stiffness:.4f} N/mm")

    print(f"\n--- Net-section stress (SZZ, panel region, z={SECTION_Z}+/-{SECTION_TOL}mm) ---")
    rail = (np.abs(coords[:, 0]) >= RAIL_INNER - TOL) & (np.abs(coords[:, 1]) >= RAIL_INNER - TOL)
    section = (np.abs(z - SECTION_Z) < SECTION_TOL) & ~rail
    section[0] = False
    n_sec = int(section.sum())
    if n_sec == 0:
        print("WARNING: no panel nodes found at this section -- check SECTION_Z/geometry.")
        net_full = net_excl = float("nan")
    else:
        vals = szz[section]
        order = np.argsort(vals)  # most negative (compressive) = highest magnitude
        cutoff = max(1, int(0.1 * n_sec))
        excl_idx = order[cutoff:-cutoff] if n_sec > 2 * cutoff else order
        net_full = vals.mean()
        net_excl = vals[excl_idx].mean() if len(excl_idx) else float("nan")
        gap_pct = 100 * (net_full - net_excl) / abs(net_excl) if net_excl else float("nan")
        print(f"nodes={n_sec}  full-section mean SZZ={net_full:.6f} MPa  "
              f"(vs applied {-PRESSURE} MPa, %diff={100*(net_full-(-PRESSURE))/PRESSURE:+.2f}%)")
        print(f"excl-top/bottom-decile mean SZZ={net_excl:.6f} MPa  "
              f"corner-sensitivity gap={gap_pct:+.2f}%")

    print("\n--- Peak von Mises (10 highest nodes) ---")
    order_vm = np.argsort(vm[1:])[::-1][:10] + 1
    for nd in order_vm[:5]:
        c = coords[nd]
        print(f"  node {nd:6d}: VM={vm[nd]:.6f}  at ({c[0]:8.3f}, {c[1]:8.3f}, {c[2]:8.3f})")

    os.makedirs("results", exist_ok=True)
    out = f"results/{os.path.basename(job)}_configB_results.json"
    res = {
        "job": job, "nodes": n, "elements": n_elem,
        "reaction_z_bottom_N": float(rb[2]), "equilibrium_err_pct": float(eq_err),
        "uz_top_mean_mm": float(uz_top_mean), "stiffness_N_per_mm": float(stiffness),
        "net_section_szz_full_MPa": float(net_full), "net_section_szz_excl_MPa": float(net_excl),
        "peak_vm_MPa": float(vm[order_vm[0]]),
    }
    with open(out, "w") as f:
        json.dump(res, f, indent=2)
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
