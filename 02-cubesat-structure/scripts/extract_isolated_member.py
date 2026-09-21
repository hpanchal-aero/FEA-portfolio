"""
Project 02, Sub-Stages 2a / 2b -- re-extraction of stored CalculiX results.

Usage:  /usr/bin/python3 scripts/extract_isolated_member.py <a|b>

  a : isolated 83 x 2 x 100 mm panel   (simulation/stage2a_panel_solve.*)
  b : isolated 8.5 x 8.5 x 100 mm rail (simulation/stage2b_rail_solve.*)

Read-only: parses <job>.inp and <job>.frd. Reference is uniaxial (axial-bar)
theory for a prismatic member under uniform top pressure with a fixed base:
    sigma = P,  delta = P*L/E.
Also recomputes the closed-form buckling checks:
    2a: plate, k = 4 (simply supported), b = 83 mm, t = 2 mm
    2b: Euler column, K = 2 (fixed-free)
Neither is an FEA buckling analysis.

The mean von Mises here uses the slab 25 < z < 75 mm (same definition as
extract_config_a.py). The definition used when 2a/2b were first run is not
known, so differences from those earlier values are printed, not reconciled.
"""

import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extract_config_a import read_frd_block, von_mises  # noqa: E402

E_MOD = 68900.0
NU = 0.33
PRESSURE = 0.1643
L = 100.0
TOL = 1e-3

MEMBERS = {
    "a": dict(
        label="Sub-stage 2a: isolated panel 83 x 2 x 100 mm",
        job="simulation/stage2a_panel_solve",
        area=83.0 * 2.0,
        out="results/stage2a_panel_results.json",
        prior=dict(nodes=81735, elements=47888, equilibrium_pct=0.0000,
                   mean_vm=0.161993, uz_top=-0.00023537, peak_vm=0.345026),
    ),
    "b": dict(
        label="Sub-stage 2b: isolated rail 8.5 x 8.5 x 100 mm",
        job="simulation/stage2b_rail_solve",
        area=8.5 * 8.5,
        out="results/stage2b_rail_results.json",
        prior=dict(nodes=18755, elements=11398, equilibrium_pct=0.0002,
                   mean_vm=0.163368, uz_top=-0.00023767, peak_vm=0.259284),
    ),
}


def read_inp(path):
    tags, xyz, n_elem, mode = [], [], 0, None
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
                tags.append(int(p[0]))
                xyz.append([float(v) for v in p[1:4]])
            elif mode == "elem":
                n_elem += 1
    tags = np.array(tags)
    if not np.array_equal(tags, np.arange(1, len(tags) + 1)):
        print("ERROR: node tags in .inp are not contiguous 1..N.")
        sys.exit(1)
    return np.array(xyz), n_elem


def buckling(member):
    if member == "a":
        b, t, k = 83.0, 2.0, 4.0
        sigma_cr = k * math.pi ** 2 * E_MOD / (12.0 * (1.0 - NU ** 2)) * (t / b) ** 2
        p_cr = sigma_cr * b * t
        desc = f"plate, k={k}, b={b}, t={t}: sigma_cr={sigma_cr:.4f} MPa"
    else:
        side, kk = 8.5, 2.0
        inertia = side ** 4 / 12.0
        p_cr = math.pi ** 2 * E_MOD * inertia / (kk * L) ** 2
        desc = f"Euler, K={kk}, I={inertia:.4f} mm^4"
    return p_cr, desc


def main():
    if len(sys.argv) != 2 or sys.argv[1] not in MEMBERS:
        print("Usage: extract_isolated_member.py <a|b>")
        sys.exit(1)
    key = sys.argv[1]
    m = MEMBERS[key]
    job = m["job"]
    inp, frd = job + ".inp", job + ".frd"
    for p in (inp, frd):
        if not os.path.isfile(p):
            print(f"ERROR: {p} not found.")
            sys.exit(1)

    force = PRESSURE * m["area"]
    delta_theory = PRESSURE * L / E_MOD

    xyz, n_elem = read_inp(inp)
    n = len(xyz)
    coords = np.vstack([np.zeros((1, 3)), xyz])
    z = coords[:, 2]
    print("=" * 70)
    print(m["label"])
    print("=" * 70)
    print(f"Job: {job}\nNodes: {n}, elements: {n_elem}, equations: {3 * n}")
    print(f"Applied: P = {PRESSURE} MPa x A = {m['area']} mm^2 -> F = {force:.6f} N")
    print(f"Theory: sigma = {PRESSURE} MPa, delta = {delta_theory:.8f} mm")

    U = read_frd_block(frd, "DISP", 3, n)
    RF = read_frd_block(frd, "FORC", 3, n)
    S = read_frd_block(frd, "STRESS", 6, n)
    vm = von_mises(S)

    bottom = np.nonzero(np.abs(z[1:] - 0.0) < TOL)[0] + 1
    top = np.nonzero(np.abs(z[1:] - L) < TOL)[0] + 1

    print("\n--- Equilibrium ---")
    rb, rt, ra = RF[bottom].sum(axis=0), RF[top].sum(axis=0), RF[1:].sum(axis=0)
    eq_pct = 100.0 * (rb[2] - force) / force
    print(f"Bottom nodes: {len(bottom)}   Top nodes: {len(top)}")
    print(f"Sum RF bottom: X={rb[0]:.6f} Y={rb[1]:.6f} Z={rb[2]:.6f} N")
    print(f"Sum RF top:    X={rt[0]:.6f} Y={rt[1]:.6f} Z={rt[2]:.6f} N")
    print(f"Sum RF all:    X={ra[0]:.6f} Y={ra[1]:.6f} Z={ra[2]:.6f} N (should be ~0)")
    print(f"Z equilibrium error (bottom vs applied): {eq_pct:+.4f}%")

    print("\n--- Top-face axial displacement ---")
    uz = U[top, 2]
    uz_pct = 100.0 * (uz.mean() + delta_theory) / delta_theory
    print(f"Uz mean={uz.mean():.8f} min={uz.min():.8f} max={uz.max():.8f} mm")
    print(f"Theory Uz = {-delta_theory:.8f} mm -> mean %diff = {uz_pct:+.4f}%")

    print("\n--- Mean von Mises, slab 25 < z < 75 mm ---")
    slab = (z > 25.0) & (z < 75.0)
    slab[0] = False
    vm_mean = vm[slab].mean()
    vm_pct = 100.0 * (vm_mean - PRESSURE) / PRESSURE
    print(f"nodes={int(slab.sum())}  mean VM={vm_mean:.6f} MPa  %diff vs {PRESSURE} = {vm_pct:+.3f}%")

    print("\n--- Peak von Mises (10 highest nodes) ---")
    order = np.argsort(vm[1:])[::-1][:10] + 1
    print(f"Peak VM = {vm[order[0]]:.6f} MPa (ratio to slab mean {vm[order[0]] / vm_mean:.2f}x)")
    for nd in order:
        c = coords[nd]
        print(f"  node {nd:6d}: VM={vm[nd]:.6f}  at ({c[0]:8.3f}, {c[1]:8.3f}, {c[2]:8.3f})")
    n_base = int((coords[order, 2] < 3.0).sum())
    print(f"Top-10 nodes with z < 3 mm: {n_base} of 10")

    p_cr, desc = buckling(key)
    print("\n--- Closed-form buckling check (isolated member; not an FEA buckling analysis) ---")
    print(f"{desc}, P_cr = {p_cr:.2f} N; applied/critical = {force / p_cr:.4f}")

    pr = m["prior"]
    print("\n--- Comparison with EARLIER TRANSFERRED VALUES (unverified) ---")
    print(f"nodes:        re-extracted {n:7d} | earlier {pr['nodes']:7d}")
    print(f"elements:     re-extracted {n_elem:7d} | earlier {pr['elements']:7d}")
    print(f"equilibrium%: re-extracted {eq_pct:+8.4f} | earlier {pr['equilibrium_pct']:+8.4f}")
    print(f"mean VM:      re-extracted {vm_mean:8.6f} | earlier {pr['mean_vm']:8.6f}  (definitions may differ)")
    print(f"Uz top mean:  re-extracted {uz.mean():.8f} | earlier {pr['uz_top']:.8f}")
    print(f"peak VM:      re-extracted {vm[order[0]]:8.6f} | earlier {pr['peak_vm']:8.6f}")

    os.makedirs("results", exist_ok=True)
    res = {
        "job": job, "nodes": n, "elements": n_elem, "equations": 3 * n,
        "applied_N": force,
        "reaction_z_bottom_N": float(rb[2]), "rf_z_top_N": float(rt[2]),
        "rf_z_all_N": float(ra[2]), "equilibrium_err_pct": float(eq_pct),
        "uz_top_mean_mm": float(uz.mean()), "uz_theory_mm": -delta_theory,
        "uz_pct_diff": float(uz_pct),
        "mean_vm_slab_MPa": float(vm_mean), "mean_vm_slab_pct_diff": float(vm_pct),
        "peak_vm_MPa": float(vm[order[0]]),
        "peak_vm_location": [float(v) for v in coords[order[0]]],
        "top10_peak_nodes_below_3mm": n_base,
        "closed_form_P_cr_N": float(p_cr),
        "applied_over_critical": float(force / p_cr),
    }
    with open(m["out"], "w") as f:
        json.dump(res, f, indent=2)
    print(f"\nWrote {m['out']}")


if __name__ == "__main__":
    main()
