"""
Project 01 - Aerospace Mounting Bracket
Diagnostic investigation of the converged verification-plate results,
per structured investigation request (see chat log). Uses ONLY
already-saved result files - no FEA re-solving.

Addresses:
  1. Whether the ~5-8% analytical discrepancy is explained by
     beam-theory idealization vs 3D solid behavior (not assumed).
  2. Root-band extraction window sensitivity ([0,6],[1,7],[3,9],[5,11]).
  3. Root cause of non-monotonic stress convergence across mesh levels.
  4. Tet-vs-hex agreement - what it does/doesn't demonstrate (discussed
     in the accompanying chat response, not computed here).
"""

import os
import sys
import numpy as np

sys.path.insert(0, ".")
import hex_lib as lib

SIM_ROOT = "../simulation/convergence_hex"
FINEST_CASE = os.path.join(SIM_ROOT, "mesh_0.3mm")

L, B, H = lib.L, lib.B, lib.H
F = 235.44  # N
E = 71700.0  # MPa
I_SECTION = B * H**3 / 12.0  # mm^4
C = H / 2.0  # mm, distance to outer fiber

def sigma_theory(x):
    """Euler-Bernoulli bending stress at outer fiber, position x (mm)."""
    M = F * (L - x)
    return M * C / I_SECTION

print("="*70)
print("TASK 1: Beam-theory idealization vs 3D solid behavior")
print("="*70)

print(f"\nSection properties check: I = {I_SECTION:.4f} mm^4 "
      f"(expected 213.3333), c = {C} mm")
print(f"sigma_theory(x=0)  = {sigma_theory(0):.3f} MPa (target: 132.4)")
print(f"sigma_theory(x=3)  = {sigma_theory(3):.3f} MPa")
print(f"sigma_theory(x=6)  = {sigma_theory(6):.3f} MPa")
print(f"sigma_theory(x=9)  = {sigma_theory(9):.3f} MPa")
print("^ NOTE: these are all LOWER than 132.4 MPa, since bending moment")
print("  decreases with x. The 132.4 MPa target is a SINGLE POINT value")
print("  at x=0, which our root-band [3,9] deliberately excludes.")

print("\n--- 1a. Tip Uz: full-face average vs centroid-node vs corner-only ---")

clean_inp = os.path.join(FINEST_CASE, "mesh_clean.inp")
nset_inp = os.path.join(FINEST_CASE, "nsets.inp")
frd_path = os.path.join(FINEST_CASE, "analysis.frd")

nodes, elements = lib.parse_clean_mesh(clean_inp)
with open(frd_path) as f:
    frd_lines = f.readlines()
disp = lib.parse_frd_block(frd_lines, "DISP", 3)
stress = lib.parse_frd_block(frd_lines, "STRESS", 6)

tip_nodes = lib.parse_nset(nset_inp, "tip")
tip_uz_all = np.array([disp[n][2] for n in tip_nodes if n in disp])
print(f"Tip face node count: {len(tip_uz_all)}")
print(f"Full-face average Uz: {tip_uz_all.mean():.6f} mm")
print(f"Full-face std dev:    {tip_uz_all.std():.6f} mm")
print(f"Full-face min/max:    {tip_uz_all.min():.6f} / {tip_uz_all.max():.6f} mm")
print("^ Nonzero std dev across the tip face means the cross-section")
print("  is NOT staying perfectly planar/rigid under load - beam theory")
print("  assumes plane sections remain plane; 3D solid allows warping.")

# Nearest-to-centroid node
centroid = np.array([L, B/2, H/2])
tip_node_list = list(tip_nodes)
tip_coords = np.array([nodes[n] for n in tip_node_list])
dists = np.linalg.norm(tip_coords - centroid, axis=1)
nearest_idx = np.argmin(dists)
nearest_node = tip_node_list[nearest_idx]
print(f"\nNearest-to-centroid tip node: {nearest_node} "
      f"at {nodes[nearest_node]}, distance={dists[nearest_idx]:.4f}mm")
print(f"Uz at nearest-centroid node: {disp[nearest_node][2]:.6f} mm")

# Corner-only tip nodes
corner_node_ids = set()
for eid, nlist in elements.items():
    corner_node_ids.update(nlist[0:8])
tip_corners = tip_nodes & corner_node_ids
tip_uz_corners = np.array([disp[n][2] for n in tip_corners if n in disp])
print(f"\nCorner-only tip nodes: {len(tip_corners)}")
print(f"Corner-only average Uz: {tip_uz_corners.mean():.6f} mm")

print("\n--- 1b. Pointwise stress comparison: FEA top-fiber profile vs "
      "beam theory M(x)c/I ---")
print("(Sampling top surface z=4mm, mid-width y=20mm, at several x stations)")

TOL = 0.05
x_stations = [1, 3, 6, 9, 12, 20, 30, 40, 50]
print(f"\n{'x (mm)':>8} {'sigma_theory (MPa)':>20} "
      f"{'FEA top-fiber VM (MPa)':>24} {'diff %':>10}")
for xs in x_stations:
    candidates = [nid for nid, (x, y, z) in nodes.items()
                  if abs(x - xs) < TOL and abs(y - B/2) < TOL and abs(z - H) < TOL]
    if not candidates:
        # relax tolerance if no exact match
        candidates = [nid for nid, (x, y, z) in nodes.items()
                      if abs(x - xs) < 0.4 and abs(y - B/2) < 0.4 and abs(z - H) < 0.4]
    if candidates and candidates[0] in stress:
        vm = lib.von_mises(stress[candidates[0]])
        th = sigma_theory(xs)
        diff_pct = (vm - th) / th * 100
        print(f"{xs:>8} {th:>20.3f} {vm:>24.3f} {diff_pct:>9.2f}%")
    else:
        print(f"{xs:>8} {sigma_theory(xs):>20.3f} {'no node found':>24}")

print("\n" + "="*70)
print("TASK 2: Root-band extraction window sensitivity")
print("="*70)

windows = [(0, 6), (1, 7), (3, 9), (5, 11)]
for xmin, xmax in windows:
    band_nodes = [nid for nid, (x, y, z) in nodes.items() if xmin <= x <= xmax]
    band_vm = [lib.von_mises(stress[n]) for n in band_nodes if n in stress]
    if band_vm:
        max_vm = max(band_vm)
        th_at_xmin = sigma_theory(xmin)
        th_at_xmax = sigma_theory(xmax)
        print(f"\nWindow [{xmin},{xmax}] mm "
              f"(represents {L-xmax:.0f}-{L-xmin:.0f}mm from free tip, "
              f"or {xmin}-{xmax}mm from fixed root):")
        print(f"  FEA max von Mises: {max_vm:.3f} MPa")
        print(f"  Beam theory range in this window: "
              f"{th_at_xmax:.3f} (at x={xmax}) to {th_at_xmin:.3f} (at x={xmin}) MPa")

print("\n" + "="*70)
print("TASK 3: Non-monotonic stress convergence - peak location tracking")
print("="*70)

LEVELS = [
    ("2.0mm", "mesh_2.0mm"), ("1.3mm", "mesh_1.3mm"), ("0.8mm", "mesh_0.8mm"),
    ("0.5mm", "mesh_0.5mm"), ("0.3mm", "mesh_0.3mm"),
]

for label, dirname in LEVELS:
    case_dir = os.path.join(SIM_ROOT, dirname)
    c_inp = os.path.join(case_dir, "mesh_clean.inp")
    f_path = os.path.join(case_dir, "analysis.frd")

    c_nodes, _ = lib.parse_clean_mesh(c_inp)
    with open(f_path) as f:
        c_frd = f.readlines()
    c_stress = lib.parse_frd_block(c_frd, "STRESS", 6)

    band_nodes = [nid for nid, (x, y, z) in c_nodes.items() if 3 <= x <= 9]
    band_results = [(nid, lib.von_mises(c_stress[nid]))
                     for nid in band_nodes if nid in c_stress]
    band_results.sort(key=lambda t: -t[1])

    top3 = band_results[:3]
    print(f"\n{label} ({len(band_results)} nodes in band):")
    for nid, vm in top3:
        x, y, z = c_nodes[nid]
        print(f"  node {nid}: VM={vm:.3f} MPa at "
              f"x={x:.3f}, y={y:.3f}, z={z:.3f}")

print("\n" + "="*70)
print("Investigation complete.")
print("="*70)
