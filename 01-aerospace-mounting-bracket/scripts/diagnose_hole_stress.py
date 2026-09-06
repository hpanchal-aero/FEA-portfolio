"""
Project 01 - Aerospace Mounting Bracket
Stage 2: Diagnose the sigma_xx/VM discrepancy at the hole boundary -
check z-location of peak stress, local stress tensor state, and
whether comparison against top-fiber nominal stress is appropriate.
"""

import sys
sys.path.insert(0, "../../scripts")
import hex_lib as lib

NSET_FILE = "nsets_hole.inp"
FRD_FILE = "analysis_hole.frd"
CLEAN_INP = "mesh_hole_clean.inp"

L, B, H = lib.L, lib.B, lib.H
F = 235.44
I_SECTION = B * H**3 / 12.0
C = H / 2.0

def sigma_theory(x):
    M = F * (L - x)
    return M * C / I_SECTION

nodes, elements = lib.parse_clean_mesh(CLEAN_INP)
hole_nodes = lib.parse_nset(NSET_FILE, "hole")

with open(FRD_FILE) as f:
    frd_lines = f.readlines()
stress = lib.parse_frd_block(frd_lines, "STRESS", 6)

hole_vm = [(nid, lib.von_mises(stress[nid])) for nid in hole_nodes if nid in stress]
hole_vm.sort(key=lambda t: -t[1])

print("Top 10 stress locations - full detail:")
print(f"{'node':>8} {'x':>7} {'y':>7} {'z':>7} {'VM':>9} "
      f"{'sxx':>9} {'syy':>9} {'szz':>9}")
for nid, vm in hole_vm[:10]:
    x, y, z = nodes[nid]
    sxx, syy, szz, sxy, syz, szx = stress[nid]
    print(f"{nid:>8} {x:>7.3f} {y:>7.3f} {z:>7.3f} {vm:>9.3f} "
          f"{sxx:>9.3f} {syy:>9.3f} {szz:>9.3f}")

# Distribution of max VM as a function of z (bin hole nodes by z)
print("\nMax VM stress at hole boundary, binned by z (thickness position):")
z_bins = {}
for nid in hole_nodes:
    if nid not in stress:
        continue
    x, y, z = nodes[nid]
    z_rounded = round(z, 1)
    vm = lib.von_mises(stress[nid])
    if z_rounded not in z_bins or vm > z_bins[z_rounded][1]:
        z_bins[z_rounded] = (nid, vm)

for z in sorted(z_bins.keys()):
    nid, vm = z_bins[z]
    sigma_nom_at_z = sigma_theory(30) * (2 * z / H - 1)  # linear through thickness
    print(f"  z={z:.2f}: max VM={vm:.3f} MPa "
          f"(local nominal beam bending stress at this z: "
          f"{sigma_nom_at_z:.3f} MPa)")

print(f"\nFor reference: top fiber (z=4) nominal = {sigma_theory(30):.3f} MPa, "
      f"mid-plane (z=2) nominal = 0 MPa (neutral axis)")
