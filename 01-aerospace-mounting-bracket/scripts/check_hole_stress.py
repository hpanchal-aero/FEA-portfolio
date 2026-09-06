"""
Project 01 - Aerospace Mounting Bracket
Stage 2: Extract peak stress at the hole boundary from the converged
solve, compare against the Peterson/Howland Kt-based prediction.

Predicted peak (per chat log derivation):
  sigma_nominal at x=30mm (verified in prior stage) = 65.6 MPa
  Kt (gross-section, d/B=0.125) ~= 2.68
  Predicted peak ~= 175.8 MPa
"""

import sys
sys.path.insert(0, "../../scripts")
import hex_lib as lib

CASE_DIR = "."
NSET_FILE = "nsets_hole.inp"
FRD_FILE = "analysis_hole.frd"

hole_nodes = lib.parse_nset(NSET_FILE, "hole")
print(f"Parsed {len(hole_nodes)} hole boundary node IDs")

with open(FRD_FILE) as f:
    frd_lines = f.readlines()

stress = lib.parse_frd_block(frd_lines, "STRESS", 6)
disp = lib.parse_frd_block(frd_lines, "DISP", 3)

print(f"Parsed {len(stress)} stress records, {len(disp)} displacement records")

hole_vm = [(nid, lib.von_mises(stress[nid])) for nid in hole_nodes if nid in stress]
hole_vm.sort(key=lambda t: -t[1])

print(f"\nTop 10 highest von Mises stress at hole boundary:")
for nid, vm in hole_vm[:10]:
    print(f"  node {nid}: VM = {vm:.3f} MPa")

max_vm = hole_vm[0][1] if hole_vm else None
print(f"\nMax von Mises at hole boundary: {max_vm:.3f} MPa")
print(f"Predicted peak (Kt=2.68 x sigma_nom=65.6 MPa): 175.8 MPa")
if max_vm:
    diff_pct = (max_vm - 175.8) / 175.8 * 100
    print(f"Percent difference: {diff_pct:.2f}%")

# Also report tip displacement as a basic sanity check against the
# original (no-hole) verification-plate result
tip_nodes = lib.parse_nset(NSET_FILE, "tip")
tip_uz = [disp[n][2] for n in tip_nodes if n in disp]
if tip_uz:
    avg_uz = sum(tip_uz) / len(tip_uz)
    print(f"\nTip Uz average: {avg_uz:.6f} mm "
          f"(sanity check vs no-hole baseline: -1.051 mm converged)")
