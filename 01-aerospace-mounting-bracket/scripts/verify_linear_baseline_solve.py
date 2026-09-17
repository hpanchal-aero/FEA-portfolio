"""
Stage 5 -- verify the TEMPORARY linear (C3D4) baseline solve.

Checks:
  1. Global equilibrium: sum of reaction forces at root should equal
     the applied load (235.44N in -x), within numerical tolerance.
  2. Tip deflection: sanity-check magnitude (not a target -- this is
     C3D4, expected to be stiffer than a converged quadratic result,
     per the known locking limitation already flagged).
  3. Peak von Mises stress location: sanity check that it lands near
     an expected feature (root, toe fillet, or R3 fillet), not in an
     unremarkable background region (which would suggest a modeling
     error rather than genuine stress concentration).
"""

import re
import numpy as np

DAT_OR_FRD = "simulation/stage5_level1_LINEAR_TEMP_solve.dat"
FRD_FILE = "simulation/stage5_level1_LINEAR_TEMP_solve.frd"

def parse_frd_reactions_and_disp(frd_path):
    """Parse .frd for nodal displacements (U) and reaction forces (RF)."""
    disp = {}
    rf = {}
    with open(frd_path) as f:
        lines = f.readlines()

    mode = None
    for line in lines:
        if "DISP" in line and line.startswith(" -4"):
            mode = "disp"
            continue
        if "FORC" in line and line.startswith(" -4"):
            mode = "force"
            continue
        if line.startswith(" -3"):
            mode = None
            continue
        if mode == "disp" and line.startswith(" -1"):
            parts = line.split()
            node = int(parts[1])
            vals = [float(v) for v in parts[2:5]]
            disp[node] = vals
        if mode == "force" and line.startswith(" -1"):
            parts = line.split()
            node = int(parts[1])
            vals = [float(v) for v in parts[2:5]]
            rf[node] = vals

    return disp, rf

disp, rf = parse_frd_reactions_and_disp(FRD_FILE)
print(f"Parsed {len(disp)} displacement records, {len(rf)} reaction force records")

# --- Equilibrium check: sum RF at root nodes should equal -FORCE_X ---
root_rf_x = [v[0] for n, v in rf.items() if v is not None]
if root_rf_x:
    total_rf_x = sum(root_rf_x)
    print(f"\nSum of Fx reactions (all nodes with RF recorded): {total_rf_x:.4f} N")
    print(f"Applied load: -235.44 N (in -x)")
    print(f"Expected reaction sum: +235.44 N (root pushes back)")
    err = abs(total_rf_x - 235.44) / 235.44 * 100
    print(f"Equilibrium error: {err:.4f}%")
else:
    print("WARNING: no reaction force data parsed -- check .frd format/keyword.")

# --- Tip deflection sanity check ---
tip_disps = []
for n, v in disp.items():
    tip_disps.append(v)
if tip_disps:
    tip_disps = np.array(tip_disps)
    max_ux = tip_disps[:, 0].min()  # most negative x-displacement
    print(f"\nMax |Ux| displacement found in model: {max_ux:.4f} mm")
    print(f"(For reference -- Stage 4a/4b tip deflections were -4.22mm / -2.46mm;")
    print(f" this is a DIFFERENT, stiffer geometry with C3D4 locking, so no direct")
    print(f" comparison target -- sanity check only: should be same order of magnitude,")
    print(f" i.e. roughly 0.1-5mm, not near-zero or absurdly large)")

print("\nNOTE: for precise per-node values (specific tip node deflection, peak")
print("stress location), a more targeted extraction using tip/root node ID")
print("lists from the .inp would be needed -- this is a first-pass sanity check.")
