"""
Project 01 - Aerospace Mounting Bracket
Diagnostic: test the biaxial-stress (Poisson/anticlastic-restraint)
hypothesis for the FEA-vs-beam-theory stress discrepancy, by
extracting raw stress tensor components (not just von Mises) at
several stations along the span, top surface, mid-width.

Beam theory assumes pure uniaxial sigma_xx = Mc/I. If the true FEA
stress state is biaxial (sigma_yy present due to width restraint),
we should see sigma_yy/sigma_xx approaching Poisson's ratio (0.33)
away from boundary effects.
"""

import os
import sys
sys.path.insert(0, ".")
import hex_lib as lib

CASE_DIR = "../simulation/convergence_hex/mesh_0.3mm"
L, B, H = lib.L, lib.B, lib.H
NU = 0.33

clean_inp = os.path.join(CASE_DIR, "mesh_clean.inp")
frd_path = os.path.join(CASE_DIR, "analysis.frd")

nodes, elements = lib.parse_clean_mesh(clean_inp)
with open(frd_path) as f:
    frd_lines = f.readlines()
stress = lib.parse_frd_block(frd_lines, "STRESS", 6)

TOL = 0.05
x_stations = [3, 6, 9, 12, 20, 30, 40]

print(f"{'x':>6} {'sigma_xx':>10} {'sigma_yy':>10} {'sigma_zz':>10} "
      f"{'sigma_xy':>10} {'ratio_yy/xx':>12} {'nu target':>10}")
for xs in x_stations:
    candidates = [nid for nid, (x, y, z) in nodes.items()
                  if abs(x - xs) < TOL and abs(y - B/2) < TOL and abs(z - H) < TOL]
    if not candidates:
        candidates = [nid for nid, (x, y, z) in nodes.items()
                      if abs(x - xs) < 0.4 and abs(y - B/2) < 0.4 and abs(z - H) < 0.4]
    if candidates and candidates[0] in stress:
        sxx, syy, szz, sxy, syz, szx = stress[candidates[0]]
        ratio = syy / sxx if abs(sxx) > 1e-6 else float('nan')
        print(f"{xs:>6} {sxx:>10.3f} {syy:>10.3f} {szz:>10.3f} "
              f"{sxy:>10.3f} {ratio:>12.4f} {NU:>10.2f}")
    else:
        print(f"{xs:>6}  no node found")
