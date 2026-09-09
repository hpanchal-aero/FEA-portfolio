"""
Stage 4b - targeted convergence check at the gusset toe (x~38, z~2).

For each mesh level, extracts:
  1. The peak |sigma_xx| among nodes very close to the toe line
     (x in [37.5, 38.5], z in [1.5, 2.5]) - tracking both value and
     exact (x,y,z) location.
  2. sigma_xx sampled at increasing distance from the toe along the
     flange (x = 38, 36, 34, 32, 30, 25, 20mm, at z=2, y~20) - to see
     whether stress a finite distance away is stable even if the toe
     value itself is not.

Run once per level, pointed at that level's own .frd/nset files via
the LEVEL_DIR argument.
"""

import sys
import numpy as np

LEVEL_DIR = sys.argv[1] if len(sys.argv) > 1 else "."
FRD_FILE = f"{LEVEL_DIR}/analysis_frame_gusset.frd"


def parse_frd(path):
    coords, stress = {}, {}
    with open(path) as f:
        lines = f.readlines()
    i, n = 0, len(lines)
    while i < n:
        line = lines[i]
        if line.startswith("    2C"):
            i += 1
            while i < n and lines[i].startswith(" -1"):
                rec = lines[i].rstrip("\n")
                nid = int(rec[3:13])
                coords[nid] = (float(rec[13:25]), float(rec[25:37]), float(rec[37:49]))
                i += 1
            continue
        if "STRESS" in line and line.startswith(" -4"):
            i += 1
            while i < n and lines[i].startswith(" -5"):
                i += 1
            while i < n and lines[i].startswith(" -1"):
                rec = lines[i].rstrip("\n")
                nid = int(rec[3:13])
                vals = [float(rec[13 + 12*k:25 + 12*k]) for k in range(6)]
                stress[nid] = tuple(vals)
                i += 1
            continue
        i += 1
    return coords, stress


print(f"Reading {FRD_FILE}")
coords, stress = parse_frd(FRD_FILE)
print(f"Parsed {len(coords)} coords, {len(stress)} stress results")

# ---- 1. Peak stress right at the toe ----
toe_candidates = []
for nid, (x, y, z) in coords.items():
    if 37.5 <= x <= 38.5 and 1.5 <= z <= 2.5:
        if nid in stress:
            toe_candidates.append((nid, x, y, z, stress[nid][0]))

toe_candidates.sort(key=lambda t: abs(t[4]), reverse=True)
print(f"\n=== PEAK STRESS AT TOE (x in [37.5,38.5], z in [1.5,2.5]) ===")
print(f"Candidates found: {len(toe_candidates)}")
print("Top 5 by |sigma_xx|:")
for nid, x, y, z, sxx in toe_candidates[:5]:
    print(f"  node {nid} at ({x:.3f},{y:.3f},{z:.3f}): sigma_xx={sxx:.4f} MPa")

if toe_candidates:
    peak = toe_candidates[0]
    print(f"\nPEAK: {peak[4]:.4f} MPa at ({peak[1]:.3f},{peak[2]:.3f},{peak[3]:.3f})")

# ---- 2. Stress trend away from the toe, along the flange (z=2, y~20) ----
print(f"\n=== STRESS TREND AWAY FROM TOE (z~2, y~20) ===")
print(f"{'x target':>10} {'sigma_xx (MPa)':>16} {'n':>4}")
for x_target in [38, 36, 34, 32, 30, 25, 20]:
    matches = []
    for nid, (x, y, z) in coords.items():
        if abs(x - x_target) < 0.5 and abs(z - 2.0) < 0.3 and abs(y - 20.0) < 1.5:
            if nid in stress:
                matches.append(stress[nid][0])
    if matches:
        print(f"{x_target:>10} {np.mean(matches):>16.4f} {len(matches):>4}")
    else:
        print(f"{x_target:>10} {'NO MATCH':>16}")
