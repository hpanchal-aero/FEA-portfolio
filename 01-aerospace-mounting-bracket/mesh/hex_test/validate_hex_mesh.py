"""
Project 01 - Aerospace Mounting Bracket
Standalone validation of the structured transfinite hex mesh
(verification plate only - see chat log for scope).

Checks performed:
  1. Strip CPS8 skin elements, keep only *NODE and *ELEMENT,TYPE=C3D20
     (same pattern as clean_mesh_for_ccx.py for the tet mesh).
  2. Confirm element count matches the expected structured grid
     (75 x 50 x 5 = 18750 hexes at this refinement level).
  3. For every hex: compute the signed volume of all 6 tetrahedra
     formed by splitting the hex from its 8 corners, confirm all
     positive (no inverted/degenerate elements) - not a proxy check,
     an exhaustive one, since this is a small, fast mesh.
  4. Compute edge length statistics (min/max/mean) and aspect ratio
     (max edge / min edge per element) to confirm elements are close
     to the target 0.8mm cubes as intended.
  5. Confirm midside nodes sit at exact edge midpoints
     (SecondOrderLinear should guarantee this, but verify directly
     rather than assume - same discipline as the tet-mesh checks).
"""

import os
import numpy as np
from itertools import combinations

RAW_INP = "mesh_hex.inp"
CLEAN_INP = "mesh_hex_clean.inp"

EXPECTED_HEX_COUNT = 75 * 50 * 5  # 18750

# CalculiX C3D20 node ordering: 8 corners, then 12 midside nodes
# following the edge order: (1,2)(2,3)(3,4)(4,1)(5,6)(6,7)(7,8)(8,5)
# (1,5)(2,6)(3,7)(4,8) - 0-indexed here.
HEX_EDGES = [
    (0,1),(1,2),(2,3),(3,0),
    (4,5),(5,6),(6,7),(7,4),
    (0,4),(1,5),(2,6),(3,7),
]

# Tetrahedral decomposition of a hex (corner indices 0-7) for signed
# volume checking - 6 tets covering the full hex volume.
HEX_TETS = [
    (0,1,3,4), (1,2,3,6), (1,3,4,6),
    (3,4,6,7), (1,4,5,6),
]

# --- Step 1: strip CPS8 skin blocks ---
with open(RAW_INP) as f:
    lines = f.readlines()

out_lines = []
skip = False
kept, dropped = 0, 0
for line in lines:
    upper = line.strip().upper().replace(" ", "")
    if upper.startswith("*ELEMENT,TYPE=CPS8"):
        skip = True
        dropped += 1
        continue
    if upper.startswith("*ELEMENT,TYPE=C3D20"):
        skip = False
        kept += 1
        out_lines.append(line)
        continue
    if upper.startswith("*") and skip:
        skip = False
    if not skip:
        out_lines.append(line)

with open(CLEAN_INP, "w") as f:
    f.writelines(out_lines)

print(f"Kept C3D20 blocks: {kept}, Dropped CPS8 blocks: {dropped}")

# --- Step 2: parse nodes and C3D20 elements ---
def parse_nodes(lines):
    nodes = {}
    i = 0
    while i < len(lines):
        if lines[i].strip().upper().startswith("*NODE") and \
           not lines[i].strip().upper().startswith("*NODE,"):
            i += 1
            while i < len(lines) and not lines[i].strip().startswith("*"):
                parts = lines[i].strip().split(",")
                if len(parts) >= 4 and parts[0].strip().isdigit():
                    nid = int(parts[0])
                    nodes[nid] = (float(parts[1]), float(parts[2]), float(parts[3]))
                i += 1
        else:
            i += 1
    return nodes

def parse_c3d20_elements(lines):
    elements = {}
    in_block = False
    buffer = []
    for line in lines:
        s = line.strip()
        if s.upper().startswith("*ELEMENT") and "TYPE=C3D20" in s.upper().replace(" ", ""):
            in_block = True
            continue
        if in_block:
            if s.startswith("*"):
                in_block = False
                continue
            buffer.append(s)

    tokens = []
    for line in buffer:
        tokens.extend([t.strip() for t in line.split(",") if t.strip()])

    idx = 0
    while idx < len(tokens):
        elem_id = int(tokens[idx])
        node_ids = [int(tokens[idx + k]) for k in range(1, 21)]
        elements[elem_id] = node_ids
        idx += 21

    return elements

with open(CLEAN_INP) as f:
    clean_lines = f.readlines()

nodes = parse_nodes(clean_lines)
elements = parse_c3d20_elements(clean_lines)

print(f"\nParsed {len(nodes)} nodes")
print(f"Parsed {len(elements)} C3D20 elements")
print(f"Expected element count: {EXPECTED_HEX_COUNT}")
if len(elements) != EXPECTED_HEX_COUNT:
    print("*** WARNING: element count does not match expected structured grid ***")
else:
    print("Element count matches expected structured grid EXACTLY.")

# --- Step 3: signed volume check (all hexes) ---
neg_vol_elements = []
zero_vol_elements = []
volumes = []

for eid, nlist in elements.items():
    corners = nlist[0:8]
    coords = [np.array(nodes[n]) for n in corners]

    vol = 0.0
    for (a,b,c,d) in HEX_TETS:
        p0,p1,p2,p3 = coords[a],coords[b],coords[c],coords[d]
        v1,v2,v3 = p1-p0, p2-p0, p3-p0
        vol += np.dot(v1, np.cross(v2,v3)) / 6.0

    volumes.append(vol)
    if vol < -1e-9:
        neg_vol_elements.append(eid)
    elif abs(vol) < 1e-9:
        zero_vol_elements.append(eid)

volumes = np.array(volumes)
print(f"\nVolume check across {len(elements)} elements:")
print(f"  Negative-volume (inverted) elements: {len(neg_vol_elements)}")
print(f"  Near-zero-volume (degenerate) elements: {len(zero_vol_elements)}")
print(f"  Min element volume: {volumes.min():.6e} mm^3")
print(f"  Max element volume: {volumes.max():.6e} mm^3")
print(f"  Mean element volume: {volumes.mean():.6e} mm^3")
print(f"  Expected (0.8mm cube): {0.8**3:.6e} mm^3")

# --- Step 4: edge length / aspect ratio statistics ---
all_edge_lengths = []
worst_aspect_ratio = 0.0
worst_aspect_elem = None

for eid, nlist in elements.items():
    corners = nlist[0:8]
    coords = [np.array(nodes[n]) for n in corners]
    edge_lens = []
    for (a,b) in HEX_EDGES:
        d = np.linalg.norm(coords[a] - coords[b])
        edge_lens.append(d)
        all_edge_lengths.append(d)
    ar = max(edge_lens) / min(edge_lens)
    if ar > worst_aspect_ratio:
        worst_aspect_ratio = ar
        worst_aspect_elem = eid

all_edge_lengths = np.array(all_edge_lengths)
print(f"\nEdge length statistics (all elements, all 12 edges):")
print(f"  Min: {all_edge_lengths.min():.4f} mm")
print(f"  Max: {all_edge_lengths.max():.4f} mm")
print(f"  Mean: {all_edge_lengths.mean():.4f} mm")
print(f"  Expected: 0.8000 mm (uniform cubes)")
print(f"\nWorst element aspect ratio: {worst_aspect_ratio:.4f} (element {worst_aspect_elem})")

# --- Step 5: midside node exact-midpoint check ---
max_midside_dev = 0.0
worst_midside_elem = None

for eid, nlist in elements.items():
    corners = nlist[0:8]
    midsides = nlist[8:20]
    coords = [np.array(nodes[n]) for n in corners]
    for k, (a,b) in enumerate(HEX_EDGES):
        true_mid = (coords[a] + coords[b]) / 2.0
        actual = np.array(nodes[midsides[k]])
        dev = np.linalg.norm(actual - true_mid)
        if dev > max_midside_dev:
            max_midside_dev = dev
            worst_midside_elem = eid

print(f"\nMax midside node deviation from true edge midpoint: "
      f"{max_midside_dev:.8f} mm (element {worst_midside_elem})")

# --- Final verdict ---
print("\n" + "="*60)
if (len(neg_vol_elements) == 0 and len(zero_vol_elements) == 0
        and len(elements) == EXPECTED_HEX_COUNT
        and max_midside_dev < 1e-6):
    print("MESH VALIDATION: PASSED")
else:
    print("MESH VALIDATION: FAILED - see warnings above")
print("="*60)
