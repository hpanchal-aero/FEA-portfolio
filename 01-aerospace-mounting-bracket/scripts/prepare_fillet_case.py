"""
Project 01 - Aerospace Mounting Bracket
Stage 3 (fillet study): strip skin elements and build the tip facial
surface for the shoulder-fillet case, reusing the exact validated
logic from Stage 2's prepare_hole_case.py (skin-element stripping,
coordinate-matched facial-surface construction).

NSETs (root, fillet_top, fillet_bottom) are already written by
build_fillet_case.py in the same Gmsh session as the mesh - no reopen
needed here, unlike Stage 2's .geo-reopen approach.

Run from mesh/fillet_study/ - operates on mesh_fillet_raw.inp in that
directory.
"""

import os

RAW_INP = "mesh_fillet_raw.inp"
CLEAN_INP = "mesh_fillet_clean.inp"
FACIAL_OUT = "tip_facial_surface_fillet.inp"

L = 60.0

# --- Step 1: strip CPS6 skin elements, keep only *NODE + C3D10 volume ---
with open(RAW_INP) as f:
    lines = f.readlines()

out_lines = []
skip = False
kept, dropped = 0, 0
for line in lines:
    upper = line.strip().upper().replace(" ", "")
    if upper.startswith("*ELEMENT,TYPE=CPS6"):
        skip = True
        dropped += 1
        continue
    if upper.startswith("*ELEMENT,TYPE=C3D10"):
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
print(f"Kept C3D10 blocks: {kept}, Dropped CPS6 blocks: {dropped}")
print(f"Wrote: {CLEAN_INP}")

# --- Step 2: build the tip facial surface (x=L plane) for DISTRIBUTING
#     coupling - identical coordinate-search logic to Stage 2 ---
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

def parse_c3d10_elements(lines):
    elements = {}
    in_block = False
    buffer = []
    for line in lines:
        s = line.strip()
        if s.upper().startswith("*ELEMENT") and "TYPE=C3D10" in s.upper().replace(" ", ""):
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
        node_ids = [int(tokens[idx + k]) for k in range(1, 11)]
        elements[elem_id] = node_ids
        idx += 11
    return elements

with open(CLEAN_INP) as f:
    clean_lines = f.readlines()

nodes = parse_nodes(clean_lines)
elements = parse_c3d10_elements(clean_lines)
print(f"\nParsed {len(nodes)} nodes, {len(elements)} C3D10 elements")

FACE_DEFS = {
    "S1": (0, 1, 2), "S2": (0, 1, 3),
    "S3": (1, 2, 3), "S4": (0, 2, 3),
}
TOL = 1e-6
found = []
for elem_id, node_ids in elements.items():
    corners = node_ids[0:4]
    coords = [nodes[n] for n in corners]
    for label, idx_triplet in FACE_DEFS.items():
        xs = [coords[i][0] for i in idx_triplet]
        if all(abs(x - L) < TOL for x in xs):
            found.append((elem_id, label))

print(f"Found {len(found)} element faces on x={L} plane")

with open(FACIAL_OUT, "w") as f:
    f.write("*SURFACE,NAME=TIP_FACE,TYPE=ELEMENT\n")
    for elem_id, label in found:
        f.write(f"{elem_id},{label}\n")
print(f"Wrote: {FACIAL_OUT}")

max_node_id = max(nodes.keys())
print(f"\nMax real node ID: {max_node_id} (use reference node ID "
      f"{max_node_id + 1000} to avoid collision)")
