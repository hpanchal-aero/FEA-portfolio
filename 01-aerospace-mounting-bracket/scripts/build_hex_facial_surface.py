"""
Project 01 - Aerospace Mounting Bracket
Build the tip facial surface for the structured hex mesh.

CalculiX C3D20 face-to-corner mapping (corner nodes 1-8):
  S1: 1,2,3,4 (face at local -z, here: bottom-ish depends on hex
      orientation - determined empirically below by coordinate check,
      not assumed)
  S2: 5,6,7,8
  S3: 1,2,6,5
  S4: 2,3,7,6
  S5: 3,4,8,7
  S6: 4,1,5,8

Rather than assume which Sn label corresponds to the tip (x=60) face
in Gmsh's corner ordering, this script checks each hex face's actual
coordinates directly - same discipline as the tet-mesh facial surface
script, cross-checked by construction here since every element's
shape is identical.
"""

import os

CLEAN_INP = "../mesh/hex_test/mesh_hex_clean.inp"
OUT_FILE = "../mesh/hex_test/tip_facial_surface_hex.inp"
TARGET_X = 60.0
TOL = 1e-6
SURFACE_NAME = "TIP_FACE"

FACE_DEFS = {
    "S1": (0,1,2,3),
    "S2": (4,5,6,7),
    "S3": (0,1,5,4),
    "S4": (1,2,6,5),
    "S5": (2,3,7,6),
    "S6": (3,0,4,7),
}

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
    lines = f.readlines()

nodes = parse_nodes(lines)
elements = parse_c3d20_elements(lines)
print(f"Parsed {len(nodes)} nodes, {len(elements)} C3D20 elements")

found = []
for eid, nlist in elements.items():
    corners = nlist[0:8]
    coords = [nodes[n] for n in corners]
    for label, idx4 in FACE_DEFS.items():
        xs = [coords[i][0] for i in idx4]
        if all(abs(x - TARGET_X) < TOL for x in xs):
            found.append((eid, label))

print(f"Found {len(found)} element faces on x={TARGET_X} plane")
expected = 50 * 5  # ny_elements x nz_elements
print(f"Expected (structured grid, 50 x 5 face grid): {expected}")

with open(OUT_FILE, "w") as f:
    f.write(f"*SURFACE,NAME={SURFACE_NAME},TYPE=ELEMENT\n")
    for eid, label in found:
        f.write(f"{eid},{label}\n")

print(f"Wrote: {OUT_FILE}")
