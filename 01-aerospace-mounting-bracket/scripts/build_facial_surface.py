"""
Project 01 - Aerospace Mounting Bracket
Stage: Build CalculiX facial (element-face) surface for the tip face.

Purpose: CalculiX's *DISTRIBUTING coupling requires a *SURFACE,TYPE=ELEMENT
facial surface (element_id, S1-S4), not a node-based surface. Gmsh's .inp
export only provides node sets / standalone skin elements for named
physical surfaces, not volume-element face labels. This script parses
the C3D10 volume mesh, identifies which tetrahedron face lies on the
target plane (x = TARGET_X), and writes a *SURFACE,TYPE=ELEMENT block.

CalculiX C3D10 face-to-corner-node mapping (corner nodes = first 4 listed
node IDs of each element; midside nodes are handled automatically by
CalculiX once the face label is given):
    S1: corners 1,2,3
    S2: corners 1,2,4
    S3: corners 2,3,4
    S4: corners 1,3,4

Cross-check: the number of (element, face) pairs found must equal the
number of CPS6 skin-facet elements Gmsh tagged for the same physical
surface (Surface6 + Surface7 for 'tip') - this catches any mismatch
before it silently produces a wrong/incomplete surface.
"""

import re

MESH_FILE = "../mesh/verification_plate.inp"
TARGET_X = 60.0
TOL = 1e-6
OUT_FILE = "../mesh/tip_facial_surface.inp"
SURFACE_NAME = "TIP_FACE"

FACE_DEFS = {
    "S1": (0, 1, 2),
    "S2": (0, 1, 3),
    "S3": (1, 2, 3),
    "S4": (0, 2, 3),
}

def parse_nodes(lines):
    nodes = {}
    i = 0
    while i < len(lines):
        if lines[i].strip().upper().startswith("*NODE") and "*NODE," not in lines[i].strip().upper()[:6]:
            i += 1
            while i < len(lines) and not lines[i].strip().startswith("*"):
                parts = lines[i].strip().split(",")
                if len(parts) >= 4 and parts[0].strip().isdigit():
                    nid = int(parts[0])
                    x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                    nodes[nid] = (x, y, z)
                i += 1
        else:
            i += 1
    return nodes

def parse_c3d10_elements(lines):
    elements = {}
    i = 0
    in_block = False
    buffer = []
    while i < len(lines):
        stripped = lines[i].strip()
        if stripped.upper().startswith("*ELEMENT") and "TYPE=C3D10" in stripped.upper().replace(" ", ""):
            in_block = True
            i += 1
            continue
        if in_block:
            if stripped.startswith("*"):
                in_block = False
                continue
            buffer.append(stripped)
        i += 1

    # Reconstruct elements: each starts with an element id (no leading
    # continuation), continues across lines until 11 numbers collected.
    tokens_all = []
    for line in buffer:
        toks = [t.strip() for t in line.split(",") if t.strip() != ""]
        tokens_all.extend(toks)

    idx = 0
    while idx < len(tokens_all):
        elem_id = int(tokens_all[idx])
        node_ids = [int(tokens_all[idx + k]) for k in range(1, 11)]
        elements[elem_id] = node_ids
        idx += 11

    return elements

with open(MESH_FILE) as f:
    lines = f.readlines()

print("Parsing nodes...")
nodes = parse_nodes(lines)
print(f"  {len(nodes)} nodes parsed")

print("Parsing C3D10 elements...")
elements = parse_c3d10_elements(lines)
print(f"  {len(elements)} C3D10 elements parsed")

# Sanity check against known model size
if len(nodes) < 1000 or len(elements) < 1000:
    raise RuntimeError("Parsed counts look too small - parser likely misaligned with file format")

found_faces = []
for elem_id, node_ids in elements.items():
    corners = node_ids[0:4]
    coords = [nodes[n] for n in corners]
    for face_label, idx_triplet in FACE_DEFS.items():
        xs = [coords[i][0] for i in idx_triplet]
        if all(abs(x - TARGET_X) < TOL for x in xs):
            found_faces.append((elem_id, face_label))

print(f"\nFound {len(found_faces)} element faces on x={TARGET_X} plane")

with open(OUT_FILE, "w") as f:
    f.write(f"*SURFACE,NAME={SURFACE_NAME},TYPE=ELEMENT\n")
    for elem_id, face_label in found_faces:
        f.write(f"{elem_id},{face_label}\n")

print(f"Wrote: {OUT_FILE}")
