"""
Project 01 - Aerospace Mounting Bracket
Stage 4a (right-angle frame): strip skin elements and build the tip
facial surface, reusing the exact validated logic from Stage 3's
prepare_fillet_case.py.

Tip face is at z=60 (Flange B's free end), not x=60 as in prior
stages - the facial-surface search below matches on z, not x.

Run from mesh/frame_study/ - operates on mesh_frame_raw.inp.
"""

RAW_INP = "mesh_frame_raw.inp"
CLEAN_INP = "mesh_frame_clean.inp"
FACIAL_OUT = "tip_facial_surface_frame.inp"

TIP_Z = 60.0

# --- Step 1: strip CPS6 skin elements ---
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

# --- Step 2: build the tip facial surface (z=TIP_Z plane) ---
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
        zs = [coords[i][2] for i in idx_triplet]
        if all(abs(z - TIP_Z) < TOL for z in zs):
            found.append((elem_id, label))

print(f"Found {len(found)} element faces on z={TIP_Z} plane")

with open(FACIAL_OUT, "w") as f:
    f.write("*SURFACE,NAME=TIP_FACE,TYPE=ELEMENT\n")
    for elem_id, label in found:
        f.write(f"{elem_id},{label}\n")
print(f"Wrote: {FACIAL_OUT}")

max_node_id = max(nodes.keys())
print(f"\nMax real node ID: {max_node_id} (use reference node ID "
      f"{max_node_id + 1000} to avoid collision)")
