"""
Project 01 - Aerospace Mounting Bracket
Stage 2 (hole study): mesh quality check for the C3D10 tet mesh
around the hole, before committing to an expensive solve.

Checks (vectorized with numpy for speed given ~824k elements):
  - element count sanity
  - signed volume of every element (negative = inverted/degenerate)
  - edge length distribution (confirms refinement gradient is as
    intended: ~0.15mm near hole, ~0.5mm far-field)
"""

import numpy as np

MESH_FILE = "mesh_hole_raw.inp"

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

print("Parsing mesh file...")
with open(MESH_FILE) as f:
    lines = f.readlines()

nodes = parse_nodes(lines)
elements = parse_c3d10_elements(lines)
print(f"Parsed {len(nodes)} nodes, {len(elements)} C3D10 elements")

elem_ids = list(elements.keys())
corner_arrays = np.array([[nodes[n] for n in elements[eid][0:4]]
                           for eid in elem_ids])  # shape (N,4,3)

p0 = corner_arrays[:, 0, :]
p1 = corner_arrays[:, 1, :]
p2 = corner_arrays[:, 2, :]
p3 = corner_arrays[:, 3, :]

v1 = p1 - p0
v2 = p2 - p0
v3 = p3 - p0
signed_vol = np.sum(v1 * np.cross(v2, v3), axis=1) / 6.0

neg_count = np.sum(signed_vol < -1e-9)
zero_count = np.sum(np.abs(signed_vol) < 1e-9)

print(f"\nSigned volume check ({len(elements)} elements):")
print(f"  Negative-volume (inverted): {neg_count}")
print(f"  Near-zero (degenerate): {zero_count}")
print(f"  Min volume: {signed_vol.min():.6e} mm^3")
print(f"  Max volume: {signed_vol.max():.6e} mm^3")
print(f"  Mean volume: {signed_vol.mean():.6e} mm^3")

# Edge length distribution (all 6 corner-to-corner edges per element)
edge_pairs = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
all_edges = []
for i, j in edge_pairs:
    d = np.linalg.norm(corner_arrays[:, i, :] - corner_arrays[:, j, :], axis=1)
    all_edges.append(d)
all_edges = np.concatenate(all_edges)

print(f"\nEdge length distribution (all corner edges, {len(all_edges)} total):")
print(f"  Min: {all_edges.min():.4f} mm")
print(f"  Max: {all_edges.max():.4f} mm")
print(f"  5th percentile: {np.percentile(all_edges, 5):.4f} mm")
print(f"  95th percentile: {np.percentile(all_edges, 95):.4f} mm")
print(f"  Target near hole: ~0.15 mm, far-field: ~0.5 mm")

print("\n" + "="*50)
if neg_count == 0 and zero_count == 0:
    print("MESH QUALITY: PASSED (no inverted/degenerate elements)")
else:
    print("MESH QUALITY: FAILED - see counts above")
print("="*50)
