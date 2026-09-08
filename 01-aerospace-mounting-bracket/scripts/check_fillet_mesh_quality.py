"""
Stage 3 (fillet) — mesh quality check.

Adapted from the Stage 2 hole-mesh quality logic (signed-volume /
Jacobian-sign check + edge-length aspect ratio), applied here to the
Stage 3 fillet baseline mesh. Numpy-vectorized for speed, per the
approach already validated on the much larger Stage 2 mesh.

Checks:
  1. Signed volume of every C3D10 tet (via its 4 corner nodes) -
     flags inverted/degenerate elements (nonpositive Jacobian).
  2. Edge-length aspect ratio (longest/shortest corner-edge per
     element) - flags elements with pathological aspect ratios,
     which matters here specifically because the mesh transitions
     from 0.15mm to 3mm over only 1.5mm at the fillet.

Does not modify the mesh file. Read-only diagnostic.
"""

import numpy as np
import re

MESH_FILE = "01-aerospace-mounting-bracket/mesh/fillet_study/mesh_fillet_raw.inp"

def read_abaqus_mesh(path):
    nodes = {}
    elements = {}  # elem_id -> list of node ids (10, for C3D10)
    mode = None
    elem_type = None

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("**"):
                continue
            if line.upper().startswith("*NODE"):
                mode = "NODE"
                continue
            if line.upper().startswith("*ELEMENT"):
                mode = "ELEMENT"
                m = re.search(r"TYPE\s*=\s*(\w+)", line, re.IGNORECASE)
                elem_type = m.group(1).upper() if m else None
                continue
            if line.startswith("*"):
                mode = None
                continue

            if mode == "NODE":
                parts = line.split(",")
                nid = int(parts[0])
                x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                nodes[nid] = (x, y, z)
            elif mode == "ELEMENT" and elem_type == "C3D10":
                parts = [p.strip() for p in line.split(",")]
                eid = int(parts[0])
                conn = [int(p) for p in parts[1:11]]
                elements[eid] = conn

    return nodes, elements

print(f"Reading mesh: {MESH_FILE}")
nodes, elements = read_abaqus_mesh(MESH_FILE)
print(f"Nodes read: {len(nodes)}")
print(f"C3D10 elements read: {len(elements)}")

if len(elements) == 0:
    print("FATAL: no C3D10 elements found - check file format / TYPE string.")
    raise SystemExit(1)

elem_ids = np.array(list(elements.keys()))
conn_array = np.array([elements[e] for e in elem_ids])  # (N, 10)

max_node_id = max(nodes.keys())
node_coords = np.zeros((max_node_id + 1, 3))
for nid, coord in nodes.items():
    node_coords[nid] = coord

# Corner nodes of a C3D10 are the first 4 entries of the connectivity
corners = conn_array[:, 0:4]  # (N, 4)
p0 = node_coords[corners[:, 0]]
p1 = node_coords[corners[:, 1]]
p2 = node_coords[corners[:, 2]]
p3 = node_coords[corners[:, 3]]

# ---- Signed volume check ----
v1 = p1 - p0
v2 = p2 - p0
v3 = p3 - p0
cross_v2v3 = np.cross(v2, v3)
signed_vol = np.einsum('ij,ij->i', v1, cross_v2v3) / 6.0

n_inverted = np.sum(signed_vol <= 0)
n_positive = np.sum(signed_vol > 0)

print("\n=== SIGNED VOLUME / JACOBIAN CHECK ===")
print(f"Positive-volume elements: {n_positive}")
print(f"Nonpositive-volume (inverted/degenerate) elements: {n_inverted}")
if n_inverted > 0:
    bad_ids = elem_ids[signed_vol <= 0]
    print(f"Inverted element IDs (first 20 shown): {bad_ids[:20]}")
    print(f"Corresponding signed volumes (first 20): {signed_vol[signed_vol <= 0][:20]}")
else:
    print("RESULT: zero inverted elements.")

# ---- Edge-length aspect ratio check (6 corner-to-corner edges per tet) ----
edge_pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
edge_lengths = np.zeros((len(elem_ids), 6))
for i, (a, b) in enumerate(edge_pairs):
    pa = node_coords[corners[:, a]]
    pb = node_coords[corners[:, b]]
    edge_lengths[:, i] = np.linalg.norm(pa - pb, axis=1)

min_edge = edge_lengths.min(axis=1)
max_edge = edge_lengths.max(axis=1)
aspect_ratio = max_edge / min_edge

print("\n=== EDGE-LENGTH ASPECT RATIO CHECK ===")
print(f"Min aspect ratio: {aspect_ratio.min():.2f}")
print(f"Max aspect ratio: {aspect_ratio.max():.2f}")
print(f"Mean aspect ratio: {aspect_ratio.mean():.2f}")
print(f"Median aspect ratio: {np.median(aspect_ratio):.2f}")

# Flag thresholds - CalculiX/general FEA guidance: AR > 10 is a
# concern, AR > 20 is a strong concern, particularly at a size-
# transition zone like this one.
thresholds = [5, 10, 20, 50]
for t in thresholds:
    n_over = np.sum(aspect_ratio > t)
    pct = 100.0 * n_over / len(aspect_ratio)
    print(f"Elements with aspect ratio > {t}: {n_over} ({pct:.3f}%)")

worst_idx = np.argsort(aspect_ratio)[::-1][:10]
print("\nWorst 10 elements by aspect ratio:")
for idx in worst_idx:
    eid = elem_ids[idx]
    ar = aspect_ratio[idx]
    centroid = node_coords[corners[idx]].mean(axis=0)
    print(f"  Element {eid}: AR={ar:.2f}, centroid=({centroid[0]:.3f}, {centroid[1]:.3f}, {centroid[2]:.3f})")

# ---- Element size sanity vs. transition zone location ----
# Report min/max edge length overall, to confirm the mesh actually
# reaches both the intended FILLET_SIZE (~0.15mm) and FAR_SIZE (~3mm).
print("\n=== ELEMENT SIZE RANGE CHECK ===")
print(f"Smallest edge in mesh: {edge_lengths.min():.4f} mm (expect near FILLET_SIZE=0.15mm)")
print(f"Largest edge in mesh:  {edge_lengths.max():.4f} mm (expect near FAR_SIZE=3.0mm)")

print("\n=== SUMMARY ===")
print(f"Total elements: {len(elem_ids)}")
print(f"Inverted elements: {n_inverted}  {'-- FATAL, do not proceed to solve' if n_inverted > 0 else '-- OK'}")
print(f"Elements with AR > 20: {np.sum(aspect_ratio > 20)}  {'-- flag for review' if np.sum(aspect_ratio > 20) > 0 else '-- OK'}")
