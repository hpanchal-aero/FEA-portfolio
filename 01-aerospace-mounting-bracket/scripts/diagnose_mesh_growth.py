"""
Stage 3 (fillet study) - diagnostic only. Does NOT generate a mesh or
run CalculiX. Buckets nodes/elements of an already-generated mesh by
approximate distance to the nearest fillet surface, to determine
whether disproportionate node growth (0.075mm -> 0.06mm) comes from
the immediate fillet-size band, the fixed 1.5mm transition band, or
elsewhere.

Distance approximation: the fillet surface is a true circular arc
(constant cross-section, extruded along y), locked geometry:
  X_TRANS = 30.0, HALF_TIP = 2.0, FILLET_R = 0.95
  Top arc center:    (X_TRANS + FILLET_R, +HALF_TIP + FILLET_R)
  Bottom arc center: (X_TRANS + FILLET_R, -HALF_TIP - FILLET_R)
Distance from a point (x, z) with z>0 to the top arc (and symmetric
for z<0) is |sqrt((x-xc)^2+(z-zc)^2) - FILLET_R|. This is exact for
points whose nearest point lies within the arc's actual angular
span, and a slight over/underestimate very close to the arc's two
tangent endpoints - noted as an approximation, not exact geometry
query.
"""

import re
import numpy as np

X_TRANS = 30.0
HALF_TIP = 2.0
FILLET_R = 0.95
TRANSITION_WIDTH = 1.5

TOP_CENTER = (X_TRANS + FILLET_R, HALF_TIP + FILLET_R)
BOTTOM_CENTER = (X_TRANS + FILLET_R, -HALF_TIP - FILLET_R)


def read_mesh(path):
    nodes = {}
    elements = {}
    mode = None
    elem_type = None
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("**"):
                continue
            if s.upper().startswith("*NODE"):
                mode = "NODE"
                continue
            if s.upper().startswith("*ELEMENT"):
                mode = "ELEMENT"
                m = re.search(r"TYPE\s*=\s*(\w+)", s, re.IGNORECASE)
                elem_type = m.group(1).upper() if m else None
                continue
            if s.startswith("*"):
                mode = None
                continue
            if mode == "NODE":
                parts = s.split(",")
                nid = int(parts[0])
                nodes[nid] = (float(parts[1]), float(parts[2]), float(parts[3]))
            elif mode == "ELEMENT" and elem_type == "C3D10":
                parts = [p.strip() for p in s.split(",")]
                eid = int(parts[0])
                conn = [int(p) for p in parts[1:11]]
                elements[eid] = conn
    return nodes, elements


def dist_to_fillet(x, z):
    if z >= 0:
        cx, cz = TOP_CENTER
    else:
        cx, cz = BOTTOM_CENTER
    return abs(np.sqrt((x - cx) ** 2 + (z - cz) ** 2) - FILLET_R)


def bucket_counts(nodes, elements, fillet_size, label):
    immediate_thresh = 3 * fillet_size

    node_ids = np.array(list(nodes.keys()))
    coords = np.array([nodes[n] for n in node_ids])
    node_dist = np.array([dist_to_fillet(x, z) for x, y, z in coords])

    n_immediate = np.sum(node_dist < immediate_thresh)
    n_transition = np.sum((node_dist >= immediate_thresh) & (node_dist <= TRANSITION_WIDTH))
    n_far = np.sum(node_dist > TRANSITION_WIDTH)
    n_total = len(node_ids)

    elem_ids = np.array(list(elements.keys()))
    conn_array = np.array([elements[e] for e in elem_ids])
    corners = conn_array[:, 0:4]
    max_node_id = max(nodes.keys())
    node_coords_arr = np.zeros((max_node_id + 1, 3))
    for nid, c in nodes.items():
        node_coords_arr[nid] = c
    centroids = node_coords_arr[corners].mean(axis=1)
    elem_dist = np.array([dist_to_fillet(x, z) for x, y, z in centroids])

    e_immediate = np.sum(elem_dist < immediate_thresh)
    e_transition = np.sum((elem_dist >= immediate_thresh) & (elem_dist <= TRANSITION_WIDTH))
    e_far = np.sum(elem_dist > TRANSITION_WIDTH)
    e_total = len(elem_ids)

    print(f"\n=== {label} (FILLET_SIZE={fillet_size}mm, immediate band < {immediate_thresh:.3f}mm) ===")
    print(f"NODES  total={n_total}")
    print(f"  immediate (<{immediate_thresh:.3f}mm):        {n_immediate:>8} ({100*n_immediate/n_total:5.2f}%)")
    print(f"  transition ({immediate_thresh:.3f}-{TRANSITION_WIDTH}mm): {n_transition:>8} ({100*n_transition/n_total:5.2f}%)")
    print(f"  far-field (>{TRANSITION_WIDTH}mm):             {n_far:>8} ({100*n_far/n_total:5.2f}%)")
    print(f"ELEMENTS total={e_total}")
    print(f"  immediate (<{immediate_thresh:.3f}mm):        {e_immediate:>8} ({100*e_immediate/e_total:5.2f}%)")
    print(f"  transition ({immediate_thresh:.3f}-{TRANSITION_WIDTH}mm): {e_transition:>8} ({100*e_transition/e_total:5.2f}%)")
    print(f"  far-field (>{TRANSITION_WIDTH}mm):             {e_far:>8} ({100*e_far/e_total:5.2f}%)")

    return {
        "n_total": n_total, "n_immediate": n_immediate, "n_transition": n_transition, "n_far": n_far,
        "e_total": e_total, "e_immediate": e_immediate, "e_transition": e_transition, "e_far": e_far,
    }


print("Reading fine mesh (0.075mm)...")
nodes_fine, elems_fine = read_mesh(
    "01-aerospace-mounting-bracket/mesh/fillet_study/convergence/level_fine/mesh_fillet_raw.inp"
)
result_fine = bucket_counts(nodes_fine, elems_fine, 0.075, "FINE (0.075mm)")

print("\nReading intermediate2 mesh (0.06mm)...")
nodes_i2, elems_i2 = read_mesh(
    "01-aerospace-mounting-bracket/mesh/fillet_study/convergence/level_intermediate2/mesh_fillet_raw.inp"
)
result_i2 = bucket_counts(nodes_i2, elems_i2, 0.06, "INTERMEDIATE2 (0.06mm)")

print("\n=== GROWTH ATTRIBUTION: fine (0.075mm) -> intermediate2 (0.06mm) ===")
for key, label in [("n_immediate", "NODES immediate"), ("n_transition", "NODES transition"),
                    ("n_far", "NODES far-field"), ("n_total", "NODES TOTAL")]:
    old = result_fine[key]
    new = result_i2[key]
    delta = new - old
    pct_of_total_growth = 100 * delta / (result_i2["n_total"] - result_fine["n_total"])
    print(f"{label:22s}: {old:>8} -> {new:>8}  (delta={delta:+7d}, {pct_of_total_growth:5.1f}% of total node growth)")

print()
for key, label in [("e_immediate", "ELEMENTS immediate"), ("e_transition", "ELEMENTS transition"),
                    ("e_far", "ELEMENTS far-field"), ("e_total", "ELEMENTS TOTAL")]:
    old = result_fine[key]
    new = result_i2[key]
    delta = new - old
    pct_of_total_growth = 100 * delta / (result_i2["e_total"] - result_fine["e_total"])
    print(f"{label:22s}: {old:>8} -> {new:>8}  (delta={delta:+7d}, {pct_of_total_growth:5.1f}% of total element growth)")
