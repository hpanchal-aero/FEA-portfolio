"""
Stage 4b (right-angle frame WITH gusset) - geometry/mesh
visualization via PyVista: full-model overview plus a mid-width
cross-section slice, to confirm the gusset shape, its attachment to
both flanges, and the (expected, now confirmed) absence of the old
Stage 4a reentrant corner as a distinct surface feature.
"""

import re
import numpy as np
import pyvista as pv

MESH_FILE = "01-aerospace-mounting-bracket/mesh/frame_gusset_study/mesh_frame_gusset_raw.inp"
NSET_FILE = "01-aerospace-mounting-bracket/mesh/frame_gusset_study/nsets_frame_gusset.inp"
OUT_PNG_OVERVIEW = "01-aerospace-mounting-bracket/figures/stage4b_gusset_mesh_overview.png"
OUT_PNG_SLICE = "01-aerospace-mounting-bracket/figures/stage4b_gusset_corner_slice.png"


def read_nodes_and_c3d10(path):
    nodes, elements = {}, []
    mode, elem_type = None, None
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("**"):
                continue
            if s.upper().startswith("*NODE"):
                mode = "NODE"; continue
            if s.upper().startswith("*ELEMENT"):
                mode = "ELEMENT"
                m = re.search(r"TYPE\s*=\s*(\w+)", s, re.IGNORECASE)
                elem_type = m.group(1).upper() if m else None
                continue
            if s.startswith("*"):
                mode = None; continue
            if mode == "NODE":
                parts = s.split(",")
                nodes[int(parts[0])] = (float(parts[1]), float(parts[2]), float(parts[3]))
            elif mode == "ELEMENT" and elem_type == "C3D10":
                parts = [p.strip() for p in s.split(",")]
                elements.append([int(p) for p in parts[1:11]])
    return nodes, elements


def read_nset(path, name):
    node_ids = []
    with open(path) as f:
        lines = f.readlines()
    capture = False
    for line in lines:
        s = line.strip()
        if s.upper().startswith("*NSET"):
            capture = (f"NSET={name}".upper() in s.upper().replace(" ", ""))
            continue
        if s.startswith("*"):
            capture = False; continue
        if capture:
            parts = [p.strip() for p in s.split(",") if p.strip()]
            node_ids.extend(int(p) for p in parts)
    return set(node_ids)


print(f"Reading {MESH_FILE}")
nodes, elements = read_nodes_and_c3d10(MESH_FILE)
print(f"Nodes: {len(nodes)}, C3D10 elements: {len(elements)}")

root_nodes = read_nset(NSET_FILE, "root")
tip_nodes = read_nset(NSET_FILE, "tip")

node_id_list = sorted(nodes.keys())
id_to_idx = {nid: i for i, nid in enumerate(node_id_list)}
points = np.array([nodes[nid] for nid in node_id_list])

cells = []
for conn in elements:
    idx = [id_to_idx[nd] for nd in conn[0:4]]
    cells.append([4] + idx)
cells = np.array(cells).flatten()
cell_types = np.full(len(elements), pv.CellType.TETRA)
grid = pv.UnstructuredGrid(cells, cell_types, points)

# ---- Overview ----
plotter = pv.Plotter(off_screen=True, window_size=(1400, 1000))
plotter.add_mesh(grid.extract_surface(algorithm='dataset_surface'), color="lightsteelblue",
                  show_edges=True, edge_color="gray", opacity=0.9)
root_pts = np.array([nodes[nid] for nid in root_nodes])
tip_pts = np.array([nodes[nid] for nid in tip_nodes])
plotter.add_points(root_pts, color="red", point_size=6, render_points_as_spheres=True,
                    label="Root (fixed, x=0)")
plotter.add_points(tip_pts, color="green", point_size=6, render_points_as_spheres=True,
                    label="Tip (loaded, z=60)")
plotter.add_legend()
plotter.add_axes()
plotter.camera_position = [(150, -120, 100), (30, 20, 20), (0, 0, 1)]
plotter.add_text("Stage 4b: Right-Angle Frame WITH Gusset", font_size=14)
plotter.screenshot(OUT_PNG_OVERVIEW)
print(f"Saved: {OUT_PNG_OVERVIEW}")

# ---- Cross-section slice at y=20 (mid-width), same view as the 4a corner slice ----
slice_mesh = grid.slice(normal=(0, 1, 0), origin=(30, 20, 20))

plotter2 = pv.Plotter(off_screen=True, window_size=(1200, 1200))
plotter2.add_mesh(slice_mesh, color="lightsteelblue", show_edges=True,
                   edge_color="black", line_width=1)
plotter2.add_axes()
plotter2.enable_parallel_projection()
plotter2.camera_position = [(45, -100, 12), (45, 20, 12), (0, 0, 1)]
plotter2.camera.parallel_scale = 30
plotter2.add_text("Stage 4b: cross-section at y=20mm - gusset triangle visible",
                   font_size=12)
plotter2.screenshot(OUT_PNG_SLICE)
print(f"Saved: {OUT_PNG_SLICE}")
