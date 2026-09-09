"""
Stage 4a (right-angle frame) - geometry/mesh visualization via PyVista.

Parses the raw C3D10 mesh directly (corner nodes only, for shape
visualization - not a stress plot), highlights the root (fixed) and
tip (loaded) NSETs, and marks the reentrant corner location.

Saves a static PNG to 01-aerospace-mounting-bracket/figures/.
"""

import re
import numpy as np
import pyvista as pv

MESH_FILE = "01-aerospace-mounting-bracket/mesh/frame_study/mesh_frame_raw.inp"
NSET_FILE = "01-aerospace-mounting-bracket/mesh/frame_study/nsets_frame.inp"
OUT_PNG = "01-aerospace-mounting-bracket/figures/stage4_frame_mesh_visualization.png"


def read_nodes_and_c3d10(path):
    nodes = {}
    elements = []
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
                elements.append(conn)
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
            capture = False
            continue
        if capture:
            parts = [p.strip() for p in s.split(",") if p.strip()]
            node_ids.extend(int(p) for p in parts)
    return set(node_ids)


print(f"Reading {MESH_FILE}")
nodes, elements = read_nodes_and_c3d10(MESH_FILE)
print(f"Nodes: {len(nodes)}, C3D10 elements: {len(elements)}")

root_nodes = read_nset(NSET_FILE, "root")
tip_nodes = read_nset(NSET_FILE, "tip")
print(f"root NSET: {len(root_nodes)} nodes, tip NSET: {len(tip_nodes)} nodes")

# Build a linear-tet UnstructuredGrid (corner nodes only - shape
# visualization, not a stress/field plot)
node_id_list = sorted(nodes.keys())
id_to_idx = {nid: i for i, nid in enumerate(node_id_list)}
points = np.array([nodes[nid] for nid in node_id_list])

cells = []
for conn in elements:
    corners = conn[0:4]
    idx = [id_to_idx[n] for n in corners]
    cells.append([4] + idx)
cells = np.array(cells).flatten()
cell_types = np.full(len(elements), pv.CellType.TETRA)

grid = pv.UnstructuredGrid(cells, cell_types, points)

surface = grid.extract_surface()

plotter = pv.Plotter(off_screen=True, window_size=(1400, 1000))
plotter.add_mesh(surface, color="lightsteelblue", show_edges=True,
                  edge_color="gray", opacity=0.9)

# Highlight root and tip node sets
root_pts = np.array([nodes[nid] for nid in root_nodes])
tip_pts = np.array([nodes[nid] for nid in tip_nodes])
plotter.add_points(root_pts, color="red", point_size=6, render_points_as_spheres=True,
                    label="Root (fixed, x=0)")
plotter.add_points(tip_pts, color="green", point_size=6, render_points_as_spheres=True,
                    label="Tip (loaded, z=60)")

# Mark the reentrant corner location (x=58, z=2, y=20 - mid-width)
corner_pt = np.array([[58.0, 20.0, 2.0]])
plotter.add_points(corner_pt, color="orange", point_size=20,
                    render_points_as_spheres=True, label="Reentrant corner")

plotter.add_legend()
plotter.add_axes()
plotter.camera_position = [(150, -120, 100), (30, 20, 20), (0, 0, 1)]
plotter.add_text("Stage 4a: Right-Angle Frame (sharp corner, no gusset)",
                  font_size=14)

plotter.screenshot(OUT_PNG)
print(f"Saved: {OUT_PNG}")
