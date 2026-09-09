"""
Stage 4a (right-angle frame) - close-up cross-section view at the
reentrant corner, sliced at y=20 (mid-width), to directly confirm
the interior corner geometry is sharp (not accidentally filleted,
gapped, or malformed by the boolean union).
"""

import re
import numpy as np
import pyvista as pv

MESH_FILE = "01-aerospace-mounting-bracket/mesh/frame_study/mesh_frame_raw.inp"
OUT_PNG = "01-aerospace-mounting-bracket/figures/stage4_frame_corner_slice.png"


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


print(f"Reading {MESH_FILE}")
nodes, elements = read_nodes_and_c3d10(MESH_FILE)

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

slice_mesh = grid.slice(normal=(0, 1, 0), origin=(30, 20, 20))
print(f"Slice bounds: {slice_mesh.bounds}")
print(f"Slice n_points: {slice_mesh.n_points}")

plotter = pv.Plotter(off_screen=True, window_size=(1200, 1200))
plotter.add_mesh(slice_mesh, color="lightsteelblue", show_edges=True,
                  edge_color="black", line_width=1)

corner_pt = (58.0, 20.0, 2.0)
plotter.add_points(np.array([corner_pt]), color="orange", point_size=25,
                    render_points_as_spheres=True, label="Reentrant corner (x=58,z=2)")

plotter.add_legend()
plotter.add_axes()

# Explicit parallel-projection camera, looking along -y, focused on
# the corner region - avoids the view_xz()+tight() conflict that
# produced a blank render previously.
plotter.enable_parallel_projection()
plotter.camera_position = [(55, -100, 5), (55, 20, 5), (0, 0, 1)]
plotter.camera.parallel_scale = 15

plotter.add_text("Stage 4a: cross-section at y=20mm, reentrant corner close-up",
                  font_size=12)

plotter.screenshot(OUT_PNG)
print(f"Saved: {OUT_PNG}")
