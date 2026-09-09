"""
Stage 4a - deformed-shape visualization, exaggerated, colored by
sigma_xx, to visually resolve the Flange A tension/compression sign
question and sanity-check the overall bending sense against the
-x tip load.
"""

import re
import numpy as np
import pyvista as pv

MESH_FILE = "01-aerospace-mounting-bracket/mesh/frame_study/mesh_frame_raw.inp"
FRD_FILE = "01-aerospace-mounting-bracket/mesh/frame_study/analysis_frame.frd"
OUT_PNG = "01-aerospace-mounting-bracket/figures/stage4_frame_deformed.png"

SCALE = 5.0  # deformation exaggeration


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


def parse_frd_disp_stress(path):
    disp, stress = {}, {}
    with open(path) as f:
        lines = f.readlines()
    i, n = 0, len(lines)
    while i < n:
        line = lines[i]
        if "DISP" in line and line.startswith(" -4"):
            i += 1
            while i < n and lines[i].startswith(" -5"):
                i += 1
            while i < n and lines[i].startswith(" -1"):
                rec = lines[i].rstrip("\n")
                nid = int(rec[3:13])
                disp[nid] = (float(rec[13:25]), float(rec[25:37]), float(rec[37:49]))
                i += 1
            continue
        if "STRESS" in line and line.startswith(" -4"):
            i += 1
            while i < n and lines[i].startswith(" -5"):
                i += 1
            while i < n and lines[i].startswith(" -1"):
                rec = lines[i].rstrip("\n")
                nid = int(rec[3:13])
                vals = [float(rec[13 + 12*k:25 + 12*k]) for k in range(6)]
                stress[nid] = tuple(vals)
                i += 1
            continue
        i += 1
    return disp, stress


print("Reading mesh + results...")
nodes, elements = read_nodes_and_c3d10(MESH_FILE)
disp, stress = parse_frd_disp_stress(FRD_FILE)

node_id_list = sorted(nodes.keys())
id_to_idx = {nid: i for i, nid in enumerate(node_id_list)}
points = np.array([nodes[nid] for nid in node_id_list])
disp_arr = np.array([disp.get(nid, (0, 0, 0)) for nid in node_id_list])
sxx_arr = np.array([stress.get(nid, (0,)*6)[0] for nid in node_id_list])
szz_arr = np.array([stress.get(nid, (0,)*6)[2] for nid in node_id_list])

deformed_points = points + SCALE * disp_arr

cells = []
for conn in elements:
    idx = [id_to_idx[nd] for nd in conn[0:4]]
    cells.append([4] + idx)
cells = np.array(cells).flatten()
cell_types = np.full(len(elements), pv.CellType.TETRA)

grid_orig = pv.UnstructuredGrid(cells, cell_types, points)
grid_def = pv.UnstructuredGrid(cells, cell_types, deformed_points)
grid_def["sigma_xx"] = sxx_arr
grid_def["sigma_zz"] = szz_arr

plotter = pv.Plotter(off_screen=True, window_size=(1400, 1100))
plotter.add_mesh(grid_orig.extract_surface(), color="lightgray", opacity=0.25,
                  label="Undeformed")
surf_def = grid_def.extract_surface()
plotter.add_mesh(surf_def, scalars="sigma_zz", cmap="coolwarm",
                  show_edges=False, label=f"Deformed ({SCALE}x), colored by sigma_zz")

plotter.add_legend()
plotter.add_axes()
plotter.camera_position = [(180, -150, 120), (30, 20, 30), (0, 0, 1)]
plotter.add_text(f"Stage 4a: deformed shape ({SCALE}x exaggeration), sigma_zz (MPa)",
                  font_size=13)
plotter.screenshot(OUT_PNG)
print(f"Saved: {OUT_PNG}")
