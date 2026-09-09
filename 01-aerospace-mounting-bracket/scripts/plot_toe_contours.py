"""
Stage 4b - contour visualizations.
  1. Full-model sigma_xx contour (deformed, exaggerated), finest
     solved level (toe_intermediate, 0.10mm) - shows the gusset,
     overall stress pattern, and toe hot spot in context.
  2. Zoomed contour at the toe region, side-by-side across all three
     solved mesh levels, SAME fixed color scale - visually shows the
     hot spot intensifying/tightening with refinement (non-
     convergence), rather than stabilizing.
"""

import re
import numpy as np
import pyvista as pv

LEVELS = [
    ("Baseline (~0.5mm)", "."),
    ("toe_fine (0.15mm)", "toe_convergence/level_toe_fine"),
    ("toe_intermediate (0.10mm)", "toe_convergence/level_toe_intermediate"),
]


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


def build_grid(mesh_file, frd_file, scale=0.0):
    nodes, elements = read_nodes_and_c3d10(mesh_file)
    disp, stress = parse_frd_disp_stress(frd_file)

    node_id_list = sorted(nodes.keys())
    id_to_idx = {nid: i for i, nid in enumerate(node_id_list)}
    points = np.array([nodes[nid] for nid in node_id_list])
    disp_arr = np.array([disp.get(nid, (0, 0, 0)) for nid in node_id_list])
    sxx_arr = np.array([stress.get(nid, (0,) * 6)[0] for nid in node_id_list])
    # von Mises from the 6 stress components (CalculiX order: xx,yy,zz,xy,yz,zx)
    s6 = np.array([stress.get(nid, (0,) * 6) for nid in node_id_list])
    sxx_, syy_, szz_, sxy_, syz_, szx_ = [s6[:, k] for k in range(6)]
    vm_arr = np.sqrt(0.5 * ((sxx_ - syy_)**2 + (syy_ - szz_)**2 + (szz_ - sxx_)**2
                            + 6*(sxy_**2 + syz_**2 + szx_**2)))

    plot_points = points + scale * disp_arr if scale else points

    cells = []
    for conn in elements:
        idx = [id_to_idx[nd] for nd in conn[0:4]]
        cells.append([4] + idx)
    cells = np.array(cells).flatten()
    cell_types = np.full(len(elements), pv.CellType.TETRA)

    grid = pv.UnstructuredGrid(cells, cell_types, plot_points)
    grid["sigma_xx"] = sxx_arr
    grid["von_mises"] = vm_arr
    return grid


# ---- Figure 1: full-model contour, finest solved level, deformed ----
print("Building full-model contour (toe_intermediate level)...")
grid_full = build_grid(
    "toe_convergence/level_toe_intermediate/mesh_frame_gusset_raw.inp",
    "toe_convergence/level_toe_intermediate/analysis_frame_gusset.frd",
    scale=5.0,
)

plotter1 = pv.Plotter(off_screen=True, window_size=(1400, 1100))
plotter1.add_mesh(grid_full.extract_surface(algorithm='dataset_surface'),
                   scalars="von_mises", cmap="turbo", clim=[0, 200],
                   show_edges=False, scalar_bar_args={"title": "von Mises (MPa)"})
plotter1.camera_position = [(180, -150, 120), (30, 20, 30), (0, 0, 1)]
plotter1.add_text("Stage 4b: Full-model von Mises stress (5x deformed), toe_intermediate mesh\n"
                   "von Mises used (not sigma_xx) so both flanges show correctly - each bends about a different axis\n"
                   "Scale capped at 200 MPa to show the global pattern; toe hot spot saturates (see separate non-convergence figures)",
                   font_size=11)
plotter1.screenshot("../../figures/stage4b_gusset_full_contour.png")
print("Saved: stage4b_gusset_full_contour.png")

# ---- Figure 2: zoomed toe region, side-by-side across levels, SAME clim ----
print("\nBuilding zoomed toe contours across all three levels...")
plotter2 = pv.Plotter(off_screen=True, shape=(1, 3), window_size=(2100, 750))

CLIM = [-300, 0]  # fixed across all three, spans the full observed peak range

for i, (label, level_dir) in enumerate(LEVELS):
    mesh_file = f"{level_dir}/mesh_frame_gusset_raw.inp"
    frd_file = f"{level_dir}/analysis_frame_gusset.frd"
    print(f"  {label}...")
    grid = build_grid(mesh_file, frd_file, scale=0.0)

    # clip to a box around the toe: x in [30,45], z in [-3,10]
    clipped = grid.clip_box(bounds=(30, 45, 0, 40, -3, 10), invert=False)

    plotter2.subplot(0, i)
    plotter2.add_mesh(clipped.extract_surface(algorithm='dataset_surface'),
                       scalars="sigma_xx", cmap="coolwarm", clim=CLIM,
                       show_edges=True, edge_color="gray", line_width=0.5,
                       scalar_bar_args={"title": "sigma_xx (MPa)"} if i == 2 else None,
                       show_scalar_bar=(i == 2))
    plotter2.camera_position = [(38, -60, 15), (38, 20, 2), (0, 0, 1)]
    plotter2.add_text(label, font_size=11)

plotter2.screenshot("../../figures/stage4b_toe_contour_comparison.png")
print("Saved: stage4b_toe_contour_comparison.png")
