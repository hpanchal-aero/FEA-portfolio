"""
Stage 5 -- PyVista visual verification, v2.
Fixes the blank zoom render from v1 and adds targeted, tag-identified
highlighting of the corner-bite void boundary (R3 arc vs gusset chord)
and the toe fillet, so we confirm geometry directly rather than by
eyeballing a flat-grey STL.
"""

import gmsh
import pyvista as pv
import numpy as np
import os

INPUT_BREP = "mesh/full_bracket_study/frame_gusset_toefillet_check.brep"
FIG_DIR = "figures"
os.makedirs(FIG_DIR, exist_ok=True)

gmsh.initialize()
gmsh.model.add("toe_fillet_visual_v2")
gmsh.model.occ.importShapes(INPUT_BREP)
gmsh.model.occ.synchronize()

# --- Identify target surfaces by bbox/area, same logic as prior diagnostics ---
surfaces = gmsh.model.getEntities(dim=2)
TOL = 0.05

arc_tag = chord_tag = toe_fillet_tag = None
for (dim, tag) in surfaces:
    bb = gmsh.model.occ.getBoundingBox(dim, tag)
    area = gmsh.model.occ.getMass(2, tag)
    xr = (round(bb[0], 2), round(bb[3], 2))
    zr = (round(bb[2], 2), round(bb[5], 2))
    if xr == (55.0, 58.0) and zr == (2.0, 5.0):
        if abs(area - 188.4956) < 0.5:
            arc_tag = tag
        elif abs(area - 169.7056) < 0.5:
            chord_tag = tag
    if abs(bb[0] - 37.59) < 0.1 and abs(bb[3] - 38.29) < 0.1:
        if abs(area - 31.4159) < 0.5:
            toe_fillet_tag = tag

print(f"Identified R3 base-fillet arc surface: tag {arc_tag}")
print(f"Identified gusset At-Bt chord surface: tag {chord_tag}")
print(f"Identified R1 toe-fillet surface: tag {toe_fillet_tag}")

if arc_tag is None or chord_tag is None or toe_fillet_tag is None:
    print("ERROR: one or more target surfaces not found -- cannot proceed with highlighting.")
    gmsh.finalize()
    raise SystemExit(1)

# --- Mesh, finer than v1 to resolve the small features properly ---
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.8)
gmsh.option.setNumber("Mesh.MeshSizeMin", 0.08)
gmsh.model.mesh.generate(2)

def extract_surface_mesh(tag):
    """Pull triangle mesh for a single gmsh surface into a PyVista PolyData."""
    elem_types, elem_tags, node_tags_flat = gmsh.model.mesh.getElements(2, tag)
    all_verts = []
    faces = []
    node_id_map = {}
    next_idx = 0
    for et, etags, ntags in zip(elem_types, elem_tags, node_tags_flat):
        n_per_elem = 3 if et == 2 else None  # type 2 = 3-node triangle
        if n_per_elem is None:
            continue
        ntags = np.array(ntags).reshape(-1, n_per_elem)
        for tri in ntags:
            face_idx = []
            for nt in tri:
                nt = int(nt)
                if nt not in node_id_map:
                    coord, _, _, _ = gmsh.model.mesh.getNode(nt)
                    all_verts.append(coord)
                    node_id_map[nt] = next_idx
                    next_idx += 1
                face_idx.append(node_id_map[nt])
            faces.append([3] + face_idx)
    verts = np.array(all_verts)
    faces_arr = np.hstack(faces)
    return pv.PolyData(verts, faces_arr)

arc_mesh = extract_surface_mesh(arc_tag)
chord_mesh = extract_surface_mesh(chord_tag)
toe_fillet_mesh = extract_surface_mesh(toe_fillet_tag)
print(f"Arc mesh: {arc_mesh.n_points} pts, {arc_mesh.n_cells} tris")
print(f"Chord mesh: {chord_mesh.n_points} pts, {chord_mesh.n_cells} tris")
print(f"Toe fillet mesh: {toe_fillet_mesh.n_points} pts, {toe_fillet_mesh.n_cells} tris")

# --- Full solid for context (grey, semi-transparent) ---
STL_OUT = "mesh/full_bracket_study/toe_fillet_visual_check_v2.stl"
gmsh.write(STL_OUT)
full_mesh = pv.read(STL_OUT)
gmsh.finalize()

# --- Figure A: corner-void region close-up (oblique 3D) ---
pl = pv.Plotter(off_screen=True, window_size=(1400, 1000))
pl.add_mesh(full_mesh, color="lightsteelblue", opacity=0.35, show_edges=False)
pl.add_mesh(arc_mesh, color="red", show_edges=True, edge_color="darkred", label="R3 base fillet arc")
pl.add_mesh(chord_mesh, color="blue", show_edges=True, edge_color="darkblue", label="Gusset At-Bt chord")
pl.add_mesh(toe_fillet_mesh, color="green", show_edges=True, edge_color="darkgreen", label="R1 toe fillet")
pl.add_legend()
pl.camera_position = [(90, -40, 40), (48, 20, 6), (0, 0, 1)]
pl.add_text("Stage 5 - Corner void (red vs blue) and toe fillet (green)", font_size=11)
outA = os.path.join(FIG_DIR, "stage5_v2_corner_void_3d.png")
pl.screenshot(outA)
pl.close()
print(f"Saved: {outA}")

# --- Figure B: tight top-down view of just the corner-void region ---
pl2 = pv.Plotter(off_screen=True, window_size=(1400, 1000))
pl2.add_mesh(full_mesh, color="lightsteelblue", opacity=0.25, show_edges=False)
pl2.add_mesh(arc_mesh, color="red", show_edges=True, edge_color="darkred")
pl2.add_mesh(chord_mesh, color="blue", show_edges=True, edge_color="darkblue")
pl2.camera_position = [(56.5, -60, 3.5), (56.5, 20, 3.5), (0, 0, 1)]
pl2.camera.parallel_projection = True
pl2.camera.parallel_scale = 8
pl2.add_text("Stage 5 - Corner void detail (arc=red, chord=blue)", font_size=11)
outB = os.path.join(FIG_DIR, "stage5_v2_corner_void_detail.png")
pl2.screenshot(outB)
pl2.close()
print(f"Saved: {outB}")

# --- Figure C: fixed toe-fillet zoom (corrected camera vs v1's blank render) ---
pl3 = pv.Plotter(off_screen=True, window_size=(1400, 1000))
pl3.add_mesh(full_mesh, color="lightsteelblue", opacity=0.25, show_edges=False)
pl3.add_mesh(toe_fillet_mesh, color="green", show_edges=True, edge_color="darkgreen")
pl3.camera_position = [(38, -60, 2.2), (38, 20, 2.2), (0, 0, 1)]
pl3.camera.parallel_projection = True
pl3.camera.parallel_scale = 4
pl3.add_text("Stage 5 - Toe fillet detail (green)", font_size=11)
outC = os.path.join(FIG_DIR, "stage5_v2_toe_fillet_detail.png")
pl3.screenshot(outC)
pl3.close()
print(f"Saved: {outC}")

print("\nDone. Inspect:")
print(f"  A. {outA} -- oblique 3D, shows how red/blue/green surfaces relate to the whole corner")
print(f"  B. {outB} -- tight view: should show a visible thin gap/notch between red (arc)")
print(f"     and blue (chord) if the corner-bite void is real and correctly located")
print(f"  C. {outC} -- should show a smooth green quarter-round blend at the toe, no spikes")
