"""
Project 03, Stage 2 -- contour figures (PyVista, off-screen) for the
2.0mm mesh (the finer of the two converged levels).

Builds a VTK unstructured grid from the CalculiX .inp/.frd (corner
nodes only, linear tets -- for visualization; the C3D10 midside nodes
are not needed), attaches displacement and von Mises fields, saves it
as .vtu (for the interactive viewer), and renders:
  figures/contours/stage2_h2p0_uy_iso.png     (deformed, x SCALE)
  figures/contours/stage2_h2p0_uy_side.png    (deformed, x SCALE)
  figures/contours/stage2_h2p0_vonmises.png   (undeformed)
  figures/contours/stage2_h2p0_vonmises_root.png (root close-up)

Usage (conda env ccx223): python3 scripts/render_stage2_contours.py
"""
import os
import numpy as np
import pyvista as pv

JOB = "simulation/stage2_ribs_h2p0"
VTU = "simulation/stage2_ribs_h2p0_fields.vtu"
OUT_DIR = "figures/contours"
SCALE = 50.0  # deformation exaggeration, stated on every deformed figure


def parse_inp_mesh(path):
    nodes, elems = [], []
    mode = None
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("*"):
                u = s.upper()
                if u.startswith("*NODE") and not u.startswith("*NODE FILE"):
                    mode = "node"
                elif u.startswith("*ELEMENT"):
                    mode = "elem"
                else:
                    mode = None
                continue
            p = [t.strip() for t in s.split(",") if t.strip()]
            if mode == "node":
                nodes.append([float(v) for v in p[1:4]])
            elif mode == "elem":
                elems.append([int(v) for v in p[1:11]])
    return np.array(nodes), np.array(elems)


def read_frd_block(path, name, ncomp, n_nodes):
    data = None
    in_block = False
    with open(path) as f:
        for line in f:
            if line.startswith(" -4"):
                in_block = (line[5:13].strip() == name)
                if in_block:
                    data = np.zeros((n_nodes + 1, ncomp))
                continue
            if in_block:
                if line.startswith(" -3"):
                    in_block = False
                elif line.startswith(" -1"):
                    node = int(line[3:13])
                    for k in range(ncomp):
                        data[node, k] = float(line[13 + 12*k:25 + 12*k])
    if data is None:
        raise RuntimeError(f"block {name} not found in {path}")
    return data


def von_mises(s):
    sxx, syy, szz, sxy, syz, szx = (s[:, i] for i in range(6))
    return np.sqrt(0.5*((sxx-syy)**2 + (syy-szz)**2 + (szz-sxx)**2)
                   + 3*(sxy**2 + syz**2 + szx**2))


def build_grid():
    xyz, elems = parse_inp_mesh(JOB + ".inp")
    n = len(xyz)
    U = read_frd_block(JOB + ".frd", "DISP", 3, n)
    S = read_frd_block(JOB + ".frd", "STRESS", 6, n)
    vm = von_mises(S)
    ne = len(elems)
    cells = np.hstack([np.full((ne, 1), 4), elems[:, :4] - 1]).ravel()
    celltypes = np.full(ne, 10, dtype=np.uint8)  # VTK_TETRA
    grid = pv.UnstructuredGrid(cells, celltypes, xyz)
    grid.point_data["U"] = U[1:]
    grid.point_data["Uy_mm"] = U[1:, 1]
    grid.point_data["Umag_mm"] = np.linalg.norm(U[1:], axis=1)
    grid.point_data["von_Mises_MPa"] = vm[1:]
    return grid


def render(mesh, scalars, bar_title, clim, out, camera, size, label, outline=None):
    pl = pv.Plotter(off_screen=True, window_size=size)
    pl.add_mesh(mesh, scalars=scalars, cmap="turbo", clim=clim,
                scalar_bar_args={"title": bar_title, "vertical": True,
                                 "n_labels": 6, "fmt": "%.2f"})
    if outline is not None and outline.n_points > 0:
        pl.add_mesh(outline, color="black", line_width=1, opacity=0.4)
    pl.add_text(label, font_size=9, position="upper_left")
    pl.add_axes()
    if camera == "iso":
        pl.camera_position = "iso"
        pl.camera.azimuth = 30
        pl.camera.elevation = 20
    elif camera == "yz":
        pl.view_yz()
    pl.reset_camera()
    pl.screenshot(out)
    pl.close()
    print(f"Saved: {out}")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    grid = build_grid()
    grid.save(VTU)
    print(f"Grid: {grid.n_points} points, {grid.n_cells} tets; saved {VTU}")

    uy_all = grid.point_data["Uy_mm"]
    vm_all = grid.point_data["von_Mises_MPa"]
    p99 = float(np.percentile(vm_all, 99.0))
    print(f"Uy range: [{uy_all.min():.4f}, {uy_all.max():.4f}] mm")
    print(f"von Mises: max={vm_all.max():.3f} MPa, 99th pct={p99:.3f}, "
          f"median={np.median(vm_all):.3f}")

    surf = grid.extract_surface(algorithm="dataset_surface")
    warped = surf.copy()
    warped.points = surf.points + SCALE * surf.point_data["U"]
    outline = surf.extract_feature_edges(feature_angle=30.0)
    dlabel = (f"Stage 2 (spar + 6 ribs), 2.0mm mesh, Uy contour\n"
              f"deformation exaggerated x{SCALE:g}; black = undeformed edges")
    uy_clim = (float(uy_all.min()), float(uy_all.max()))

    render(warped, "Uy_mm", "Uy [mm]", uy_clim,
           os.path.join(OUT_DIR, "stage2_h2p0_uy_iso.png"),
           "iso", (1600, 1000), dlabel, outline)
    render(warped, "Uy_mm", "Uy [mm]", uy_clim,
           os.path.join(OUT_DIR, "stage2_h2p0_uy_side.png"),
           "yz", (800, 1500), dlabel, outline)

    vlabel = (f"Stage 2 (spar + 6 ribs), 2.0mm mesh, von Mises (undeformed)\n"
              f"color scale capped at 99th pct ({p99:.1f} MPa); peak {vm_all.max():.1f} MPa "
              f"is a local, unconverged corner value")
    render(surf, "von_Mises_MPa", "von Mises [MPa]", (0.0, p99),
           os.path.join(OUT_DIR, "stage2_h2p0_vonmises.png"),
           "iso", (1600, 1000), vlabel)

    root = surf.clip_box(bounds=(-30.0, 30.0, -22.0, 22.0, -1.0, 40.0), invert=False)
    if root.n_points == 0:
        print("WARNING: root clip is empty -- skipping root close-up")
    else:
        rmax = float(root.point_data["von_Mises_MPa"].max())
        rlabel = (f"Root region close-up (x within +/-30 mm, z < 40 mm), von Mises, 2.0mm mesh\n"
                  f"local max in view {rmax:.1f} MPa (unconverged corner peak)")
        render(root, "von_Mises_MPa", "von Mises [MPa]", (0.0, rmax),
               os.path.join(OUT_DIR, "stage2_h2p0_vonmises_root.png"),
               "iso", (1400, 900), rlabel)


if __name__ == "__main__":
    main()
