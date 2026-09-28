"""
Rebuilds the spar + 6 ribs via SEQUENTIAL pairwise fusion (confirmed to
give exactly 1 volume at every step, unlike the real script's single
6-way simultaneous fuse). Writes STEP (for later reuse) and STL (for
quick PyVista visualization only -- NOT the analysis mesh). Does not
touch build_stage2_ribs.py.
"""
import gmsh
import math
import os
import sys

SPAN = 1000.0
ROOT_CHORD, TIP_CHORD = 300.0, 150.0
SPAR_CHORD_FRAC = 0.30
RIB_STATIONS = [0.0, 200.0, 400.0, 600.0, 800.0, 1000.0]
RIB_T = 1.5
N_AIRFOIL_PTS = 40

STEP_FILE_SPAR = "mesh/stage1_spar/stage1_spar.step"
OUT_DIR = "mesh/stage2_ribs"
os.makedirs(OUT_DIR, exist_ok=True)


def chord_at(z):
    return ROOT_CHORD + (TIP_CHORD - ROOT_CHORD) * (z / SPAN)


def naca0012_half_thickness(x_over_c, t=0.12):
    return 5 * t * (0.2969 * math.sqrt(x_over_c) - 0.1260 * x_over_c
                     - 0.3516 * x_over_c**2 + 0.2843 * x_over_c**3
                     - 0.1015 * x_over_c**4)


def airfoil_outline(chord, x_shift):
    xs = [chord * (1 - math.cos(math.pi * i / (2 * (N_AIRFOIL_PTS - 1))))
          for i in range(N_AIRFOIL_PTS)]
    upper = [(x - x_shift, chord * naca0012_half_thickness(x / chord)) for x in xs]
    lower = [(x - x_shift, -chord * naca0012_half_thickness(x / chord)) for x in reversed(xs[:-1])]
    return upper + lower[:-1]


def build_rib(occ, z):
    c = chord_at(z)
    x_shift = SPAR_CHORD_FRAC * c
    outline = airfoil_outline(c, x_shift)
    pt_tags = [occ.addPoint(x, y, z - RIB_T / 2.0) for x, y in outline]
    line_tags = [occ.addLine(pt_tags[i], pt_tags[(i + 1) % len(pt_tags)])
                 for i in range(len(pt_tags))]
    loop = occ.addCurveLoop(line_tags)
    surf = occ.addPlaneSurface([loop])
    occ.synchronize()
    extruded = occ.extrude([(2, surf)], 0, 0, RIB_T)
    occ.synchronize()
    return [t for d, t in extruded if d == 3][0]


def main():
    gmsh.initialize()
    gmsh.model.add("visualize_progressive")
    occ = gmsh.model.occ

    spar_tags = occ.importShapes(STEP_FILE_SPAR)
    current = [t for d, t in spar_tags if d == 3][0]
    occ.synchronize()

    for z in RIB_STATIONS:
        rib = build_rib(occ, z)
        fused, _ = occ.fuse([(3, current)], [(3, rib)])
        occ.synchronize()
        assert len(fused) == 1, f"unexpected {len(fused)} volumes after z={z}"
        current = fused[0][1]

    final_vol = occ.getMass(3, current)
    print(f"Final sequentially-fused assembly volume: {final_vol:.4f} mm^3")

    step_path = os.path.join(OUT_DIR, "stage2_ribs_PREVIEW.step")
    gmsh.write(step_path)
    print(f"Saved: {step_path}")

    # Coarse surface mesh purely for visualization (NOT the analysis mesh)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 5.0)
    gmsh.model.mesh.generate(2)
    stl_path = os.path.join(OUT_DIR, "stage2_ribs_PREVIEW.stl")
    gmsh.write(stl_path)
    print(f"Saved: {stl_path}")

    gmsh.finalize()


if __name__ == "__main__":
    main()
