"""
Standalone diagnostic (does not modify build_stage2_ribs.py). Fuses
ribs onto the spar ONE AT A TIME, checking the resulting volume count
after each addition, to isolate exactly which rib (if any) breaks the
assembly into >1 volume.
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
    gmsh.model.add("progressive_fuse")
    occ = gmsh.model.occ

    spar_tags = occ.importShapes(STEP_FILE_SPAR)
    current = [t for d, t in spar_tags if d == 3][0]
    occ.synchronize()

    for z in RIB_STATIONS:
        rib = build_rib(occ, z)
        fused, _ = occ.fuse([(3, current)], [(3, rib)])
        occ.synchronize()
        n_vols = len(fused)
        print(f"After fusing rib z={z:.0f}mm: {n_vols} volume(s) "
              f"{'<-- BROKE HERE' if n_vols != 1 else 'OK'}")
        if n_vols != 1:
            for d, t in fused:
                bbox = gmsh.model.getBoundingBox(d, t)
                vol = occ.getMass(d, t)
                print(f"    volume tag {t}: mass={vol:.4f} mm^3, "
                      f"bbox z=[{bbox[2]:.4f},{bbox[5]:.4f}]")
            break
        current = fused[0][1]

    gmsh.finalize()


if __name__ == "__main__":
    main()
