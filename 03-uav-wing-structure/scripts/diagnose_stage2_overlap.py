"""
Standalone diagnostic (read-only w.r.t. real build script) -- does NOT
modify build_stage2_ribs.py. Imports the Stage 1 spar, builds all 6
ribs exactly as build_stage2_ribs.py does, then computes the actual
boolean intersection volume between the spar and EACH rib individually
(on copies, so nothing is consumed) to find which rib(s), if any, have
zero or near-zero geometric overlap with the spar.
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


def main():
    gmsh.initialize()
    gmsh.model.add("diagnose")
    occ = gmsh.model.occ

    spar_tags = occ.importShapes(STEP_FILE_SPAR)
    spar_vol = [t for d, t in spar_tags if d == 3][0]
    occ.synchronize()
    spar_bbox = gmsh.model.getBoundingBox(3, spar_vol)
    print(f"SPAR bbox: x=[{spar_bbox[0]:.4f},{spar_bbox[3]:.4f}] "
          f"y=[{spar_bbox[1]:.4f},{spar_bbox[4]:.4f}] "
          f"z=[{spar_bbox[2]:.4f},{spar_bbox[5]:.4f}]")
    print()

    for z in RIB_STATIONS:
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
        rib_vol = [t for d, t in extruded if d == 3][0]

        rib_bbox = gmsh.model.getBoundingBox(3, rib_vol)
        print(f"Rib z={z:.0f}mm bbox: x=[{rib_bbox[0]:.4f},{rib_bbox[3]:.4f}] "
              f"y=[{rib_bbox[1]:.4f},{rib_bbox[4]:.4f}] "
              f"z=[{rib_bbox[2]:.4f},{rib_bbox[5]:.4f}]")

        # non-destructive: copy spar and this rib, intersect the copies
        spar_copy = occ.copy([(3, spar_vol)])
        rib_copy = occ.copy([(3, rib_vol)])
        occ.synchronize()
        try:
            inter, _ = occ.intersect(spar_copy, rib_copy, removeObject=True, removeTool=True)
            occ.synchronize()
            if len(inter) == 0:
                print(f"  -> INTERSECTION VOLUME: 0 (NO OVERLAP) <-- PROBLEM RIB")
            else:
                vol = sum(occ.getMass(3, t) for d, t in inter)
                print(f"  -> intersection volume with spar: {vol:.4f} mm^3")
        except Exception as e:
            print(f"  -> intersection FAILED: {e}  <-- PROBLEM RIB")
        print()

    gmsh.finalize()


if __name__ == "__main__":
    main()
