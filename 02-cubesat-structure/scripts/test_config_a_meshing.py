#!/usr/bin/python3
"""
Config A meshing feasibility test (mesh generation ONLY).

Loads the cleaned single-solid STEP, generates a linear-tet mesh at several
uniform sizes, and reports node/element counts. No .inp is written and no
solver is run, so this cannot exhaust memory in CalculiX.

Purpose: determine whether the failures seen with config_a_base_frame.brep
were caused by the stray z=50 section face in the compound, or by a genuine
thin-wall meshing problem.
"""
import time
from netgen.occ import OCCGeometry

STEP_FILE = "mesh/config_a_study/config_a_clean.step"
SIZES = [3.0, 2.0, 1.5]

for h in SIZES:
    geo = OCCGeometry(STEP_FILE)
    n_faces = len(list(geo.shape.faces))
    t0 = time.time()
    try:
        mesh = geo.GenerateMesh(maxh=h)
        dt = time.time() - t0
        n_pts = len(mesh.Points())
        n_tet = len(mesh.Elements3D())
        print(f"RESULT maxh={h}: OK, faces_in_geo={n_faces}, "
              f"linear_nodes={n_pts}, tets={n_tet}, time={dt:.1f}s")
    except Exception as e:
        dt = time.time() - t0
        print(f"RESULT maxh={h}: FAILED after {dt:.1f}s, faces_in_geo={n_faces}, "
              f"{type(e).__name__}: {e}")
