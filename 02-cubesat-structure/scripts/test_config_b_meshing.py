"""
Config B assembly meshing feasibility test (mesh generation ONLY).

Usage: /usr/bin/python3 scripts/test_config_b_meshing.py <shape> <level_pct> [maxh ...]

Loads mesh/config_b_study/config_b_assembly_<shape>_<level>pct.step, generates
a linear-tet mesh at each given uniform size (default: just 3.0), and reports
node/element counts and intersection-warning presence. No .inp is written and
no solver is run.
"""
import sys
import time
from netgen.occ import OCCGeometry

if len(sys.argv) < 3:
    print("Usage: test_config_b_meshing.py <shape> <level_pct> [maxh ...]")
    sys.exit(1)
shape, level = sys.argv[1], sys.argv[2]
sizes = [float(a) for a in sys.argv[3:]] or [3.0]

step_file = f"mesh/config_b_study/config_b_assembly_{shape}_{level}pct.step"

for h in sizes:
    geo = OCCGeometry(step_file)
    shp = geo.shape
    n_solids = len(list(shp.solids))
    n_faces = len(list(shp.faces))
    print(f"Geometry: {step_file}  type={shp.type}  solids={n_solids}  faces={n_faces}")
    t0 = time.time()
    try:
        mesh = geo.GenerateMesh(maxh=h)
        dt = time.time() - t0
        n_pts = len(mesh.Points())
        n_tet = len(mesh.Elements3D())
        print(f"RESULT maxh={h}: OK, linear_nodes={n_pts}, tets={n_tet}, time={dt:.1f}s")
    except Exception as e:
        dt = time.time() - t0
        print(f"RESULT maxh={h}: FAILED after {dt:.1f}s, {type(e).__name__}: {e}")
