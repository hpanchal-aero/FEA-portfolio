#!/usr/bin/env python3
"""Survey an IGES/STEP geometry with the Gmsh/OpenCASCADE kernel: extent, unit scale, surface statistics.

usage: inspect_geometry.py <file.igs> [reference_semispan_m = 29.382]

Read-only. The unit of the imported model is not assumed: ratios of the model's y-extent to the reference
semispan are printed so the scale (m, mm or inch) can be identified.
"""
import resource
import sys
import time

import gmsh
import numpy as np

path = sys.argv[1]
ref = float(sys.argv[2]) if len(sys.argv) > 2 else 29.382

gmsh.initialize()
gmsh.model.add("survey")
t0 = time.time()
dimtags = gmsh.model.occ.importShapes(path, highestDimOnly=True)
gmsh.model.occ.synchronize()
t_imp = time.time() - t0

print(f"\nfile: {path}")
print(f"import + synchronize: {t_imp:.1f} s")
print("imported top-level entities by dimension:", {d: sum(1 for x in dimtags if x[0] == d) for d in (0, 1, 2, 3)})
print("model entities by dimension:", {d: len(gmsh.model.getEntities(d)) for d in (0, 1, 2, 3)})

surf = gmsh.model.getEntities(2)
if not surf:
    sys.exit("ABORT: no surfaces in the model")
bb = np.array([gmsh.model.getBoundingBox(2, t) for _, t in surf])          # xmin ymin zmin xmax ymax zmax
area = np.array([gmsh.model.occ.getMass(2, t) for _, t in surf])
tag = np.array([t for _, t in surf])
lo, hi = bb[:, :3].min(axis=0), bb[:, 3:].max(axis=0)
print(f"\nmodel bounding box (model units):  x {lo[0]:.4f} .. {hi[0]:.4f}   y {lo[1]:.4f} .. {hi[1]:.4f}   z {lo[2]:.4f} .. {hi[2]:.4f}")
print(f"extent: dx {hi[0] - lo[0]:.4f}  dy {hi[1] - lo[1]:.4f}  dz {hi[2] - lo[2]:.4f}")
print(f"ratio ymax / {ref} m = {hi[1] / ref:.4f}     ratio (ymax - ymin) / {ref} m = {(hi[1] - lo[1]) / ref:.4f}")
print("  (1 -> metres, ~39.37 -> inches, ~1000 -> millimetres, ~2 -> full span in metres)")
print(f"surfaces: {len(surf)}   total area {area.sum():.4e} (model units^2)")
print(f"surfaces touching y <= ymin + 1% of dy: {int((bb[:, 1] <= lo[1] + 0.01 * (hi[1] - lo[1])).sum())}")

print("\n12 largest surfaces by area:")
print("  tag      area        xmin      xmax      ymin      ymax      zmin      zmax")
for i in np.argsort(-area)[:12]:
    b = bb[i]
    print(f"{tag[i]:5d} {area[i]:12.4e} " + " ".join(f"{v:9.3f}" for v in (b[0], b[3], b[1], b[4], b[2], b[5])))

for axis, name in ((1, "y"), (0, "x"), (2, "z")):
    c = 0.5 * (bb[:, axis] + bb[:, axis + 3])
    edges = np.linspace(c.min(), c.max(), 9)
    cnt, _ = np.histogram(c, bins=edges)
    ar, _ = np.histogram(c, bins=edges, weights=area)
    print(f"\nsurface-centre {name}: bin edge | surfaces | area share")
    for k in range(8):
        print(f"  {edges[k]:10.3f} .. {edges[k + 1]:10.3f} | {cnt[k]:5d} | {100 * ar[k] / area.sum():5.1f} %")

print(f"\npeak resident memory {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024:.0f} MiB")
gmsh.finalize()
