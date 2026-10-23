#!/usr/bin/env python3
"""Registration test: distance from CalculiX wingbox nodes to the faces of an IGES aerodynamic surface.

usage: check_wbt_registration.py <aero.igs> <deck.inp> [factor = 25.4] [n_nodes = 200]

Node coordinates [m] are multiplied by `factor` to express them in the units of the model as imported by
Gmsh/OpenCASCADE. Only nodes with y >= 3.5 m (outside the fuselage half-width) are used. A random subset
(fixed seed) is projected on every face. Distances are reported in metres (model distance / factor).
Control: the same nodes shifted by +0.5 m in z. Read-only.
"""
import sys
import time

import gmsh
import numpy as np

igs, deck = sys.argv[1], sys.argv[2]
factor = float(sys.argv[3]) if len(sys.argv) > 3 else 25.4
n = int(sys.argv[4]) if len(sys.argv) > 4 else 200

pts, kw = [], None
for ln in open(deck):
    s = ln.strip()
    if not s or s.startswith("**"):
        continue
    if s.startswith("*"):
        kw = s.split(",")[0].upper()
        continue
    if kw == "*NODE":
        t = s.split(",")
        if len(t) >= 4:
            pts.append([float(t[1]), float(t[2]), float(t[3])])
pts = np.array(pts)
out = pts[pts[:, 1] >= 3.5]
print(f"deck nodes {len(pts)}; with y >= 3.5 m: {len(out)}")
rng = np.random.default_rng(0)
sel = out[rng.choice(len(out), size=min(n, len(out)), replace=False)]

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("reg")
gmsh.model.occ.importShapes(igs, highestDimOnly=True)
gmsh.model.occ.synchronize()
faces = [t for _, t in gmsh.model.getEntities(2)]
print(f"faces in {igs}: {len(faces)};  factor {factor}")


def dist_m(p_m):
    p = np.asarray(p_m) * factor
    best = np.inf
    for t in faces:
        c, _ = gmsh.model.getClosestPoint(2, t, list(p))
        best = min(best, float(np.linalg.norm(np.array(c) - p)))
    return best / factor


def report(label, d):
    d = np.array(d)
    q = np.percentile(d, [5, 25, 50, 75, 95])
    print(f"{label}: n={len(d)}  pct5/25/50/75/95 [m] = " + " / ".join(f"{v:.4f}" for v in q) + f"   max {d.max():.4f}")
    for tol in (0.005, 0.02, 0.1):
        print(f"    fraction within {tol * 1000:5.0f} mm: {100 * np.mean(d <= tol):5.1f} %")


t0 = time.time()
d0 = [dist_m(p) for p in sel]
report("deck nodes vs aerodynamic surface", d0)
print(f"  elapsed {time.time() - t0:.0f} s")
sub = sel[:60]
report("control: same nodes shifted +0.5 m in z", [dist_m(p + np.array([0.0, 0.0, 0.5])) for p in sub])
print(f"total elapsed {time.time() - t0:.0f} s")
gmsh.finalize()
