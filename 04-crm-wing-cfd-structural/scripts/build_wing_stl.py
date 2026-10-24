#!/usr/bin/env python3
"""Full-span closed wing surface and STL for snappyHexMesh (decision D11, option A).

usage: build_wing_stl.py [input.brep = geometry/wing_alone.brep] [output_prefix = geometry/wing_full] [sewing_tolerance_m = 0.0015]

The half-wing BREP (12 faces, model units = metres x 25.4, free edges only on the y = 0 root loop) is mirrored in the plane y = 0,
the two halves are sewn at the root, a triangular surface mesh is generated and written as a binary STL in metres with outward normals.
Sizing: 8 mm within 10 mm of the spanwise edge curves (leading edge, trailing edge, tip), growing linearly to 150 mm at 0.5 m
distance (Distance + Threshold fields); curvature-based refinement (24 elements per 2 pi); maximum 150 mm.
"""
import os
import sys

import gmsh
import numpy as np

src = sys.argv[1] if len(sys.argv) > 1 else "geometry/wing_alone.brep"
prefix = sys.argv[2] if len(sys.argv) > 2 else "geometry/wing_full"
tol_m = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0015
S = 1.0 / 25.4
SIZE_MAX, EDGE_MIN, SIZE_FLOOR, CURV = 0.150, 0.008, 0.004, 24
D_MIN, D_MAX = 0.010, 0.500

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("half")
gmsh.model.occ.importShapes(src, highestDimOnly=True)
gmsh.model.occ.synchronize()
faces = gmsh.model.getEntities(2)
if len(faces) != 12:
    sys.exit(f"ABORT: expected 12 faces in {src}, found {len(faces)}")
copy = gmsh.model.occ.copy(faces)
gmsh.model.occ.mirror(copy, 0, 1, 0, 0)
gmsh.model.occ.synchronize()
print(f"half wing faces {len(faces)}; after mirroring {len(gmsh.model.getEntities(2))}")
gmsh.write(prefix + "_unsewn.brep")
gmsh.clear()

for k, v in (("Geometry.OCCSewFaces", 1), ("Geometry.OCCFixDegenerated", 0), ("Geometry.OCCFixSmallEdges", 0),
             ("Geometry.OCCFixSmallFaces", 0), ("Geometry.OCCMakeSolids", 0), ("Geometry.Tolerance", tol_m / S)):
    gmsh.option.setNumber(k, v)
gmsh.model.add("full")
gmsh.model.occ.importShapes(prefix + "_unsewn.brep", highestDimOnly=True)
gmsh.model.occ.synchronize()
fcs = gmsh.model.getEntities(2)
crv = [t for _, t in gmsh.model.getEntities(1)]
free, multi = [], []
for c in crv:
    up, _ = gmsh.model.getAdjacencies(1, c)
    if len(up) == 1:
        free.append(c)
    elif len(up) > 2:
        multi.append(c)
area = sum(gmsh.model.occ.getMass(2, t) for _, t in fcs) * S * S
print(f"sewn at {tol_m * 1000:.2f} mm: faces {len(fcs)}, curves {len(crv)}, free edges {len(free)}, "
      f"edges shared by 3+ faces {len(multi)}, CAD-kernel area {area:.3f} m2")
for c in free[:10]:
    bb = [v * S for v in gmsh.model.getBoundingBox(1, c)]
    print(f"  free curve {c}: length {gmsh.model.occ.getMass(1, c) * S:.4f} m, y {bb[1]:.4f}..{bb[4]:.4f}")
if free or multi:
    sys.exit("ABORT: the mirrored wing is not closed (free or multiply-shared edges)")
gmsh.write(prefix + ".brep")

edge_curves = []
for c in crv:
    bb = [v * S for v in gmsh.model.getBoundingBox(1, c)]
    if (bb[4] - bb[1]) > 2.0 or bb[1] > 29.0 or bb[4] < -29.0:
        edge_curves.append(c)
print(f"size field: {len(edge_curves)} edge curves (spanwise leading/trailing edges, tip) of {len(crv)} curves")
fd = gmsh.model.mesh.field.add("Distance")
gmsh.model.mesh.field.setNumbers(fd, "CurvesList", edge_curves)
gmsh.model.mesh.field.setNumber(fd, "Sampling", 400)
ft = gmsh.model.mesh.field.add("Threshold")
gmsh.model.mesh.field.setNumber(ft, "InField", fd)
gmsh.model.mesh.field.setNumber(ft, "SizeMin", EDGE_MIN / S)
gmsh.model.mesh.field.setNumber(ft, "SizeMax", SIZE_MAX / S)
gmsh.model.mesh.field.setNumber(ft, "DistMin", D_MIN / S)
gmsh.model.mesh.field.setNumber(ft, "DistMax", D_MAX / S)
gmsh.model.mesh.field.setAsBackgroundMesh(ft)
gmsh.option.setNumber("Mesh.MeshSizeMax", SIZE_MAX / S)
gmsh.option.setNumber("Mesh.MeshSizeMin", SIZE_FLOOR / S)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", CURV)
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.Algorithm", 6)
gmsh.model.mesh.generate(2)


def triangles():
    nt, co, _ = gmsh.model.mesh.getNodes()
    nt = np.asarray(nt, dtype=np.int64)
    xyz = np.zeros((nt.max() + 1, 3))
    xyz[nt] = np.asarray(co).reshape(-1, 3)
    parts = []
    for ty, nn in zip(*gmsh.model.mesh.getElements(2)[::2]):
        if ty == 2:
            parts.append(xyz[np.asarray(nn, dtype=np.int64).reshape(-1, 3)])
    return len(nt), np.concatenate(parts)


def stats(tri):
    a = 0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1).sum() * S * S
    v = np.einsum("ij,ij->i", tri[:, 0], np.cross(tri[:, 1], tri[:, 2])).sum() / 6.0 * S ** 3
    return a, v


nn, tri = triangles()
marea, vol = stats(tri)
print(f"surface mesh: {nn} nodes, {len(tri)} triangles, mesh area {marea:.3f} m2 (CAD kernel {area:.3f} m2), signed volume {vol:.3f} m3")
if vol < 0:
    gmsh.model.mesh.reverse()
    nn, tri = triangles()
    marea, vol = stats(tri)
    print(f"element orientation reversed: signed volume now {vol:.3f} m3")
if vol <= 0:
    sys.exit("ABORT: signed volume is not positive after the orientation fix")

gmsh.option.setNumber("Mesh.ScalingFactor", S)
gmsh.option.setNumber("Mesh.Binary", 1)
gmsh.write(prefix + ".stl")
print(f"wrote {prefix}.stl ({os.path.getsize(prefix + '.stl') / 1e6:.1f} MB, metres), {prefix}.brep and {prefix}_unsewn.brep (model units)")
gmsh.finalize()
