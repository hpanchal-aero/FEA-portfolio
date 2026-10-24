#!/usr/bin/env python3
"""Layer coverage on a wing patch from the mesh: first-cell thickness h = V_owner / A_face per patch face.

usage: check_layer_coverage.py <case_dir> <patch.vtk> [patch=wing]
needs <case>/0/Vc (postProcess -func writeCellVolumes), <case>/constant/polyMesh/{owner,boundary}
and the patch exported with foamToVTK (face order = patch face order).
"""
import re
import sys

import numpy as np
import pyvista as pv

case, vtk = sys.argv[1], sys.argv[2]
patch = sys.argv[3] if len(sys.argv) > 3 else "wing"


def read_list(path, integer):
    lines = open(path).read().split("\n")
    for i in range(len(lines) - 1):
        if lines[i].strip().isdigit() and lines[i + 1].strip() == "(":
            n = int(lines[i])
            if lines[i + 2 + n].strip() != ")":
                sys.exit(f"ABORT: unexpected list end in {path}")
            return np.array(lines[i + 2:i + 2 + n], dtype=np.int64 if integer else float)
    sys.exit(f"ABORT: no list found in {path}")


V = read_list(f"{case}/0/Vc", False)
owner = read_list(f"{case}/constant/polyMesh/owner", True)
txt = open(f"{case}/constant/polyMesh/boundary").read()
m = re.search(r"\n\s*" + patch + r"\s*\{(.*?)\}", txt, re.S)
if not m:
    sys.exit(f"ABORT: patch {patch} not found in boundary file")
nF = int(re.search(r"nFaces\s+(\d+);", m.group(1)).group(1))
sF = int(re.search(r"startFace\s+(\d+);", m.group(1)).group(1))
mesh = pv.read(vtk)
print(f"boundary file: patch {patch} nFaces {nF} startFace {sF}; VTK patch faces {mesh.n_cells}; cells {len(V)}; owner entries {len(owner)}")
if mesh.n_cells != nF:
    sys.exit("ABORT: VTK face count differs from the boundary file")
A = mesh.compute_sizes(length=False, area=True, volume=False).cell_data["Area"] if hasattr(mesh, "compute_sizes") else mesh.compute_cell_sizes(length=False, area=True, volume=False).cell_data["Area"]
h = V[owner[sF:sF + nF]] / A
cc = mesh.cell_centers().points
nz = mesh.compute_normals(cell_normals=True, point_normals=False).cell_data["Normals"][:, 2]
lay = h < 0.002
print(f"first-cell thickness h = V/A [mm]: median {1000*np.median(h):.3f}  10% {1000*np.percentile(h,10):.3f}  90% {1000*np.percentile(h,90):.3f}")
print(f"faces with h < 2 mm (layered): {int(lay.sum())} of {nF} ({100.0*lay.mean():.1f} %)")
for lo, hi in ((0, 0.0011), (0.0011, 0.002), (0.002, 0.005), (0.005, 0.015), (0.015, 1e9)):
    k = (h >= lo) & (h < hi)
    print(f"   h {1000*lo:6.1f}..{min(1000*hi, 9999):6.1f} mm: {int(k.sum()):7d} faces")
print("by face area (refinement level proxy):  faces, % layered")
for name, lo, hi in (("level 7 (~6e-3 m2)", 3e-3, 1e9), ("level 8 (~1.5e-3)", 7.5e-4, 3e-3), ("level 9 (~4e-4)", 1.9e-4, 7.5e-4), ("level 10+ (<1.9e-4)", 0, 1.9e-4)):
    k = (A >= lo) & (A < hi)
    if k.any():
        print(f"   {name:22s}: {int(k.sum()):7d}  {100.0*lay[k].mean():5.1f} %")
print("by y-range:")
for lo, hi in ((0, 3.05), (3.05, 10.9), (10.9, 20), (20, 28), (28, 29.6)):
    k = (cc[:, 1] >= lo) & (cc[:, 1] < hi)
    if k.any():
        print(f"   y {lo:5.2f}..{hi:5.2f}: {int(k.sum()):7d} faces  {100.0*lay[k].mean():5.1f} % layered")
print("by normal z-component (patch normals point out of the fluid, i.e. into the wing; nz<-0.3 is the upper surface):")
for name, k in (("nz < -0.3 (upper)", nz < -0.3), ("nz > +0.3 (lower)", nz > 0.3), ("|nz| <= 0.3 (LE/TE/tip)", np.abs(nz) <= 0.3)):
    if k.any():
        print(f"   {name:26s}: {int(k.sum()):7d} faces  {100.0*lay[k].mean():5.1f} % layered")
u = ~lay
if u.any():
    print(f"unlayered faces: {int(u.sum())}, centre x {cc[u,0].min():.2f}..{cc[u,0].max():.2f} (median {np.median(cc[u,0]):.2f}), y {cc[u,1].min():.2f}..{cc[u,1].max():.2f} (median {np.median(cc[u,1]):.2f}), z {cc[u,2].min():.2f}..{cc[u,2].max():.2f}")
