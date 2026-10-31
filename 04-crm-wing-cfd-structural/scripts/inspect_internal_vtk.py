#!/usr/bin/python3
"""Inspect an internal-mesh VTK written by foamToVTK: arrays, cell types, Mach statistics.

usage: inspect_internal_vtk.py <file.vtk>
Uses the solver's own Ma cell array.  foamToVTK splits polyhedral cells into several pieces, so
statistics are volume-weighted, and OpenFOAM cells are counted through the cellID array.
The free-stream check reads the inlet patch file of the same time next to the internal file.
"""
import os
import sys
import time

import numpy as np
import pyvista as pv

path = sys.argv[1]
t0 = time.time()
m = pv.read(path)
print(f"read in {time.time() - t0:.1f} s: {type(m).__name__}, pieces {m.n_cells}, points {m.n_points}")
print("bounds", [round(float(b), 2) for b in m.bounds])
u, c = np.unique(m.celltypes, return_counts=True)
print("piece types (VTK id: count):", {int(a): int(b) for a, b in zip(u, c)})
print("cell arrays:", list(m.cell_data.keys()))

M = np.asarray(m.cell_data["Ma"])
cid = np.asarray(m.cell_data["cellID"]).astype(np.int64)
P = np.asarray(m.cell_centers().points)
t1 = time.time()
V = np.abs(np.asarray(m.compute_cell_sizes(length=False, area=False, volume=True).cell_data["Volume"]))
print(f"volumes in {time.time() - t1:.1f} s: total {V.sum():.4g} m3, OpenFOAM cells {np.unique(cid).size}")

order = np.argsort(M)
cw = np.cumsum(V[order])
pc = [float(M[order][np.searchsorted(cw, q / 100.0 * cw[-1])]) for q in (1, 50, 99)]
print("volume-weighted Mach percentiles 1/50/99:", [round(x, 4) for x in pc], " min", round(float(M.min()), 4), " max", round(float(M.max()), 4))
i = int(np.argmax(M))
print(f"max Mach {M[i]:.4f} at piece centre x {P[i, 0]:.2f}  y {P[i, 1]:.2f}  z {P[i, 2]:.2f}  (cellID {cid[i]})")

for thr in (1.0, 1.2, 1.4):
    k = M > thr
    if not k.any():
        print(f"Mach > {thr}: none")
        continue
    nc = np.unique(cid[k]).size
    print(f"Mach > {thr}: OpenFOAM cells {nc}, volume {V[k].sum():.4g} m3 ({100 * V[k].sum() / V.sum():.5f}% of domain), "
          f"box x {P[k, 0].min():.2f}..{P[k, 0].max():.2f}  y {P[k, 1].min():.2f}..{P[k, 1].max():.2f}  z {P[k, 2].min():.2f}..{P[k, 2].max():.2f}")

base = os.path.dirname(path)
t = os.path.basename(path)[:-4].rsplit("_", 1)[1]
inl = os.path.join(base, "inlet", f"inlet_{t}.vtk")
if os.path.exists(inl):
    q = pv.read(inl)
    print("inlet patch arrays:", list(q.cell_data.keys()))
    if "Ma" in q.cell_data:
        a = np.asarray(q.cell_data["Ma"])
        print(f"inlet Ma: min {a.min():.4f}  mean {a.mean():.4f}  max {a.max():.4f}")
else:
    print("no inlet patch file at", inl)
