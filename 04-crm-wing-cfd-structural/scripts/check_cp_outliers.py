#!/usr/bin/python3
"""Find isolated Cp outliers (single-face anomalies) on a CFD wing surface.

usage: check_cp_outliers.py "label1=wing1.vtk" ["label2=wing2.vtk" ...]
Cp = (p - 21662.7) / 10955.9.  For every face the deviation  dev = Cp - Cp_smooth  is computed, where
Cp_smooth is the cell->point->cell average of Cp over the faces sharing vertices.  An isolated anomaly of
size d appears with dev of about 0.75 d.  Prints counts above 0.05 / 0.1 / 0.2, the 12 largest outliers,
their distribution over the span and their face areas compared with the median face area.
"""
import sys

import numpy as np
import pyvista as pv

PINF, Q = 21662.7, 10955.9

for arg in sys.argv[1:]:
    label, path = arg.split("=", 1)
    m = pv.read(path)
    if not isinstance(m, pv.PolyData):
        m = m.extract_surface()
    cp = (np.asarray(m.cell_data["p"]) - PINF) / Q
    m.cell_data["Cp"] = cp
    sm = np.asarray(m.cell_data_to_point_data().point_data_to_cell_data().cell_data["Cp"])
    dev = cp - sm
    C = np.asarray(m.cell_centers().points)
    A = np.asarray(m.compute_cell_sizes(length=False, area=True, volume=False).cell_data["Area"])
    print(f"\n=== {label}: faces {m.n_cells}, median face area {np.median(A):.3e} m2, Cp range {cp.min():.3f}..{cp.max():.3f}")
    for thr in (0.05, 0.1, 0.2):
        k = np.abs(dev) > thr
        print(f"|dev| > {thr}: {int(k.sum())} faces ({100 * k.mean():.3f}%)")
    k = np.abs(dev) > 0.1
    if k.any():
        yb = np.arange(0, 35, 5)
        h, _ = np.histogram(C[k, 1], bins=yb)
        print("faces |dev|>0.1 by span (5 m bins from y=0):", dict(zip(yb[:-1].tolist(), h.tolist())))
        print(f"median area of those faces {np.median(A[k]):.3e} m2 (all faces {np.median(A):.3e})")
    print("12 largest |dev| (x, y, z, Cp, Cp_smooth, dev, area):")
    for i in np.argsort(-np.abs(dev))[:12]:
        print(f"  {C[i, 0]:7.2f} {C[i, 1]:6.2f} {C[i, 2]:6.2f}  {cp[i]:+.3f}  {sm[i]:+.3f}  {dev[i]:+.3f}  {A[i]:.2e}")
