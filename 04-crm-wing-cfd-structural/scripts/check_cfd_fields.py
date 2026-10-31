#!/usr/bin/env python3
"""Field-extreme diagnostics for the Euler wing run.

Usage: check_cfd_fields.py [time] [case_dir]

Reads VTK/euler_<time>.vtk and the patch VTK files written by foamToVTK, and
reports extremes (with locations) of Ma, p, T, rho and |U|, percentiles, the
state of the far-field patches relative to freestream, and Cp on the wing.
"""
import glob
import sys
import numpy as np
import pyvista as pv

time = sys.argv[1] if len(sys.argv) > 1 else "600"
case = sys.argv[2] if len(sys.argv) > 2 else "simulation/cfd_v1/euler"

P_INF = 21662.7     # Pa
Q_INF = 10955.9     # Pa
T_INF = 216.65      # K
MA_INF = 0.85

mesh = pv.read(sorted(glob.glob(f"{case}/VTK/*_{time}.vtk"))[0])
cd = mesh.cell_data
ctr = mesh.cell_centers().points
print(f"internal VTK: {mesh.n_cells} cells (polyhedra are decomposed), "
      f"{mesh.n_points} points")
print("cell arrays:", list(cd.keys()))

fields = {
    "Ma": np.asarray(cd["Ma"]),
    "p [Pa]": np.asarray(cd["p"]),
    "T [K]": np.asarray(cd["T"]),
    "rho": np.asarray(cd["rho"]),
    "|U| [m/s]": np.linalg.norm(np.asarray(cd["U"]), axis=1),
}

print("\n=== extremes with locations (x y z in m) ===")
for name, a in fields.items():
    i0, i1 = int(np.argmin(a)), int(np.argmax(a))
    print(f"{name:11s} min {a[i0]:12.5g} at ({ctr[i0,0]:8.2f} {ctr[i0,1]:8.2f} {ctr[i0,2]:8.2f})"
          f" | max {a[i1]:12.5g} at ({ctr[i1,0]:8.2f} {ctr[i1,1]:8.2f} {ctr[i1,2]:8.2f})")

print("\n=== percentiles ===")
pcts = [0.01, 1, 50, 99, 99.99]
print(f"{'':11s} " + " ".join(f"{p:>11g}%" for p in pcts))
for name in ("Ma", "p [Pa]", "T [K]"):
    v = np.percentile(fields[name], pcts)
    print(f"{name:11s} " + " ".join(f"{x:12.5g}" for x in v))
ma = fields["Ma"]
print(f"fraction of cells with Ma > 1.0: {np.mean(ma > 1.0):.4f}, "
      f"> 1.5: {np.mean(ma > 1.5):.5f}, > 2.0: {np.mean(ma > 2.0):.6f}")

print("\n=== far-field patches (p/pInf, Ma, T/Tinf) ===")
for patch in ("inlet", "outlet", "top", "bottom", "farfieldY", "symmetry"):
    g = pv.read(f"{case}/VTK/{patch}/{patch}_{time}.vtk")
    p = np.asarray(g.cell_data["p"]); t = np.asarray(g.cell_data["T"])
    m = np.asarray(g.cell_data["Ma"])
    print(f"{patch:10s} n={g.n_cells:6d} p/pInf [{p.min()/P_INF:.4f} {p.mean()/P_INF:.4f} {p.max()/P_INF:.4f}]"
          f"  Ma [{m.min():.4f} {m.mean():.4f} {m.max():.4f}]"
          f"  T/Tinf [{t.min()/T_INF:.4f} {t.mean()/T_INF:.4f} {t.max()/T_INF:.4f}]")

print("\n=== wing surface ===")
w = pv.read(f"{case}/VTK/wing/wing_{time}.vtk")
pw = np.asarray(w.cell_data["p"])
cp = (pw - P_INF) / Q_INF
wc = w.cell_centers().points
j0, j1 = int(np.argmin(cp)), int(np.argmax(cp))
print(f"wing faces: {w.n_cells}")
print(f"Cp min {cp[j0]:.4f} at ({wc[j0,0]:.2f} {wc[j0,1]:.2f} {wc[j0,2]:.2f})")
print(f"Cp max {cp[j1]:.4f} at ({wc[j1,0]:.2f} {wc[j1,1]:.2f} {wc[j1,2]:.2f})")
print(f"Cp percentiles 1/50/99: {np.percentile(cp,[1,50,99])}")
mw = np.asarray(w.cell_data["Ma"])
k = int(np.argmax(mw))
print(f"wing-surface Ma max {mw[k]:.4f} at ({wc[k,0]:.2f} {wc[k,1]:.2f} {wc[k,2]:.2f})")
