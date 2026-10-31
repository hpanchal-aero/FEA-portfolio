#!/usr/bin/env python3
"""Pressure and temperature extremes of CFD internal-field VTK files and the number of cells at the limiter bounds.

usage: check_limiter_fields.py <internal.vtk> [<internal.vtk> ...]     (run from the project root; writes results/limiter_fields.txt)
Bounds (system/fvConstraints): limitPressure 0.1 p_inf = 2166.27 Pa to 2 p_inf = 43325.4 Pa; limitTemperature 100 K to 400 K.
A cell counts as 'at a bound' when it is within 0.5 Pa (pressure) or 0.01 K (temperature) of it.
The VTK files hold a decomposition of the OpenFOAM mesh (about 3 VTK cells per OpenFOAM cell), so the counts are VTK cells.
"""
import gc
import os
import sys

import numpy as np
import pyvista as pv

PINF = 21662.7
PMIN, PMAX = 0.1 * PINF, 2.0 * PINF
TMIN, TMAX = 100.0, 400.0
OUT = "results/limiter_fields.txt"


def arr(m, name):
    if name in m.cell_data:
        return np.asarray(m.cell_data[name])
    return np.asarray(m.point_data[name])


if len(sys.argv) < 2:
    sys.exit(__doc__)

lines = ["Limiter bounds: p 2166.27 .. 43325.4 Pa, T 100 .. 400 K; 'at bound' = within 0.5 Pa / 0.01 K (VTK cells)"]
for path in sys.argv[1:]:
    m = pv.read(path)
    p = arr(m, "p")
    T = arr(m, "T")
    label = os.path.relpath(path, "simulation/cfd_v1/runs")
    line = (f"{label:46s} cells {len(p):8d}  p {p.min():9.1f} .. {p.max():9.1f}  T {T.min():7.2f} .. {T.max():7.2f}"
            f"  at p floor {int((p <= PMIN + 0.5).sum()):4d}  p ceiling {int((p >= PMAX - 0.5).sum()):4d}"
            f"  T floor {int((T <= TMIN + 0.01).sum()):4d}  T ceiling {int((T >= TMAX - 0.01).sum()):4d}")
    print(line, flush=True)
    lines.append(line)
    del m, p, T
    gc.collect()

os.makedirs("results", exist_ok=True)
with open(OUT, "w") as f:
    f.write("\n".join(lines) + "\n")
print(f"written {OUT}")
