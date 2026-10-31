#!/usr/bin/env python3
"""Force-coefficient history and convergence metrics for the Euler wing run.

Usage: cfd_force_history.py [forceCoeffs.dat] [window_iterations]

Reads an OpenFOAM 11 forceCoeffs.dat (columns: Time Cm Cd Cl Cl(f) Cl(r)),
prints the history every 50 iterations, and reports mean, range and drift of
Cl, Cd, Cm over the trailing window. A partially written last line is ignored.
"""
import sys
import numpy as np

path = sys.argv[1] if len(sys.argv) > 1 else (
    "simulation/cfd_v1/euler/postProcessing/forceCoeffsCompressible/0/forceCoeffs.dat")
win = int(sys.argv[2]) if len(sys.argv) > 2 else 100

rows = []
with open(path) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split()
        if len(parts) != 6:
            continue
        try:
            rows.append([float(x) for x in parts])
        except ValueError:
            continue

d = np.array(rows)
if d.shape[0] < 2:
    sys.exit("not enough rows in " + path)

it, cm, cd, cl = d[:, 0], d[:, 1], d[:, 2], d[:, 3]

print(f"file: {path}")
print(f"rows: {len(it)}  last iteration: {int(it[-1])}")
print(f"{'iter':>6} {'Cl':>10} {'Cd':>10} {'Cm':>10}")
for k in range(0, len(it), 50):
    print(f"{int(it[k]):6d} {cl[k]:10.5f} {cd[k]:10.5f} {cm[k]:10.5f}")
print(f"{int(it[-1]):6d} {cl[-1]:10.5f} {cd[-1]:10.5f} {cm[-1]:10.5f}  (last)")

n = min(win, len(it))
sel = slice(len(it) - n, None)
print(f"\ntrailing window: last {n} iterations")
for name, y in (("Cl", cl), ("Cd", cd), ("Cm", cm)):
    ys = y[sel]
    mean = ys.mean()
    rng = ys.max() - ys.min()
    slope = np.polyfit(it[sel], ys, 1)[0] * 100.0  # change per 100 iterations
    rel = 100.0 * rng / abs(mean) if mean != 0 else float("nan")
    print(f"{name}: mean {mean:10.5f}  range {rng:9.5f} ({rel:6.2f}% of mean)"
          f"  drift {slope:+.5f} per 100 it")
