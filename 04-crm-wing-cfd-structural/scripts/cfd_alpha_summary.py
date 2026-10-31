#!/usr/bin/env python3
"""Summary of the two-stage alpha runs (first-order stage, then TVD stage).

Usage: cfd_alpha_summary.py <tag> [<tag> ...]
Cases are read from simulation/cfd_v1/runs/<tag>. Prints force coefficients at
selected iterations, trailing-window statistics, pressure-floor statistics per
stage from the logs, the last residuals, and (for two or more tags) a rough
linear estimate of the angle of attack giving CL = 0.5.
"""
import glob
import os
import re
import sys
import numpy as np

ROOT = "simulation/cfd_v1/runs"
TARGET_CL = 0.5
WINDOW = 100


def read_alpha(case):
    for line in open(f"{case}/0/U"):
        m = re.match(r"^angleOfAttack\s+(-?[\d.]+)\s*;", line)
        if m:
            return float(m.group(1))
    return float("nan")


def read_forces(case):
    pattern = f"{case}/postProcessing/forceCoeffsCompressible/*/forceCoeffs.dat"
    files = sorted(glob.glob(pattern),
                   key=lambda f: float(os.path.basename(os.path.dirname(f))))
    rows = {}
    for f in files:
        for line in open(f):
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            if len(p) != 6:
                continue
            try:
                v = [float(x) for x in p]
            except ValueError:
                continue
            rows[v[0]] = v
    return np.array([rows[k] for k in sorted(rows)])


def limiter_stats(log):
    if not os.path.exists(log):
        return None
    t = None
    nmin = nmax = 0
    lastmin = None
    lastval = None
    iters = 0
    for line in open(log):
        if line.startswith("Time ="):
            t = line.split()[2]
            iters += 1
        elif line.startswith("limitPressure: min"):
            nmin += 1
            lastmin = t
            lastval = float(line.split()[2])
        elif line.startswith("limitPressure: max"):
            nmax += 1
    return iters, nmin, nmax, lastmin, lastval


def last_residuals(case):
    files = sorted(glob.glob(f"{case}/postProcessing/residuals/*/residuals.dat"),
                   key=lambda f: float(os.path.basename(os.path.dirname(f))))
    if not files:
        return None
    last = None
    for line in open(files[-1]):
        if not line.startswith("#") and line.strip():
            last = line.split()
    return last


summary = []
for tag in sys.argv[1:]:
    case = f"{ROOT}/{tag}"
    alpha = read_alpha(case)
    d = read_forces(case)
    print(f"\n==================== {tag}  alpha = {alpha} deg")
    if d.shape[0] < 2:
        print("no force data yet")
        continue
    it, cm, cd, cl = d[:, 0], d[:, 1], d[:, 2], d[:, 3]
    print(f"rows {len(it)}, last iteration {int(it[-1])}")
    print(f"{'iter':>6} {'Cl':>9} {'Cd':>9} {'Cm':>9}")
    for target in range(100, int(it[-1]) + 1, 100):
        k = np.where(it == target)[0]
        if k.size:
            k = int(k[0])
            print(f"{target:6d} {cl[k]:9.4f} {cd[k]:9.4f} {cm[k]:9.4f}")
    n = min(WINDOW, len(it))
    sel = slice(len(it) - n, None)
    for name, y in (("Cl", cl), ("Cd", cd), ("Cm", cm)):
        slope = np.polyfit(it[sel], y[sel], 1)[0] * 100.0
        print(f"trailing {n} it: {name} mean {y[sel].mean():9.4f} "
              f"drift {slope:+.4f} per 100 it")
    for stage in ("log.stage1", "log.stage2"):
        s = limiter_stats(f"{case}/{stage}")
        if s is None:
            print(f"{stage}: not present")
        else:
            iters, nmin, nmax, lastmin, lastval = s
            print(f"{stage}: {iters} iterations; pressure floor limited in {nmin}, "
                  f"ceiling in {nmax}; last floor-limited iteration {lastmin}; "
                  f"its pre-limit min p {lastval}")
    r = last_residuals(case)
    if r:
        print("last residuals (Time p Ux Uy Uz h):", " ".join(r))
    summary.append((alpha, cl[sel].mean()))

if len(summary) >= 2:
    summary.sort()
    (a1, c1), (a2, c2) = summary[0], summary[-1]
    if c2 != c1:
        a_est = a1 + (TARGET_CL - c1) * (a2 - a1) / (c2 - c1)
        print(f"\nrough linear estimate: CL = {TARGET_CL} at alpha = {a_est:.2f} deg "
              f"(from trailing-window mean Cl at {a1} and {a2} deg; "
              f"indicative only, Cl is still drifting)")
