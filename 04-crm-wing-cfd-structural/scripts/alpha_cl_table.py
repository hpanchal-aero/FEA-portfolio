#!/usr/bin/env python3
"""Force coefficients of the first-order upwind runs at selected iterations, and the linear estimate of the
angle of attack that gives CL = 0.5 from the two bracketing undeflected runs.

usage: alpha_cl_table.py            (writes results/alpha_cl_upwind.txt, run from the project root)
Reads simulation/cfd_v1/runs/<tag>/postProcessing/**/forceCoeffs.dat (columns: Time Cm Cd Cl Cl(f) Cl(r))
and the angle of attack from <tag>/0/U (the angleOfAttack entry).
"""
import glob
import os
import re

ROOT = "simulation/cfd_v1/runs"
TAGS = ["u_L8_p0p36", "u_p0p3", "u_p0p36", "u_p0p5", "u_L10_p0p36", "u_defl_p0p36", "u_defl_a2p07"]
ITS = (1148, 1300, 1500)
OUT = "results/alpha_cl_upwind.txt"


def alpha(tag):
    for ln in open(f"{ROOT}/{tag}/0/U"):
        m = re.match(r"\s*angleOfAttack\s+([-0-9.eE+]+)\s*;", ln)
        if m:
            return float(m.group(1))
    return float("nan")


def read_force(tag):
    rows = {}
    files = glob.glob(f"{ROOT}/{tag}/postProcessing/**/forceCoeffs.dat", recursive=True)
    for f in sorted(files, key=os.path.getmtime):
        for ln in open(f):
            if ln.lstrip().startswith("#") or not ln.strip():
                continue
            v = ln.split()
            try:
                rows[int(float(v[0]))] = (float(v[3]), float(v[2]), float(v[1]))  # Cl, Cd, Cm
            except (ValueError, IndexError):
                pass
    return rows


lines = []
lines.append("First-order upwind runs: Cl, Cd, Cm at selected iterations (forceCoeffs.dat; Aref 191.87 m2, lRef 7.01 m)")
lines.append(f"{'run':14s} {'alpha':>6s} {'last it':>7s} | " + " | ".join(f"it {i:>4d}: Cl      Cd      Cm" for i in ITS))
res = {}
for tag in TAGS:
    rows = read_force(tag)
    a = alpha(tag)
    last = max(rows) if rows else 0
    res[tag] = (a, rows)
    cells = []
    for i in ITS:
        if i in rows:
            cl, cd, cm = rows[i]
            cells.append(f"{cl:8.5f} {cd:7.5f} {cm:8.5f}")
        else:
            cells.append(f"{'-':>8s} {'-':>7s} {'-':>8s}")
    lines.append(f"{tag:14s} {a:6.2f} {last:7d} | " + " | ".join(cells))

a1, r1 = res["u_p0p3"]
a2, r2 = res["u_p0p5"]
if 1500 in r1 and 1500 in r2:
    c1, c2 = r1[1500][0], r2[1500][0]
    slope = (c2 - c1) / (a2 - a1)
    a_t = a1 + (0.5 - c1) * (a2 - a1) / (c2 - c1)
    lines.append("")
    lines.append(f"Undeflected wing, iteration 1500: CL({a1:.2f} deg) = {c1:.5f}, CL({a2:.2f} deg) = {c2:.5f}")
    lines.append(f"  lift-curve slope between them {slope:.4f} per degree; linear estimate of alpha for CL = 0.5: {a_t:.3f} deg")
else:
    lines.append("")
    lines.append("iteration 1500 missing in u_p0p3 or u_p0p5: no linear estimate")

os.makedirs("results", exist_ok=True)
with open(OUT, "w") as f:
    f.write("\n".join(lines) + "\n")
print("\n".join(lines))
print(f"written {OUT}")
