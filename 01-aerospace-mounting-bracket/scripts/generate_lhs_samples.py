"""
Stage 6 -- Latin Hypercube sample generation for the parametric sweep.

Generates 25 design points across the three free variables, using
scipy's Latin Hypercube sampler with a FIXED, documented seed for
exact reproducibility. The resulting points are written to a CSV that
is committed to the repository -- this CSV, not the seed alone, is
the authoritative record of which 25 designs were actually studied,
independent of any future change to scipy's LHS implementation or
this script.

Ranges (locked per Stage 6 sign-off):
    t         (flange thickness):     3.0  to  5.0  mm
    L_gusset  (nominal gusset leg):  15.0  to 25.0  mm
    r_toe     (toe fillet radius):    1.0  to  3.0  mm

Fixed (not sampled): R3 = 3.0mm, hole geometry.

SEED = 42 (arbitrary, fixed, documented here for reproducibility --
re-running this script with this seed reproduces the exact same 25
points, independent of the committed CSV).
"""

import numpy as np
from scipy.stats import qmc
import csv
import os

SEED = 42
N_SAMPLES = 25

T_RANGE = (3.0, 5.0)
L_GUSSET_RANGE = (15.0, 25.0)
R_TOE_RANGE = (1.0, 3.0)

OUT_CSV = "/results/stage6_parametric/lhs_design_points.csv"


def generate_samples():
    sampler = qmc.LatinHypercube(d=3, seed=SEED)
    unit_samples = sampler.random(n=N_SAMPLES)  # shape (25, 3), each column in [0,1)

    lower_bounds = [T_RANGE[0], L_GUSSET_RANGE[0], R_TOE_RANGE[0]]
    upper_bounds = [T_RANGE[1], L_GUSSET_RANGE[1], R_TOE_RANGE[1]]
    scaled = qmc.scale(unit_samples, lower_bounds, upper_bounds)

    return scaled


def main():
    samples = generate_samples()

    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    with open(OUT_CSV, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["design_id", "t", "L_gusset", "r_toe"])
        for i, (t, L_gusset, r_toe) in enumerate(samples, start=1):
            writer.writerow([f"lhs_{i:02d}", f"{t:.6f}", f"{L_gusset:.6f}", f"{r_toe:.6f}"])

    print(f"Generated {N_SAMPLES} LHS design points (seed={SEED}).")
    print(f"Saved to: {OUT_CSV}")
    print(f"\n{'design_id':<10}{'t (mm)':<12}{'L_gusset (mm)':<16}{'r_toe (mm)':<12}")
    print("-" * 50)
    for i, (t, L_gusset, r_toe) in enumerate(samples, start=1):
        print(f"lhs_{i:02d}    {t:<12.4f}{L_gusset:<16.4f}{r_toe:<12.4f}")


if __name__ == "__main__":
    main()
