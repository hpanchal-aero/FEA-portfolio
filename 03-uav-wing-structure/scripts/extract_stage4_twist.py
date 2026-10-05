#!/usr/bin/env python3
"""
Stage 4: nose-up rotation (twist) of the rib planes from a static .frd.
Run with the system python (numpy only):
  /usr/bin/python3 scripts/extract_stage4_twist.py <job1> [<job2> ...]
Each <job> needs <job>.frd and <job>.inp (static run, exactly one DISP block).

At each rib plane (solid original nodes with |z - z_rib| <= 0.76 mm; root rib uses z >= 0),
fit Uy = a + b*x over the chord by least squares. theta = -b is the rotation about the
spanwise axis, positive = nose-up (leading edge, x < 0, moves +y relative to the trailing edge).
The ribs are not rigid, so the rms fit residual is reported; the rotation is a fitted
section rotation, not a rigid-body rotation.
"""
import sys
import os
import numpy as np

RIB_Z = [0.0, 200.0, 400.0, 600.0, 800.0, 1000.0]
HALF = 0.76


def n_orig_nodes(inp):
    on = False
    mx = 0
    n = 0
    with open(inp) as f:
        for line in f:
            if not on:
                if line.strip().upper() == "*NODE":
                    on = True
                continue
            if line.startswith("*") and not line.startswith("**"):
                break
            if line.strip() and not line.startswith("**"):
                n += 1
                mx = max(mx, int(line.split(",")[0]))
    return n, mx


def read_disp(path):
    ids, vals = [], []
    in_block = False
    nblk = 0
    with open(path) as f:
        for line in f:
            if line.startswith(" -4"):
                in_block = (line[5:13].strip() == "DISP")
                if in_block:
                    nblk += 1
                continue
            if in_block:
                if line.startswith(" -3"):
                    in_block = False
                elif line.startswith(" -1"):
                    ids.append(int(line[3:13]))
                    vals.append([float(line[13 + 12 * k:25 + 12 * k]) for k in range(3)])
    if not ids:
        print(f"ERROR: DISP block not found in {path}")
        sys.exit(1)
    ids = np.array(ids)
    arr = np.zeros((ids.max() + 1, 3))
    arr[ids] = np.array(vals)
    return arr, nblk


def read_frd_nodes(path):
    ids, xyz = [], []
    on = False
    with open(path) as f:
        for line in f:
            if line.startswith("    2C"):
                on = True
                continue
            if on:
                if line.startswith(" -3"):
                    break
                if line.startswith(" -1"):
                    ids.append(int(line[3:13]))
                    xyz.append([float(line[13 + 12 * k:25 + 12 * k]) for k in range(3)])
    return np.array(ids), np.array(xyz)


def analyse(job):
    frd, inp = job + ".frd", job + ".inp"
    for p in (frd, inp):
        if not os.path.exists(p):
            print(f"ERROR: missing {p}")
            sys.exit(1)
    U, nblk = read_disp(frd)
    if nblk != 1:
        print(f"ERROR: {frd} has {nblk} DISP blocks, expected exactly 1 (static run)")
        sys.exit(1)
    ids, xyz = read_frd_nodes(frd)
    n0, mx = n_orig_nodes(inp)
    if n0 != mx:
        print(f"WARNING: .inp node count {n0} != max id {mx}; using max id as N0")
    orig = ids <= mx
    print(f"\n=== {job}")
    print(f"frd nodes {len(ids)}, original IDs present {int(orig.sum())}, N0 = {mx}")
    print(" z_rib   n_nodes  x-range [mm]        Uy(x=0) [mm]   theta [rad]   theta [deg]   rms resid [mm]")
    out = []
    for zc in RIB_Z:
        lo, hi = (0.0, HALF) if zc == 0.0 else (zc - HALF, zc + HALF)
        m = orig & (xyz[:, 2] >= lo - 1e-9) & (xyz[:, 2] <= hi + 1e-9)
        if m.sum() < 10:
            print(f"{zc:6.0f}   too few nodes ({int(m.sum())})")
            out.append(np.nan)
            continue
        x = xyz[m, 0]
        uy = U[ids[m], 1]
        M = np.vstack([np.ones_like(x), x]).T
        coef, *_ = np.linalg.lstsq(M, uy, rcond=None)
        rms = float(np.sqrt(np.mean((uy - M @ coef) ** 2)))
        th = -coef[1]
        out.append(th)
        print(f"{zc:6.0f}   {int(m.sum()):6d}   [{x.min():7.2f},{x.max():7.2f}]   {coef[0]:11.6f}   "
              f"{th:11.3e}   {np.degrees(th):11.5f}   {rms:12.6f}")
    return np.array(out)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    jobs = sys.argv[1:]
    res = [analyse(j) for j in jobs]
    if len(jobs) > 1:
        print("\n=== theta [rad], nose-up positive, by rib station")
        print(" z_rib  " + "  ".join(f"{os.path.basename(j)[-26:]:>26s}" for j in jobs))
        for i, zc in enumerate(RIB_Z):
            print(f"{zc:6.0f}  " + "  ".join(f"{r[i]:26.3e}" for r in res))


if __name__ == "__main__":
    main()
