#!/usr/bin/python3
"""Per-section rigid kinematics of a CalculiX static result, with reconstruction check.

usage: section_kinematics.py <job.frd> <out.csv> [bin_width_m = 1.0] [wing.stl]

Per spanwise bin: uz = a + b (x - xbar) + c (y - ybar) (least squares), mean ux, uy.
twist = -b (negative = nose-down washout).  A root point (y = 0, zero motion) is prepended.
Reconstruction: u(x,y) = [ux(y), uy(y), a(y) + b(y) (x - xbar(y))] with linear interpolation in y.
Nodes outside the wing envelope (x < 20 m or z < 2 m) are listed and excluded from the fits.
"""
import csv
import sys

import numpy as np

frd, out = sys.argv[1], sys.argv[2]
width = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
stl = sys.argv[4] if len(sys.argv) > 4 else None
XMIN, ZMIN = 20.0, 2.0


def read_frd(path):
    coords, disp, block = {}, {}, None
    with open(path) as f:
        for line in f:
            if line.startswith("    2C"):
                block = "coord"
            elif line.startswith(" -4  DISP"):
                block = "disp"
            elif line.startswith(" -4  "):
                block = None
            elif line.startswith(" -3"):
                block = None
            elif line.startswith(" -1") and block:
                node = int(line[3:13])
                vals = [float(line[13 + 12 * i:25 + 12 * i]) for i in range(3)]
                (coords if block == "coord" else disp)[node] = vals
    return coords, disp


coords, disp = read_frd(frd)
ids = np.array(sorted(coords))
missing = [i for i in ids if i not in disp]
print(f"frd nodes {len(coords)}, with displacement {len(disp)}, without {len(missing)}")
if missing:
    sys.exit("ABORT: nodes without displacement")
X = np.array([coords[i] for i in ids])
U = np.array([disp[i] for i in ids])
odd = (X[:, 0] < XMIN) | (X[:, 2] < ZMIN)
print(f"nodes outside wing envelope (x<{XMIN} or z<{ZMIN}): {odd.sum()}")
for k in np.where(odd)[0][:6]:
    print(f"  node {ids[k]}  xyz {X[k]}  u {U[k]}")
W = ~odd
Xw, Uw = X[W], U[W]

ymax = Xw[:, 1].max()
nb = int(np.ceil(ymax / width))
rows = []
for k in range(nb):
    m = (Xw[:, 1] >= k * width) & (Xw[:, 1] < (k + 1) * width + (1e-9 if k == nb - 1 else 0))
    if m.sum() < 6:
        continue
    xb, yb = Xw[m, 0].mean(), Xw[m, 1].mean()
    A = np.column_stack([np.ones(m.sum()), Xw[m, 0] - xb, Xw[m, 1] - yb])
    (a, b, c), *_ = np.linalg.lstsq(A, Uw[m, 2], rcond=None)
    rows.append((yb, m.sum(), xb, Uw[m, 0].mean(), Uw[m, 1].mean(), a, b, c))
R = np.array(rows)
yc, n, xbar, ux, uy, a, b, c = R.T
twist = -np.degrees(np.arctan(b))

print("\n  y_c[m]    n   ux[m]    uy[m]    uz[m]   twist[deg]  dz/dy")
for k in range(len(yc)):
    print(f"{yc[k]:8.3f} {int(n[k]):5d} {ux[k]:8.4f} {uy[k]:8.4f} {a[k]:8.4f} {twist[k]:9.3f} {c[k]:9.4f}")

with open(out, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["y_c", "n", "xbar", "ux", "uy", "uz", "b_slope", "c_slope", "twist_deg"])
    for k in range(len(yc)):
        w.writerow([f"{yc[k]:.6f}", int(n[k]), f"{xbar[k]:.6f}", f"{ux[k]:.8f}", f"{uy[k]:.8f}",
                    f"{a[k]:.8f}", f"{b[k]:.10f}", f"{c[k]:.10f}", f"{twist[k]:.6f}"])
print(f"\nwrote {out}")

# reconstruction with a zero-motion root point
Y = np.concatenate([[0.0], yc])
f_ux = np.concatenate([[0.0], ux]); f_uy = np.concatenate([[0.0], uy])
f_a = np.concatenate([[0.0], a]); f_b = np.concatenate([[0.0], b]); f_xb = np.concatenate([[xbar[0]], xbar])
y = Xw[:, 1]
mx = np.interp(y, Y, f_ux); my = np.interp(y, Y, f_uy)
mz = np.interp(y, Y, f_a) + np.interp(y, Y, f_b) * (Xw[:, 0] - np.interp(y, Y, f_xb))
tip = y >= ymax - 1.0
print("\nreconstruction of wing-node displacement from section kinematics")
for name, mod, ref in (("ux", mx, Uw[:, 0]), ("uy", my, Uw[:, 1]), ("uz", mz, Uw[:, 2])):
    e = mod - ref
    print(f"  {name}: rms err {np.sqrt((e**2).mean()):.4f} m   max |err| {np.abs(e).max():.4f} m   "
          f"(ref max |u| {np.abs(ref).max():.4f} m)")
print(f"  tip bin (y >= {ymax - 1.0:.2f}) mean uz: FEA {Uw[tip, 2].mean():.4f}   model {mz[tip].mean():.4f}")

if stl:
    rec = np.dtype([("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
    with open(stl, "rb") as f:
        f.read(80)
        nt = int(np.frombuffer(f.read(4), dtype="<u4")[0])
        t = np.fromfile(f, dtype=rec, count=nt)
    V = t["v"].reshape(-1, 3).astype(np.float64)
    V = V[V[:, 1] >= 0]
    print("\nframe check, deck wing nodes vs STL vertices (y >= 0), every 3rd bin")
    print("  y_bin     deck x range         STL x range          deck z range     STL z range")
    for k in range(0, nb, 3):
        lo, hi = k * width, (k + 1) * width
        md = (Xw[:, 1] >= lo) & (Xw[:, 1] < hi)
        ms = (V[:, 1] >= lo) & (V[:, 1] < hi)
        if md.sum() and ms.sum():
            print(f"  {lo:4.0f}-{hi:<4.0f} {Xw[md,0].min():7.3f}..{Xw[md,0].max():7.3f}   "
                  f"{V[ms,0].min():7.3f}..{V[ms,0].max():7.3f}   "
                  f"{Xw[md,2].min():6.3f}..{Xw[md,2].max():6.3f}  {V[ms,2].min():6.3f}..{V[ms,2].max():6.3f}")
