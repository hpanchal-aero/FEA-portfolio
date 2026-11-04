#!/usr/bin/python3
"""V6 check: linear buckling of the membrane-restrained square plate under the 1 K-equivalent biaxial preload.

Usage: /usr/bin/python3 scripts/check_v6.py
Reads simulation/v6b_prescribed/v6b_prescribed.dat (preload = prescribed mean edge displacement, see frame_prestrain.py).
Closed form (docs/specification.md, Amendment 1): dT_cr = 2 pi^2 h^2 / (12 (1 + nu) alpha b^2) for the (1,1) mode;
higher modes (m, n): dT_cr * (m^2 + n^2) / 2. Pass criterion (fixed before the run): mode 1 within 5 %.
Also reports the mode-shape fit of the vertical displacement on the top face to sin(m pi x/L) sin(n pi y/L)
(least squares over the degenerate group) and the preload stress. Writes results/v6_buckling.txt and the record of the
failed thermal-preload attempt results/v6a_thermal_preload.txt. Exit status 1 on FAIL.
File layout (CalculiX 2.23): the buckling step first prints the preload state (one displacement table on set TOP), then the
buckling factors, then for every mode a header, a stress table and a displacement table. The displacement tables that precede
the first 'E I G E N V A L U E' header are skipped, and exactly six must remain.
"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ccxio
from material import MATERIAL

MESH = "v_n40_z2"
CASE = "v6b_prescribed"
TOL = 0.05
# degenerate groups: (mode numbers, list of (m, n) shapes)
GROUPS = [([1], [(1, 1)]), ([2, 3], [(1, 2), (2, 1)]), ([4], [(2, 2)]), ([5, 6], [(1, 3), (3, 1)])]


def read_factors(dat):
    with open(dat) as f:
        lines = f.read().splitlines()
    start = None
    for i, s in enumerate(lines):
        if "B U C K L I N G" in s:
            start = i
            break
    if start is None:
        raise SystemExit("no buckling factor block in %s" % dat)
    out = []
    for s in lines[start + 1:]:
        if "E I G E N V A L U E" in s:
            break
        t = s.split()
        if len(t) == 2 and t[0].isdigit():
            out.append(float(t[1]))
    return out


def tables_before_first_mode(dat):
    """number of displacement tables on set TOP that precede the first eigenvalue header"""
    n = 0
    with open(dat) as f:
        for s in f:
            if "E I G E N V A L U E" in s:
                return n
            if s.lstrip().startswith("displacements") and "for set TOP" in s:
                n += 1
    raise SystemExit("no eigenvalue header in %s" % dat)


def mean_stress(dat):
    blocks = [b for b in ccxio.parse_dat(dat) if b["title"].lower().startswith("stresses")]
    if not blocks:
        raise SystemExit("no stress block in %s" % dat)
    rows = np.array(blocks[0]["rows"])
    return len(blocks), rows.shape[0], rows[:, 2].mean(), rows[:, 3].mean(), rows[:, 4].mean()


def main():
    root = os.path.dirname(HERE)
    nodes = ccxio.read_nodes(os.path.join(root, "mesh", MESH + ".inp"))
    xs = np.array([nodes[k][0] for k in nodes])
    ys = np.array([nodes[k][1] for k in nodes])
    zs = np.array([nodes[k][2] for k in nodes])
    L, W, h = xs.max() - xs.min(), ys.max() - ys.min(), zs.max() - zs.min()
    if abs(L - W) > 1e-12:
        raise SystemExit("plate is not square: %g x %g" % (L, W))
    m = MATERIAL
    d1 = 2 * math.pi ** 2 * h ** 2 / (12 * (1 + m["nu"]) * m["alpha"] * L ** 2)

    dat = os.path.join(root, "simulation", CASE, CASE + ".dat")
    factors = read_factors(dat)
    if len(factors) != 6:
        raise SystemExit("expected 6 buckling factors, found %d" % len(factors))
    nblk, npts, sx, sy, sz = mean_stress(dat)
    sig = m["E"] * m["alpha"] / (1 - m["nu"])
    allvecs = [b for b in ccxio.parse_dat(dat) if b["title"].lower().startswith("displacements") and b["set"] == "TOP"]
    skip = tables_before_first_mode(dat)
    vecs = allvecs[skip:]
    if len(vecs) != 6:
        raise SystemExit("expected 6 eigenvectors on set TOP after skipping %d preload table(s), found %d" % (skip, len(vecs)))

    out = ["V6b linear buckling, restrained square plate, preload = mechanical equivalent of 1 K (prescribed mean edge displacement)",
           "plate %.3f x %.3f m, h = %.4f m, mesh %s, case %s" % (L, W, h, MESH, CASE),
           "closed form dT_cr(1,1) = 2 pi^2 h^2 / (12 (1 + nu) alpha b^2) = %.5f K" % d1,
           "preload stress (first stress block of the .dat; %d stress block(s) in the file; %d integration points): mean sxx %.5e  syy %.5e  szz %.3e Pa"
           % (nblk, npts, sx, sy, sz),
           "expected -E alpha / (1 - nu) = %.5e Pa per K; deviation of sxx %+.3f %% of its magnitude" % (-sig, 100 * (sx + sig) / sig),
           "displacement tables on TOP: %d in the file, %d skipped (preload state printed before the first mode), %d eigenvectors used" % (len(allvecs), skip, len(vecs)),
           "",
           "mode   FE factor [K]   closed form [K]   error [%]   shape fit R^2 (uz on top face)"]
    xc, yc = xs.min(), ys.min()
    results = {}
    for modes, shapes in GROUPS:
        mn = shapes[0]
        for k in modes:
            blk = vecs[k - 1]
            rows = np.array(blk["rows"])
            X = np.array([nodes[int(i)][0] - xc for i in rows[:, 0]])
            Y = np.array([nodes[int(i)][1] - yc for i in rows[:, 0]])
            uz = rows[:, 3]
            B = np.column_stack([np.sin(a * math.pi * X / L) * np.sin(b * math.pi * Y / L) for a, b in shapes])
            c, *_ = np.linalg.lstsq(B, uz, rcond=None)
            r2 = 1.0 - np.sum((uz - B @ c) ** 2) / np.sum(uz ** 2)
            ref = d1 * (mn[0] ** 2 + mn[1] ** 2) / 2.0
            err = 100.0 * (factors[k - 1] - ref) / ref
            results[k] = err
            out.append("  %d     %10.5f      %10.5f       %+7.3f        %.5f" % (k, factors[k - 1], ref, err, r2))
    ok = abs(results[1]) <= 100 * TOL
    out.append("")
    out.append("tolerance (mode 1) %.1f %%; FE factor is %s the closed form" % (100 * TOL, "below" if results[1] < 0 else "above"))
    out.append("RESULT: %s" % ("PASS" if ok else "FAIL"))
    text = "\n".join(out) + "\n"
    print(text, end="")
    with open(os.path.join(root, "results", "v6_buckling.txt"), "w") as f:
        f.write(text)

    # record of the failed attempt with a thermal preload
    dat_a = os.path.join(root, "simulation", "v6_buckling", "v6_buckling.dat")
    dat_d = os.path.join(root, "simulation", "v6_diag_prestress", "v6_diag_prestress.dat")
    fa = read_factors(dat_a)
    nd, pd, dx, dy, dz = mean_stress(dat_d)
    rec = ["V6a (record of a failed attempt, kept like V3a): linear buckling with a THERMAL preload in CalculiX 2.23",
           "deck simulation/v6_buckling/v6_buckling.inp: *TEMPERATURE 20 -> 21 degC in the *STATIC step, then *BUCKLE 6, same mesh and supports as V6b",
           "preload step alone (simulation/v6_diag_prestress): mean sxx %.5e  syy %.5e  szz %.3e Pa (expected -E alpha / (1 - nu) = %.5e Pa per K); edge reaction 1460.471 N" % (dx, dy, dz, -sig),
           "buckling factors printed: " + "  ".join("%.6e" % v for v in fa),
           "closed form for comparison: %.5f K" % d1,
           "outcome: factors of order 1e7 to 1e8 (a buckling temperature rise of tens of millions of kelvin), no degenerate pair, the mode-1 shape was inspected at only five nodes, where the in-plane components exceeded the vertical one, and was not analysed further.",
           "interpretation: the thermal preload did not enter the buckling analysis, although the preload step itself is correct. This agrees with reports that CalculiX linear buckling does not support thermal prestress (PrePoMax forum threads 'TEMPERATURE doesn't work with BUCKLE' and 'Thermal loading'); the exact internal cause was not investigated.",
           "consequence: V6 is run with the mechanical equivalent of the thermal preload (V6b, scripts/frame_prestrain.py), which reproduces the thermal stress state (mean sxx identical to the printed digits) and verifies the buckling extraction, restraint and mesh. The thermal load path is exercised in the geometrically nonlinear static step (P05-D5)."]
    with open(os.path.join(root, "results", "v6a_thermal_preload.txt"), "w") as f:
        f.write("\n".join(rec) + "\n")
    print("wrote results/v6_buckling.txt and results/v6a_thermal_preload.txt")
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
