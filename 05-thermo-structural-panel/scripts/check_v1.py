#!/usr/bin/python3
"""V1: 1D steady conduction through the thickness, fixed face temperatures.

Usage: /usr/bin/python3 scripts/check_v1.py JOB.dat MESH.inp T_TOP T_BOTTOM

Reference: T(z) = T_bottom + (T_top - T_bottom) z / h (linear), and the total heat flow through
the face q A = k (T_top - T_bottom) / h * L * W.
Tolerance (fixed in docs/specification.md before the run): 0.5 % of the temperature difference.
The heat-flow comparison is a supplementary check with the same tolerance.
Writes results/v1_conduction.txt and exits with status 1 on failure.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ccxio
from material import MATERIAL

TOL = 0.005


def main():
    if len(sys.argv) != 5:
        raise SystemExit(__doc__)
    dat, mesh = sys.argv[1], sys.argv[2]
    t_top, t_bot = float(sys.argv[3]), float(sys.argv[4])
    par = ccxio.read_params(mesh)
    L, W, h = par["L"], par["W"], par["h"]
    nodes = ccxio.read_nodes(mesh)
    blocks = ccxio.parse_dat(dat)
    tb = [b for b in blocks if b["title"].startswith("temperatures") and b["set"] == "COL_CEN"]
    qb = [b for b in blocks if b["title"].startswith("total heat generation") and b["set"] == "TOP"]
    if not tb or not qb:
        raise SystemExit("missing temperature or heat-flow block in %s" % dat)
    rows = tb[-1]["rows"]
    z = np.array([nodes[int(r[0])][2] for r in rows])
    T = np.array([r[1] for r in rows])
    order = np.argsort(z)
    z, T = z[order], T[order]
    dT = t_top - t_bot
    T_exact = t_bot + dT * z / h
    err = np.abs(T - T_exact)
    err_pct = 100.0 * err.max() / abs(dT)
    q_exact = MATERIAL["k"] * dT / h
    Q_exact = q_exact * L * W
    Q_fe = qb[-1]["rows"][0][0]
    Q_err_pct = 100.0 * abs(abs(Q_fe) - abs(Q_exact)) / abs(Q_exact)

    lines = ["V1 steady conduction through the thickness",
             "job %s, mesh %s" % (os.path.basename(dat), os.path.basename(mesh)),
             "T_top %.3f degC, T_bottom %.3f degC, h = %g m, k = %g W/(m K)" % (t_top, t_bot, h, MATERIAL["k"]),
             "  z [mm]   T_FE [degC]   T_exact [degC]   |error| [K]"]
    for zi, ti, te, ei in zip(z, T, T_exact, err):
        lines.append("  %6.3f   %11.5f   %13.5f   %.2e" % (zi * 1e3, ti, te, ei))
    lines += ["max temperature error %.3e K = %.3e %% of the temperature difference (tolerance %.1f %%)" % (err.max(), err_pct, 100 * TOL),
              "heat flux  exact %.1f W/m2" % q_exact,
              "total heat flow through the face: FE %.2f W, exact %.2f W, error %.3e %% (tolerance %.1f %%)" % (Q_fe, Q_exact, Q_err_pct, 100 * TOL)]
    ok = err_pct <= 100 * TOL and Q_err_pct <= 100 * TOL
    lines.append("RESULT: %s" % ("PASS" if ok else "FAIL"))
    text = "\n".join(lines) + "\n"
    sys.stdout.write(text)
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    os.makedirs(os.path.join(root, "results"), exist_ok=True)
    with open(os.path.join(root, "results", "v1_conduction.txt"), "w") as f:
        f.write(text)
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
