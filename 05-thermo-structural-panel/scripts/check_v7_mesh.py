#!/usr/bin/python3
"""V7 criterion 4 (docs/specification.md, Amendment 3): mesh check.

Usage: /usr/bin/python3 scripts/check_v7_mesh.py
Compares the total centre deflection W/h of the A = 0.1 h run on the 40 x 40 x 2 mesh (results/v7_n40_a0p1_history.csv) with the
20 x 20 x 2 run (results/v7_a0p1_history.csv). Criterion: |W40/W20 - 1| <= 2 % at tau = 1.5. Reported also at other loads and for
the membrane force. Writes results/v7_mesh_check.txt. Exit status 1 on FAIL.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
root = os.path.dirname(HERE)
TOL = 0.02


def load(name):
    p = os.path.join(root, "results", name)
    if not os.path.exists(p):
        raise SystemExit("missing %s" % p)
    d = np.genfromtxt(p, delimiter=",", names=True)
    return d["tau"], d["W_over_h"], d["N_over_Ncr"]


def main():
    t20, w20, n20 = load("v7_a0p1_history.csv")
    t40, w40, n40 = load("v7_n40_a0p1_history.csv")
    out = ["V7 mesh check, A = 0.1 h: in-plane 40 x 40 against 20 x 20 (both nz = 2, C3D20R)",
           "  tau    W20/h    W40/h   change[%]   N20/Ncr   N40/Ncr   change[%]"]
    chg = {}
    for tau in (0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5):
        a, b = np.interp(tau, t20, w20), np.interp(tau, t40, w40)
        c, d = np.interp(tau, t20, n20), np.interp(tau, t40, n40)
        chg[tau] = 100.0 * (b / a - 1.0)
        out.append("  %.2f  %7.4f  %7.4f  %+9.3f  %9.4f  %8.4f  %+9.3f" % (tau, a, b, chg[tau], c, d, 100.0 * (d / c - 1.0)))
    ok = abs(chg[1.5]) <= 100 * TOL
    out.append("criterion 4: change of W at tau = 1.5 is %+.3f %% (tolerance %.0f %%)" % (chg[1.5], 100 * TOL))
    out.append("RESULT: %s" % ("PASS" if ok else "FAIL"))
    text = "\n".join(out) + "\n"
    print(text, end="")
    with open(os.path.join(root, "results", "v7_mesh_check.txt"), "w") as f:
        f.write(text)
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
