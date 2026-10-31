#!/usr/bin/env python3
"""Spanwise load, deflection and elastic twist of several CalculiX static runs of the uCRM-9 coarse deck.
usage: plot_oneway_results.py <out.png> label=deck.inp [label=deck.inp ...]
Each deck needs its <job>.dat (displacements printed for NALL) next to it.
Twist per 1 m bin: plane fit uz = a + b (x - xm) + c (y - ym); nose-up twist = -b (positive when the leading edge rises relative to the trailing edge).
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import inspect_deck_loads as idl


def read_disp(dat):
    disp = {}
    mode = None
    with open(dat) as f:
        for ln in f:
            s = ln.strip()
            if not s:
                continue
            if s[0].isalpha():
                mode = "U" if s.startswith("displacements") else None
                continue
            if mode == "U":
                t = ln.split()
                if len(t) == 4:
                    try:
                        disp[int(t[0])] = (float(t[1]), float(t[2]), float(t[3]))
                    except ValueError:
                        pass
    return disp


def analyse(path):
    nodes, block_of, cload, _, _ = idl.parse(path)
    disp = read_disp(path[:-4] + ".dat")
    ids = [i for i in nodes if block_of[i] == "NALL" and i in disp]
    X = np.array([nodes[i] for i in ids])
    U = np.array([disp[i] for i in ids])
    lid = [i for i in cload if i in nodes]
    LY = np.array([nodes[i][1] for i in lid])
    LFz = np.array([cload[i][2] for i in lid])
    edges = np.arange(0.0, 30.0, 1.0)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (X[:, 1] >= lo) & (X[:, 1] < hi)
        ml = (LY >= lo) & (LY < hi)
        load = LFz[ml].sum() / 1000.0 / (hi - lo)
        if m.sum() < 12:
            rows.append((0.5 * (lo + hi), load, np.nan, np.nan))
            continue
        xm = X[m, 0].mean()
        ym = X[m, 1].mean()
        D = np.column_stack([np.ones(m.sum()), X[m, 0] - xm, X[m, 1] - ym])
        coef = np.linalg.lstsq(D, U[m, 2], rcond=None)[0]
        rows.append((ym, load, coef[0], np.degrees(-coef[1])))
    # tip bin: nodes with y >= ymax - 1.0
    ymax = X[:, 1].max()
    mt = X[:, 1] >= ymax - 1.0
    xm = X[mt, 0].mean()
    ym = X[mt, 1].mean()
    D = np.column_stack([np.ones(mt.sum()), X[mt, 0] - xm, X[mt, 1] - ym])
    coef = np.linalg.lstsq(D, U[mt, 2], rcond=None)[0]
    tip = dict(uz_max=U[:, 2].max(), tip_mean=U[X[:, 1] >= ymax - 0.001, 2].mean(), twist_tip=np.degrees(-coef[1]), ymax=ymax)
    return np.array(rows), tip


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    out = sys.argv[1]
    cases = []
    for a in sys.argv[2:]:
        label, path = a.split("=", 1)
        rows, tip = analyse(path)
        cases.append((label, rows, tip))
        print("%-10s tip mean uz %.4f m   max uz %.4f m   tip twist (last 1 m) %.3f deg nose-up   ymax %.3f" % (
            label, tip["tip_mean"], tip["uz_max"], tip["twist_tip"], tip["ymax"]))
    labels = [c[0] for c in cases]
    print("\n   y   " + " ".join("load_%-9s" % l for l in labels) + " | " + " ".join("uz_%-10s" % l for l in labels) + " | " + " ".join("tw_%-10s" % l for l in labels))
    n = len(cases[0][1])
    for k in range(n):
        if k % 2 == 1 and k not in (n - 1,):
            continue
        y = cases[0][1][k, 0]
        print("%6.2f " % y + " ".join("%13.2f" % c[1][k, 1] for c in cases) + " | " + " ".join("%13.4f" % c[1][k, 2] for c in cases) + " | " + " ".join("%13.3f" % c[1][k, 3] for c in cases))
    ref = cases[0]
    for lab, rows, tip in cases[1:]:
        r = tip["tip_mean"] / ref[2]["tip_mean"] - 1.0
        print("\n%s vs %s: tip mean uz %+.2f %%, tip twist difference %+.3f deg" % (lab, ref[0], 100 * r, tip["twist_tip"] - ref[2]["twist_tip"]))

    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.4))
    styles = [dict(color="k", lw=2.0), dict(color="C0", lw=2.0), dict(color="C1", lw=1.6, ls="--")]
    for (lab, rows, tip), st in zip(cases, styles):
        ax[0].plot(rows[:, 0], rows[:, 1], label=lab, **st)
        ax[1].plot(rows[:, 0], rows[:, 2], label="%s (tip %.3f m)" % (lab, tip["tip_mean"]), **st)
        ax[2].plot(rows[:, 0], rows[:, 3], label=lab, **st)
    ax[0].set_ylabel("Fz applied per unit span [kN/m]")
    ax[1].set_ylabel("uz [m]")
    ax[2].set_ylabel("elastic twist, nose-up positive [deg]")
    ax[0].set_title("Aerodynamic load per unit span")
    ax[1].set_title("Vertical deflection")
    ax[2].set_title("Elastic twist (plane fit per 1 m)")
    for a_ in ax:
        a_.set_xlabel("y [m]")
        a_.grid(alpha=0.3)
        a_.legend(fontsize=8)
    fig.suptitle("uCRM-9 coarse wingbox, 1g: dataset load vs CFD-mapped load (M 0.85, first-order upwind CFD loads)")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print("\nsaved", out)


if __name__ == "__main__":
    main()
