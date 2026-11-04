#!/usr/bin/python3
"""Baseline elastic stress check (docs/specification.md, Amendment 11).

Usage: /usr/bin/python3 scripts/post_baseline_stress.py [DAT]
Reads the stress blocks (element, integration point, sxx, syy, szz, sxy, sxz, syz in Pa) of the baseline .dat, computes the von Mises stress
at every integration point and reports, per increment, the maximum, its element, integration point and the element centroid.
Stop rule: 1.267 x (maximum von Mises) >= 660 MPa at any increment -> exit status 1 (ask before continuing).
Writes results/p05_baseline_stress.txt and results/p05_baseline_stress_history.csv.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ccxio
from material import MATERIAL

SY = MATERIAL["sigma_y"]
BOUND = 1.0 / 0.789
STOP = 0.75 * SY


def main():
    root = os.path.dirname(HERE)
    dat = sys.argv[1] if len(sys.argv) > 1 else os.path.join(root, "simulation", "p05_base_struct", "p05_base_struct.dat")
    mesh = os.path.join(root, "mesh", "v7_n20_z2_a0p1.inp")
    par = ccxio.read_params(os.path.join(root, "mesh", "v6c_n20_z2.inp"))
    h = par["h"]
    nodes = ccxio.read_nodes(mesh)
    eids, conn = ccxio.read_elements(mesh)
    cen = dict((int(e), np.mean([nodes[int(n)] for n in c[:8]], axis=0)) for e, c in zip(eids, conn))
    blocks = [b for b in ccxio.parse_dat(dat) if b["title"].lower().startswith("stresses") and b["set"] == "PANEL"]
    blocks.sort(key=lambda b: b["time"])
    if not blocks:
        raise SystemExit("no stress blocks found in %s" % dat)
    out = ["Baseline elastic stress check: %d stress blocks, yield stress %.0f MPa, stop bound %.1f MPa on 1.267 x max von Mises" % (len(blocks), SY / 1e6, STOP / 1e6),
           "   t[s]  max vM [MPa]  bound [MPa]  element ip  layer   centroid (x, y, z) [m]      min(sxx,syy) [MPa]  max|szz| [MPa]"]
    hist = ["time_s,max_vM_MPa,bound_MPa,element,ip,min_sxx_syy_MPa,max_abs_szz_MPa"]
    worst = 0.0
    for b in blocks:
        r = np.array(b["rows"])
        if r.shape[1] != 8:
            raise SystemExit("stress rows have %d columns, expected 8" % r.shape[1])
        sxx, syy, szz, sxy, sxz, syz = [r[:, k] for k in range(2, 8)]
        vm = np.sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2) + 3.0 * (sxy ** 2 + sxz ** 2 + syz ** 2))
        j = int(np.argmax(vm))
        e = int(r[j, 0])
        c = cen[e]
        smin = min(sxx.min(), syy.min())
        line = "  %6.3f %12.3f %12.3f %8d %3d  %-5s (%.4f, %.4f, %.5f) %20.3f %15.3f" % (
            b["time"], vm[j] / 1e6, BOUND * vm[j] / 1e6, e, int(r[j, 1]), "top" if c[2] > h / 2 else "bot", c[0], c[1], c[2],
            smin / 1e6, np.abs(szz).max() / 1e6)
        hist.append("%.6f,%.6f,%.6f,%d,%d,%.6f,%.6f" % (b["time"], vm[j] / 1e6, BOUND * vm[j] / 1e6, e, int(r[j, 1]), smin / 1e6, np.abs(szz).max() / 1e6))
        if abs(b["time"] / 2.5 - round(b["time"] / 2.5)) < 1e-6:
            out.append(line)
        worst = max(worst, BOUND * vm[j])
    imax = int(np.argmax([float(h_.split(",")[1]) for h_ in hist[1:]]))
    out.append("overall: max von Mises %s MPa at t = %s s; largest bound %.1f MPa = %.1f %% of yield, %.1f %% of the %.0f MPa stop bound"
               % (hist[1:][imax].split(",")[1], hist[1:][imax].split(",")[0], worst / 1e6, 100 * worst / SY, 100 * worst / STOP, STOP / 1e6))
    triggered = worst >= STOP
    out.append("STOP RULE: %s" % ("TRIGGERED, ask before continuing" if triggered else "not triggered, the elastic model stands for the baseline"))
    print("\n".join(out))
    with open(os.path.join(root, "results", "p05_baseline_stress.txt"), "w") as f:
        f.write("\n".join(out) + "\n")
    with open(os.path.join(root, "results", "p05_baseline_stress_history.csv"), "w") as f:
        f.write("\n".join(hist) + "\n")
    if triggered:
        sys.exit(1)


if __name__ == "__main__":
    main()
