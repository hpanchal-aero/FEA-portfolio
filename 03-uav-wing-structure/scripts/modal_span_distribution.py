#!/usr/bin/env python3
"""
Spanwise distribution of buckling-mode displacement energy (sum |u|^2 per 100 mm span bin) from a CalculiX .frd.
Usage: /usr/bin/python3 scripts/modal_span_distribution.py L100 L126 L006
Reads simulation/stage4_s4_k8_<ID>_closeedgefacnone_t1_n3p8_buckle.frd.
The file holds 1 preload block (factor 0) followed by the 6 buckling modes; the 100CL lines give the factors,
which are asserted below. Self-contained .frd parser: coordinates block (2C) and DISP blocks, 12-char fields.
"""
import sys

import numpy as np

NMODES = 3          # modes reported
BIN = 100.0
ZMAX = 1000.0


def dte(z, x):
    """Distance to the trailing edge at span station z (same definition as collect_stage3_results.dte)."""
    return 0.7 * (300.0 - 0.15 * z) - x


def parse(path):
    coords, disp, factors = {}, [], []
    mode = None          # "coord" | "disp" | None
    cur = None
    with open(path, "r", errors="replace") as f:
        for ln in f:
            if ln.startswith("  100CL"):
                factors.append(float(ln.split()[2]))
            elif ln.startswith("    2C"):
                mode, cur = "coord", coords
            elif ln.startswith(" -4  DISP"):
                mode = "disp"
                cur = {}
                disp.append(cur)
            elif ln.startswith(" -4 "):
                mode, cur = None, None
            elif ln.startswith(" -3"):
                mode, cur = None, None
            elif ln.startswith(" -1") and mode is not None:
                nid = int(ln[3:13])
                cur[nid] = [float(ln[13 + 12 * k: 25 + 12 * k]) for k in range(3)]
    return coords, disp, factors


def main():
    ids = sys.argv[1:]
    if not ids:
        sys.exit("usage: modal_span_distribution.py L100 L126 ...")
    edges = np.arange(0.0, ZMAX + BIN, BIN)
    for L in ids:
        path = "simulation/stage4_s4_k8_%s_closeedgefacnone_t1_n3p8_buckle.frd" % L
        coords, disp, factors = parse(path)
        print("\n=== %s: %d nodes, %d DISP blocks, factors %s ===" % (L, len(coords), len(disp), factors))
        if len(coords) == 0 or len(disp) != 7 or len(factors) != 7 or factors[0] != 0.0:
            sys.exit("ABORT: expected 7 DISP blocks / 7 factors with the first = 0 (preload) in %s" % path)
        if any(factors[i] <= 0 for i in range(1, 7)):
            sys.exit("ABORT: non-positive buckling factor in %s" % path)
        for m in range(1, NMODES + 1):
            d = disp[m]
            nid = np.array([n for n in d if n in coords])
            if len(nid) < 0.99 * len(d):
                sys.exit("ABORT: %d of %d displaced nodes lack coordinates" % (len(d) - len(nid), len(d)))
            xyz = np.array([coords[n] for n in nid])
            u2 = np.array([sum(v * v for v in d[n]) for n in nid])
            frac, _ = np.histogram(xyz[:, 2], bins=edges, weights=u2)
            frac = frac / u2.sum()
            ip = int(np.argmax(u2))
            zc = float((xyz[:, 2] * u2).sum() / u2.sum())
            print(" mode %d (factor %.4f): peak |u| at x=%.1f y=%.1f z=%.1f (%.0f mm ahead of TE) ; centroid z = %.0f mm"
                  % (m, factors[m], xyz[ip, 0], xyz[ip, 1], xyz[ip, 2], dte(xyz[ip, 2], xyz[ip, 0]), zc))
            print("   energy %% per 100 mm bin (z=0..1000): " + " ".join("%4.1f" % (100 * x) for x in frac))
            print("   cumulative z<200: %.1f %%   z<400: %.1f %%" % (100 * frac[:2].sum(), 100 * frac[:4].sum()))


if __name__ == "__main__":
    main()
