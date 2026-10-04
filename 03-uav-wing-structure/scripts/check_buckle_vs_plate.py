#!/usr/bin/env python3
"""
Plausibility check of a skin-bay buckling mode against simply-supported flat-plate theory.

  sigma_cr_FEA   = buckling factor x mean compressive sigma_zz of the static preload over the
                   mode's hot region (|U| > 50% of the mode peak)
  sigma_cr_plate = k pi^2 E / (12 (1-nu^2)) (t/b)^2,  k = min_m (m b/a + a/(m b))^2,
                   b = width of the hot region in x (2nd-98th percentile), a = rib spacing 200 mm

NOT a verification: curvature, the real edge conditions and stress variation across the panel
are ignored; agreement is expected only to ~10-25%. The static job must be the same mesh/model.

Usage (system Python):
  /usr/bin/python3 scripts/check_buckle_vs_plate.py <buckle_job> <static_job> <factor> [block=2]
  (block 2 = buckling mode 1; block 1 is the static preload)
"""
import math
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_buckle_modes import read_modes
from extract_stage3_skin import read_frd, read_frd_nodes

E_MOD, NU, A_BAY = 73100.0, 0.33, 200.0


def skin_t(inp):
    with open(inp) as f:
        for line in f:
            if line.upper().startswith("*SHELL SECTION"):
                return float(next(f).split(",")[0])
    sys.exit("ERROR: *SHELL SECTION not found")


def main():
    bjob, sjob, factor = sys.argv[1], sys.argv[2], float(sys.argv[3])
    block = int(sys.argv[4]) if len(sys.argv) > 4 else 2
    t = skin_t(sjob + ".inp")
    ids, xyz, blocks = read_modes(bjob + ".frd")
    blk = blocks[block - 1]
    bid = np.array([b[0] for b in blk])
    U = np.array([b[1] for b in blk])
    mag = np.linalg.norm(U, axis=1)
    hot = bid[mag > 0.5 * mag.max()]
    o = np.argsort(ids)
    pos = np.searchsorted(ids[o], hot)
    P = xyz[o][pos]

    sid, sxyz = read_frd_nodes(sjob + ".frd")
    so = np.argsort(sid)
    spos = np.searchsorted(sid[so], hot)
    dmax = np.abs(sxyz[so][spos] - P).max()
    print(f"coordinate match between buckle and static node sets (hot nodes): max diff {dmax:.2e} mm")
    if dmax > 1e-6:
        sys.exit("ERROR: node numbering differs between jobs -- STOP.")

    S = read_frd(sjob + ".frd", "STRESS", 6)
    sz = S[hot, 2]
    sig_app = -float(sz.mean())
    lo, hi = np.percentile(-sz, [10, 90])
    xlo, xhi = np.percentile(P[:, 0], [2, 98])
    zlo, zhi = np.percentile(P[:, 2], [2, 98])
    b = xhi - xlo
    ks = [(m * b / A_BAY + A_BAY / (m * b)) ** 2 for m in range(1, 8)]
    k = min(ks)
    m_best = ks.index(k) + 1
    base = math.pi ** 2 * E_MOD / (12 * (1 - NU ** 2)) * (t / b) ** 2
    sig_plate = k * base
    sig_fea = factor * sig_app
    print(f"skin t = {t} mm; hot region: {len(hot)} nodes, x {xlo:.0f}..{xhi:.0f} (b = {b:.1f} mm), "
          f"z {zlo:.0f}..{zhi:.0f}")
    print(f"static compressive sigma_zz over hot region: mean {sig_app:.2f} MPa (10-90%: {lo:.2f}..{hi:.2f})")
    print(f"FEA critical stress = factor {factor} x {sig_app:.2f} = {sig_fea:.1f} MPa")
    print(f"SS plate theory: k = {k:.3f} (m = {m_best}), sigma_cr = {sig_plate:.1f} MPa "
          f"(a = {A_BAY:.0f}, b = {b:.1f})")
    print(f"FEA / plate = {sig_fea / sig_plate:.3f}  (range using 10-90% stress spread: "
          f"{factor * lo / sig_plate:.2f}..{factor * hi / sig_plate:.2f})")


if __name__ == "__main__":
    main()
