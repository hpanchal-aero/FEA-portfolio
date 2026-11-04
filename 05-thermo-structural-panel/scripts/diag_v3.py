#!/usr/bin/python3
"""Diagnostic for V3: element stresses along the centre row of elements, bottom layer.

Usage: /usr/bin/python3 scripts/diag_v3.py JOB.dat MESH.inp

Prints, for selected elements along x (index 0 = at the x = 0 edge, 19 = at the plate centre) in the
element row closest to y = W/2 and the bottom through-thickness layer, the range of sigma_x, sigma_z
and sigma_xz over the 8 integration points of each element.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ccxio


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    dat, mesh = sys.argv[1], sys.argv[2]
    par = ccxio.read_params(mesh)
    L, W, h = par["L"], par["W"], par["h"]
    nodes = ccxio.read_nodes(mesh)
    eids, conn = ccxio.read_elements(mesh)
    blocks = ccxio.parse_dat(dat)
    sel = [b for b in blocks if b["title"].lower().startswith("stresses") and b["set"] == "PANEL"]
    if not sel:
        raise SystemExit("no stress block for set PANEL")
    s = np.array(sel[-1]["rows"])
    cen = np.array([np.mean([nodes[int(n)] for n in c[:8]], axis=0) for c in conn])
    iy = int(np.argmin(np.abs(cen[:, 1] - W / 2)))
    ycen = cen[iy, 1]
    row = np.nonzero((np.abs(cen[:, 1] - ycen) < 1e-9) & (cen[:, 2] < h / 2))[0]
    row = row[np.argsort(cen[row, 0])]
    print("job %s, row of elements at y = %.4f m, bottom layer (z < h/2), %d elements along x" % (os.path.basename(dat), ycen, len(row)))
    print("  idx  x_centroid[mm]      sxx min        sxx max        szz min        szz max     max|sxz|   [Pa]")
    for k in [0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 19]:
        if k >= len(row):
            continue
        e = int(eids[row[k]])
        v = s[s[:, 0] == e]
        sxx, szz, sxz = v[:, 2], v[:, 4], v[:, 6]
        print("  %3d  %12.3f   %13.5e  %13.5e  %13.5e  %13.5e  %10.3e"
              % (k, cen[row[k], 0] * 1e3, sxx.min(), sxx.max(), szz.min(), szz.max(), np.abs(sxz).max()))


if __name__ == "__main__":
    main()
