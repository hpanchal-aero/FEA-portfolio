#!/usr/bin/python3
"""Project 05: where is the maximum von Mises stress? Reads a CalculiX .dat that holds *EL PRINT S for ELSET PANEL.

Usage: /usr/bin/python3 scripts/locate_max_stress.py DAT_FILE BASE_MESH_INP
BASE_MESH_INP is the base mesh of the case (e.g. mesh/p5_h3_n20.inp); it supplies the undeformed element centroids.
Prints: the maximum over the pulse and its time, the 8 largest integration-point values of the last increment with
element centroid, distance from the plate centre and from the nearest edge, and the maxima by region
(interior / edge strip within 0.05 m / hot-spot region within 0.10 m of the centre) and by integration-point plane.
Integration points 1-4 and 5-8 are the lower and upper plane of the 2 x 2 x 2 rule (CalculiX ordering, element local z
along global z for this mesh); this convention is NOT independently verified here.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ccxio


def von_mises(s):
    sxx, syy, szz, sxy, sxz, syz = [s[:, k] for k in range(6)]
    return np.sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2) + 3.0 * (sxy ** 2 + sxz ** 2 + syz ** 2))


def main():
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    dat, mesh = sys.argv[1], sys.argv[2]
    par = ccxio.read_params(mesh)
    L, W, h = par["L"], par["W"], par["h"]
    nodes = ccxio.read_nodes(mesh)
    eids, conn = ccxio.read_elements(mesh)
    cen = {}
    for e, cn in zip(eids, conn):
        try:
            pts = np.array([nodes[int(c)] for c in list(cn)[:8]])
        except KeyError as err:
            raise SystemExit("element connectivity does not hold node ids as assumed: %s" % err)
        cen[int(e)] = pts.mean(axis=0)
    blocks = ccxio.parse_dat(dat)
    sb = sorted([b for b in blocks if b["title"].lower().startswith("stresses") and b["set"] == "PANEL"], key=lambda b: b["time"])
    if not sb:
        raise SystemExit("no stress blocks for PANEL in " + dat)
    peaks = []
    for b in sb:
        r = np.array(b["rows"])
        peaks.append(float(von_mises(r[:, 2:8]).max()) / 1e6)
    peaks = np.array(peaks)
    times = np.array([b["time"] for b in sb])
    print("plate %.3f x %.3f x %.4f m; %d stress increments, last at t = %g s" % (L, W, h, len(sb), times[-1]))
    print("maximum von Mises over the pulse %.3f MPa at t = %g s; last increment %.3f MPa (global maximum is the last: %s)" % (
        peaks.max(), times[int(np.argmax(peaks))], peaks[-1], bool(np.argmax(peaks) == len(peaks) - 1)))
    r = np.array(sb[-1]["rows"])
    vm = von_mises(r[:, 2:8]) / 1e6
    el = r[:, 0].astype(int)
    ip = r[:, 1].astype(int)
    cx = np.array([cen[e][0] for e in el])
    cy = np.array([cen[e][1] for e in el])
    rc = np.hypot(cx - L / 2, cy - W / 2)
    de = np.minimum(np.minimum(cx, L - cx), np.minimum(cy, W - cy))
    print("\n8 largest integration-point values at t = %g s (stresses in MPa)" % times[-1])
    print("  elem  ip plane   x_c    y_c   r_centre  d_edge     vM     sxx     syy     szz     sxy     sxz     syz")
    for k in np.argsort(vm)[::-1][:8]:
        s = r[k, 2:8] / 1e6
        print("%6d %3d %-5s %6.3f %6.3f %8.3f %7.3f %7.2f %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f" % (
            el[k], ip[k], "lower" if ip[k] <= 4 else "upper", cx[k], cy[k], rc[k], de[k], vm[k], *s))
    print("\nmaximum von Mises by region at t = %g s" % times[-1])
    for label, mask in (("interior (d_edge > 0.05 m)", de > 0.05), ("edge strip (d_edge <= 0.05 m)", de <= 0.05),
                        ("hot-spot region (r < 0.10 m)", rc < 0.10), ("lower plane (ip 1-4)", ip <= 4), ("upper plane (ip 5-8)", ip > 4)):
        print("  %-32s %8.2f MPa  (%d integration points)" % (label, vm[mask].max(), int(mask.sum())))


if __name__ == "__main__":
    main()
