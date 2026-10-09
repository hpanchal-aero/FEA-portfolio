#!/usr/bin/env python3
"""
List the smallest skin pieces of an imprint face list (mesh/stage3_skin/stage3_skin_faces[_rs<N>].json),
with centroid, shared flag and the distance of the centroid to the nearest rib centre plane.
Usage: /usr/bin/python3 scripts/diag_rs_slivers.py <faces.json> <rib_spacing_mm> [<n_rows>]
"""
import json
import sys
import numpy as np


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(2)
    path, spacing = sys.argv[1], float(sys.argv[2])
    nrow = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    P = json.load(open(path))
    a = np.array([p["area"] for p in P])
    order = np.argsort(a)
    print(f"\n{path}: {len(P)} pieces; area < 0.5 mm^2: {int((a < 0.5).sum())}, "
          f"< 1.0: {int((a < 1.0).sum())}, < 5: {int((a < 5).sum())}")
    print(" rank  area[mm2]  shared      cx        cy        cz   |cz - nearest rib centre|  list idx")
    for r, k in enumerate(order[:nrow], 1):
        p = P[k]
        zr = round(p["cz"] / spacing) * spacing
        print(f"{r:5d}  {p['area']:9.4f}  {str(p['shared']):6s}  {p['cx']:8.3f}  {p['cy']:8.3f}  {p['cz']:8.3f}"
              f"  {abs(p['cz'] - zr):10.3f}  {k:8d}")


if __name__ == "__main__":
    main()
