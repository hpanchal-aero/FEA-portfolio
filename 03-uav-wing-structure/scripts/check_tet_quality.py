#!/usr/bin/env python3
"""
Tet quality of a raw Netgen mesh (.npz with 'points'/'P' and 'tets', mesh_stage3_skin_raw.py output).
Metric: mean ratio q = 12 (3V)^(2/3) / sum(edge^2)  (= 1 for a regular tetrahedron, -> 0 for flat tets),
and the longest/shortest edge ratio. Reports the distribution and the worst tets.
Usage: /usr/bin/python3 scripts/check_tet_quality.py <raw.npz> [<raw2.npz> ...]
"""
import sys
import numpy as np


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    for path in sys.argv[1:]:
        d = np.load(path)
        key = "P" if "P" in d.files else "points"
        P, T = d[key], d["tets"]
        print(f"\n=== {path}: {key} {P.shape}, tets {T.shape}, tet index range {T.min()}..{T.max()}")
        if T.ndim != 2 or T.shape[1] != 4:
            print("ERROR: tets is not (N,4)")
            sys.exit(1)
        if T.max() >= len(P) or T.min() < 0:
            print("ERROR: tet indices out of range")
            sys.exit(1)
        A, B, C, D = (P[T[:, k]] for k in range(4))
        V = np.abs(np.einsum("ij,ij->i", np.cross(B - A, C - A), D - A)) / 6.0
        pairs = [(A, B), (A, C), (A, D), (B, C), (B, D), (C, D)]
        l2 = np.stack([np.sum((x - y) ** 2, axis=1) for x, y in pairs], axis=1)
        q = 12.0 * (3.0 * V) ** (2.0 / 3.0) / l2.sum(axis=1)
        er = np.sqrt(l2.max(axis=1) / l2.min(axis=1))
        pct = lambda a, p: float(np.percentile(a, p))
        print(f"  volume [mm^3]: min {V.min():.6f}, median {np.median(V):.4f}, max {V.max():.3f}")
        print(f"  mean ratio q : min {q.min():.4f}, p0.1 {pct(q, 0.1):.4f}, p1 {pct(q, 1):.4f}, "
              f"p5 {pct(q, 5):.4f}, median {np.median(q):.4f}")
        print(f"  tets with q < 0.02: {int((q < 0.02).sum())}, < 0.05: {int((q < 0.05).sum())}, "
              f"< 0.10: {int((q < 0.10).sum())}  (of {len(q)})")
        print(f"  edge ratio longest/shortest: median {np.median(er):.2f}, p99 {pct(er, 99):.2f}, max {er.max():.2f}")
        cen = (A + B + C + D) / 4.0
        print("  worst 5 tets by q:   q        volume      centroid (x, y, z)")
        for k in np.argsort(q)[:5]:
            print(f"                    {q[k]:.5f}  {V[k]:.6f}   ({cen[k,0]:.2f}, {cen[k,1]:.2f}, {cen[k,2]:.2f})")


if __name__ == "__main__":
    main()
