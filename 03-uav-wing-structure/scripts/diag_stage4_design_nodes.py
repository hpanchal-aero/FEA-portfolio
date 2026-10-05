#!/usr/bin/env python3
"""
Stage 4: where is the 'design' skin stress maximum? Uses collect_stage3_results.stress_stats / dte
(definitions unchanged: design = expanded skin nodes with dTE >= 3 mm and z > 10 mm).
Usage (system python): /usr/bin/python3 scripts/diag_stage4_design_nodes.py <job> [<job2> ...]
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collect_stage3_results import stress_stats, dte

RIBS = np.array([0.0, 200.0, 400.0, 600.0, 800.0, 1000.0])


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    for job in sys.argv[1:]:
        st, (vm, eid, ex, so, coords) = stress_stats(job)
        v = vm[eid]
        d = dte(ex[:, 2], ex[:, 0])
        des = (d >= 3.0) & (ex[:, 2] > 10.0)
        print(f"\n=== {job}: design max {v[des].max():.3f} MPa, design nodes {int(des.sum())}")
        idx = np.nonzero(des)[0]
        top = idx[np.argsort(-v[idx])[:12]]
        print("  rank   VM [MPa]    x       y       z     dTE    |z-rib|")
        for r, k in enumerate(top, 1):
            dr = np.abs(ex[k, 2] - RIBS).min()
            print(f"  {r:4d}  {v[k]:8.3f}  {ex[k,0]:7.1f} {ex[k,1]:7.2f} {ex[k,2]:7.1f}  {d[k]:6.1f}  {dr:6.1f}")
        for zc in (995.0, 990.0, 950.0):
            mk = des & (ex[:, 2] < zc)
            print(f"  design max excluding z >= {zc:.0f}: {v[mk].max():.3f} MPa")
        for dc in (6.0, 15.0):
            mk = des & (d >= dc)
            print(f"  design max with dTE >= {dc:.0f} (all z>10): {v[mk].max():.3f} MPa")
        mk = des & (ex[:, 2] < 990.0) & (d >= 6.0)
        print(f"  design max with z < 990 and dTE >= 6: {v[mk].max():.3f} MPa")


if __name__ == "__main__":
    main()
