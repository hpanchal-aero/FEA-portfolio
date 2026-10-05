#!/usr/bin/env python3
"""
Stage 4: skin stress comparison against a baseline job, using the Stage 3 definitions unchanged
(collect_stage3_results.stress_stats / dte): design skin stress = max VM over expanded skin nodes with
dTE >= 3 mm and z > 10 mm. Skin stress = CalculiX expanded-node values (interpretation verified near
the root only). Yield 290 MPa is ASSUMED (uncited).

Usage (system python):
  /usr/bin/python3 scripts/compare_stage4_stress.py <baseline_job> <job2> [<job3> ...]
Jobs are paths without extension. Needs <job>.inp and <job>.frd. Baseline and jobs must share a mesh.
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collect_stage3_results import stress_stats, dte

YIELD = 290.0


def load(job):
    stats, (vm, eid, ex, so, coords) = stress_stats(job)
    return stats, vm[eid], eid, ex


def groups(ex):
    n = len(ex)
    return (("all", np.ones(n, bool)), ("upper", ex[:, 1] > 0.0), ("lower", ex[:, 1] <= 0.0))


def report(job, stats, v, ex):
    d = dte(ex[:, 2], ex[:, 0])
    des = (d >= 3.0) & (ex[:, 2] > 10.0)
    print(f"\n=== {job}")
    print(f"  stress_stats: p99 {stats['skin_p99']:.3f} | p99.9 {stats['skin_p99_9']:.3f} | max {stats['skin_max']:.3f} | "
          f"design(dTE>=3) {stats['skin_max_dTE_ge3']:.3f} | dTE>=15 {stats['skin_max_dTE_ge15']:.3f} | "
          f"solid-only max {stats['solid_only_max']:.3f}")
    for lab, g in groups(ex):
        mk = des & g
        k = int(np.argmax(np.where(mk, v, -1.0)))
        print(f"  {lab:6s} design max {v[k]:8.3f} MPa at (x={ex[k,0]:.1f}, y={ex[k,1]:.2f}, z={ex[k,2]:.1f}) "
              f"| p99 {np.percentile(v[g], 99):7.3f} | design nodes {int(mk.sum())} "
              f"| FoS at design max {YIELD / v[k]:.2f}")
    pan = (ex[:, 1] > 0.0) & (ex[:, 0] >= 40.0) & (ex[:, 0] <= 160.0) & (ex[:, 2] >= 10.0) & (ex[:, 2] <= 400.0)
    print(f"  upper aft panel (x 40-160, z 10-400, n={int(pan.sum())}): VM mean {v[pan].mean():.3f}, "
          f"p95 {np.percentile(v[pan], 95):.3f}, max {v[pan].max():.3f}")


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(2)
    jobs = sys.argv[1:]
    s0, v0, id0, ex0 = load(jobs[0])
    report(jobs[0], s0, v0, ex0)
    for j in jobs[1:]:
        s, v, ids, ex = load(j)
        report(j, s, v, ex)
        if not (np.array_equal(ids, id0) and np.allclose(ex, ex0, atol=1e-9)):
            print("  WARNING: node sets differ from baseline; node-by-node comparison skipped")
            continue
        d = dte(ex[:, 2], ex[:, 0])
        des = (d >= 3.0) & (ex[:, 2] > 10.0)
        dv = v - v0
        kmax = int(np.argmax(np.where(des, dv, -1e30)))
        kmin = int(np.argmin(np.where(des, dv, 1e30)))
        print(f"  vs baseline over design nodes (n={int(des.sum())}): mean |dVM| {np.abs(dv[des]).mean():.3f} MPa, "
              f"rms {np.sqrt(np.mean(dv[des] ** 2)):.3f}")
        print(f"    largest increase {dv[kmax]:+.3f} MPa at (x={ex[kmax,0]:.1f}, y={ex[kmax,1]:.2f}, z={ex[kmax,2]:.1f}) "
              f"[{v0[kmax]:.2f} -> {v[kmax]:.2f}]")
        print(f"    largest decrease {dv[kmin]:+.3f} MPa at (x={ex[kmin,0]:.1f}, y={ex[kmin,1]:.2f}, z={ex[kmin,2]:.1f}) "
              f"[{v0[kmin]:.2f} -> {v[kmin]:.2f}]")
        for lab, g in groups(ex):
            mk = des & g
            print(f"    {lab:6s} design-node mean VM {v0[mk].mean():.3f} -> {v[mk].mean():.3f} MPa")


if __name__ == "__main__":
    main()
