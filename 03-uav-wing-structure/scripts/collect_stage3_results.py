#!/usr/bin/env python3
"""
Collect the Stage 3 results from the local solver output into tables (results/*.csv, .json).

  results/stage3_static_summary.csv    limit-load static results per job (tip deflection from the
                                       per-job JSON; skin / solid stress statistics from the .frd)
  results/stage3_buckling_summary.csv  buckling factors (6 modes) per job, from the .dat files
  results/stage3_te_refinement.csv     trailing-edge refinement study at rib z=800
  results/stage3_summary.json          all of the above

Skin stress = CalculiX expanded-node values (interpretation verified near the root only).
"design" skin stress = max VM over skin nodes with dTE >= 3 mm and z > 10 mm; peaks closer to the
trailing edge at the rib stations are a mesh-dependent singularity (see the TE refinement table).
Yield 290 MPa is an ASSUMED, uncited value.

Run with the system Python (needs the local simulation/ and mesh/stage3_skin/*.npz files):
  /usr/bin/python3 scripts/collect_stage3_results.py
"""
import csv
import json
import os
import re
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from extract_stage3_skin import parse_inp, read_frd, read_frd_nodes, von_mises

YIELD = 290.0
REQ_BUCKLE = 1.5
ROW_RE = re.compile(r"^\s*(\d+)\s+([-+]?\d*\.?\d+[EeDd][+-]?\d+)\s*$")
RIBS = (200.0, 400.0, 600.0, 800.0, 1000.0)


def dte(z, x):
    return 0.7 * (300.0 - 0.15 * z) - x


def stress_stats(job):
    xyz, solid, shell, root, cload, t = parse_inp(job + ".inp")
    n0 = len(xyz)
    coords = np.vstack([np.zeros((1, 3)), xyz])
    S = read_frd(job + ".frd", "STRESS", 6)
    vm = von_mises(S)
    fid, fxyz = read_frd_nodes(job + ".frd")
    m = fid > n0
    eid, ex = fid[m], fxyz[m]
    v = vm[eid]
    d = dte(ex[:, 2], ex[:, 0])
    far = ex[:, 2] > 10.0
    so = np.setdiff1d(np.unique(solid), np.unique(shell))
    return dict(
        skin_p99=float(np.percentile(v, 99)), skin_p99_9=float(np.percentile(v, 99.9)),
        skin_max=float(v.max()), skin_max_z_gt10=float(v[far].max()),
        skin_max_dTE_ge3=float(v[(d >= 3.0) & far].max()),
        skin_max_dTE_ge15=float(v[(d >= 15.0) & far].max()),
        solid_only_max=float(vm[so].max()),
    ), (vm, eid, ex, so, coords)


def tag_of(job):
    return os.path.basename(job)


def static_table():
    jobs = [("s2p5_k8", t) for t in ("0p05", "0p5", "1", "1p5")] + [("s4_k8", "1"), ("s4_k6", "1")]
    rows = []
    for mesh, t in jobs:
        job = f"simulation/stage3_skin_{mesh}_t{t}_n3p8"
        jf = f"results/stage3_skin_{mesh}_t{t}_n3p8_results.json"
        if not (os.path.isfile(job + ".frd") and os.path.isfile(jf)):
            print(f"  skip {job} (missing files)"); continue
        j = json.load(open(jf))
        st, _ = stress_stats(job)
        row = dict(job=tag_of(job), mesh=mesh, skin_t_mm=j["t_skin_mm"], lift_N=j["lift_N"], n=j["n"],
                   W_solid_N=j["W_solid_N"], W_skin_N=j["W_skin_N"],
                   tip_uy_mm=j["uy_spar_axis_mm"], beam_B_lo_mm=j["bracket_lo_mm"],
                   beam_B_hi_mm=j["bracket_hi_mm"],
                   tip_vs_B_hi_pct=100 * (j["uy_spar_axis_mm"] / j["bracket_hi_mm"] - 1.0))
        row.update(st)
        row["FoS_design_skin_assumed_yield"] = YIELD / st["skin_max_dTE_ge3"]
        rows.append(row)
        print(f"  static {tag_of(job)} done")
    return rows


def buckling_table():
    rows = []
    for mesh, t in (("s4_k8", "0p5"), ("s4_k8", "1"), ("s4_k8", "1p5"), ("s4_k6", "0p5"), ("s4_k6", "1")):
        job = f"simulation/stage3_skin_{mesh}_t{t}_n3p8_buckle"
        if not os.path.isfile(job + ".dat"):
            print(f"  skip {job} (missing .dat)"); continue
        lines = open(job + ".dat").readlines()
        hi = next(i for i, l in enumerate(lines) if "BUCKLINGFACTOROUTPUT" in "".join(l.upper().split()))
        fac = []
        for l in lines[hi:hi + 60]:
            m = ROW_RE.match(l)
            if m:
                fac.append(float(m.group(2).replace("D", "E").replace("d", "E")))
            elif fac and l.strip() == "":
                break
        crit = min(fac)
        row = dict(job=tag_of(job), mesh=mesh, skin_t_mm={"0p5": 0.5, "1": 1.0, "1p5": 1.5}[t])
        for k, f in enumerate(fac, 1):
            row[f"mode{k}"] = f
        row.update(critical=crit, requirement=REQ_BUCKLE, margin_pct=100 * (crit / REQ_BUCKLE - 1.0),
                   result="PASS" if crit >= REQ_BUCKLE else "FAIL")
        rows.append(row)
        print(f"  buckling {tag_of(job)}: critical {crit:.4f}")
    return rows


def te_table():
    cases = [("baseline", None, "stage3_skin_s4_k8_t1_n3p8", "raw_s4_k8"),
             ("te_h=1", 1.0, "stage3_skin_s4_k8_te1_z800_t1_n3p8", "raw_s4_k8_te1_z800"),
             ("te_h=0.5", 0.5, "stage3_skin_s4_k8_te0p5_z800_t1_n3p8", "raw_s4_k8_te0p5_z800"),
             ("te_h=0.25", 0.25, "stage3_skin_s4_k8_te0p25_z800_t1_n3p8", "raw_s4_k8_te0p25_z800")]
    rows = []
    for label, h, job, raw in cases:
        job = "simulation/" + job
        rawp = f"mesh/stage3_skin/{raw}.npz"
        if not (os.path.isfile(job + ".frd") and os.path.isfile(rawp)):
            print(f"  skip {label} (missing files)"); continue
        d = np.load(rawp)
        P, tris = d["P"], d["tris"]
        c0 = np.array([0.7 * (300.0 - 0.15 * 800.0), 0.23, 800.0])
        e = np.concatenate([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [0, 2]]])
        pa, pb = P[e[:, 0]], P[e[:, 1]]
        near = (np.linalg.norm(pa - c0, axis=1) < 6.0) & (np.linalg.norm(pb - c0, axis=1) < 6.0)
        L = np.linalg.norm(pa[near] - pb[near], axis=1)
        _, (vm, eid, ex, so, coords) = stress_stats(job)
        v = vm[eid]
        d800 = dte(ex[:, 2], ex[:, 0])
        w = np.abs(ex[:, 2] - 800.0) < 2.0
        row = dict(case=label, te_edge_size_requested_mm=h, local_median_edge_mm=float(np.median(L)),
                   local_min_edge_mm=float(L.min()), n_local_edges=int(near.sum()))
        win = w & (d800 < 6.0)
        row["skin_peak_rib800"] = float(v[win].max())
        sw = (np.abs(coords[so, 2] - 800.0) < 2.0) & (dte(coords[so, 2], coords[so, 0]) < 6.0)
        row["solid_only_peak_rib800"] = float(vm[so][sw].max())
        for lo, hi in ((1.0, 3.0), (3.0, 6.0), (6.0, 15.0), (15.0, 40.0)):
            mm = w & (d800 >= lo) & (d800 < hi)
            row[f"ring_{lo:g}_{hi:g}mm"] = float(v[mm].max())
        for zr in (200.0, 400.0, 600.0, 1000.0):
            mm = (np.abs(ex[:, 2] - zr) < 2.0) & (dte(ex[:, 2], ex[:, 0]) < 6.0)
            row[f"skin_peak_rib{zr:g}_control"] = float(v[mm].max())
        rows.append(row)
        print(f"  TE refinement {label} done")
    return rows


def write_csv(path, rows):
    if not rows:
        print(f"  (no rows for {path})"); return
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {path} ({len(rows)} rows)")


def show(rows, cols):
    print("  " + " | ".join(f"{c}" for c in cols))
    for r in rows:
        print("  " + " | ".join((f"{r[c]:.4g}" if isinstance(r.get(c), float) else str(r.get(c))) for c in cols))


def main():
    os.makedirs("results", exist_ok=True)
    print("static jobs:")
    st = static_table()
    print("buckling jobs:")
    bu = buckling_table()
    print("TE refinement:")
    te = te_table()
    write_csv("results/stage3_static_summary.csv", st)
    write_csv("results/stage3_buckling_summary.csv", bu)
    write_csv("results/stage3_te_refinement.csv", te)
    json.dump(dict(static=st, buckling=bu, te_refinement=te, yield_assumed_MPa=YIELD),
              open("results/stage3_summary.json", "w"), indent=1)
    print("wrote results/stage3_summary.json\n")
    show(st, ["mesh", "skin_t_mm", "tip_uy_mm", "tip_vs_B_hi_pct", "skin_p99", "skin_max_dTE_ge3", "skin_max_dTE_ge15", "skin_max", "solid_only_max"])
    show(bu, ["mesh", "skin_t_mm", "critical", "margin_pct", "result"])
    show(te, ["case", "local_median_edge_mm", "skin_peak_rib800", "solid_only_peak_rib800", "ring_1_3mm", "ring_3_6mm", "ring_6_15mm", "ring_15_40mm", "skin_peak_rib600_control"])


if __name__ == "__main__":
    main()
