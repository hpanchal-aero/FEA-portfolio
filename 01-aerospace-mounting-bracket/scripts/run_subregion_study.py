"""
Stage 6 sub-region study -- secondary sample restricted to the
high-t region (t in [4.5, 5.0]) where both coarse spot-checks
(lhs_06, lhs_25) satisfied the displacement and stress constraints,
while every one of the 14 original fixed-scheme (fine-mesh) designs
violated displacement.

Deliberately meshed at the COARSE sizing scheme from the start
(hole_h=0.5, r3_h=0.5, toe_h=0.25 -- Stage 5's L1 values), since the
fine L3 scheme is empirically unsolvable on this hardware for designs
in this size/mass range (9 of 25 original designs failed at up to
10.5GB RAM). This is a SEPARATE, coarser-resolution sub-study, kept
explicitly distinct from the original 14-point fine-mesh dataset --
not merged into one claimed-resolution result.

Design variables: t fixed to [4.5, 5.0] (narrowed); L_gusset and
r_toe retain their full original ranges [15,25] and [1,3], since
their interaction with feasibility in this t-region is not yet known.

Fixed seed=777 (distinct from the original seed=42) for reproducibility.

Runs geometry construction + validity checks + coarse meshing only.
Does NOT solve -- solving is done afterward via run_parametric_sweep.py
per design, same as the main study.
"""

import csv
import os
import sys
from scipy.stats import qmc

sys.path.insert(0, os.path.dirname(__file__))
from build_parametric_bracket import build_bracket_geometry
from mesh_stage6_design import build_mesh

SEED = 777
N_SAMPLES = 9

T_RANGE = (4.5, 5.0)
L_GUSSET_RANGE = (15.0, 25.0)
R_TOE_RANGE = (1.0, 3.0)

# Coarse sizing -- same as the lhs_06/lhs_25 spot-checks
HOLE_H = 0.5
R3_H = 0.5
TOE_H = 0.25

SUBREGION_CSV = "results/stage6_parametric/lhs_subregion_points.csv"
GEOM_CSV = "results/stage6_parametric/geometry_validation_results.csv"
BREP_DIR = "mesh/stage6_parametric"
INP_DIR = "simulation/stage6_parametric"

GEOM_FIELDNAMES = ["design_id", "t", "L_gusset", "r_toe", "effective_leg",
                    "pass", "mass", "failure_stage", "failure_reason"]


def generate_samples():
    sampler = qmc.LatinHypercube(d=3, seed=SEED)
    unit = sampler.random(n=N_SAMPLES)
    lower = [T_RANGE[0], L_GUSSET_RANGE[0], R_TOE_RANGE[0]]
    upper = [T_RANGE[1], L_GUSSET_RANGE[1], R_TOE_RANGE[1]]
    scaled = qmc.scale(unit, lower, upper)

    os.makedirs(os.path.dirname(SUBREGION_CSV), exist_ok=True)
    with open(SUBREGION_CSV, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["design_id", "t", "L_gusset", "r_toe"])
        points = []
        for i, (t, L_gusset, r_toe) in enumerate(scaled, start=1):
            design_id = f"sub_{i:02d}"
            writer.writerow([design_id, f"{t:.6f}", f"{L_gusset:.6f}", f"{r_toe:.6f}"])
            points.append((design_id, t, L_gusset, r_toe))
    print(f"Generated {N_SAMPLES} sub-region LHS points (seed={SEED}), t in {T_RANGE}")
    print(f"Saved: {SUBREGION_CSV}")
    return points


def append_geometry_results(new_rows):
    existing = {}
    if os.path.exists(GEOM_CSV):
        with open(GEOM_CSV, newline="") as f:
            for r in csv.DictReader(f):
                existing[r["design_id"]] = r
    for r in new_rows:
        existing[r["design_id"]] = r
    with open(GEOM_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=GEOM_FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        for design_id in sorted(existing.keys()):
            writer.writerow(existing[design_id])
    print(f"Appended {len(new_rows)} row(s) to {GEOM_CSV}")


def main():
    points = generate_samples()

    geom_rows = []
    mesh_summary = []

    for design_id, t, L_gusset, r_toe in points:
        print(f"\n{'#'*70}\n# {design_id}: t={t:.4f}, L_gusset={L_gusset:.4f}, r_toe={r_toe:.4f}\n{'#'*70}")

        success, brep_path, mass, diag = build_bracket_geometry(t, L_gusset, r_toe, design_id, BREP_DIR)

        row = {
            "design_id": design_id, "t": t, "L_gusset": L_gusset, "r_toe": r_toe,
            "effective_leg": L_gusset - 3.0, "pass": success, "mass": mass,
            "failure_stage": diag.get("failure_stage", ""),
            "failure_reason": diag.get("failure_reason", "; ".join(diag.get("prevalidity_reasons", []))),
        }
        geom_rows.append(row)

        if not success:
            print(f"[{design_id}] Geometry FAILED -- skipping mesh step.")
            mesh_summary.append((design_id, "geometry_failed", None, None))
            continue

        inp_out = os.path.join(INP_DIR, f"stage6_{design_id}_mesh.inp")
        mesh_diag = build_mesh(design_id, t, L_gusset, r_toe, HOLE_H, R3_H, TOE_H, brep_path, inp_out)

        if mesh_diag["success"]:
            mesh_summary.append((design_id, "ok", mesh_diag["n_nodes"], mesh_diag["n_elements"]))
        else:
            mesh_summary.append((design_id, f"mesh_failed:{mesh_diag.get('failure_stage')}", None, None))

    append_geometry_results(geom_rows)

    print(f"\n{'='*70}\nSUB-REGION SUMMARY\n{'='*70}")
    print(f"{'design_id':<10}{'status':<25}{'nodes':<10}{'elements'}")
    for design_id, status, nodes, elems in mesh_summary:
        print(f"{design_id:<10}{status:<25}{str(nodes):<10}{str(elems)}")

    n_ok = sum(1 for _, s, _, _ in mesh_summary if s == "ok")
    print(f"\n{n_ok}/{len(mesh_summary)} designs passed geometry + coarse meshing.")
    print("NO SOLVING performed. Solve individually via run_parametric_sweep.py --design-id sub_XX")


if __name__ == "__main__":
    main()
