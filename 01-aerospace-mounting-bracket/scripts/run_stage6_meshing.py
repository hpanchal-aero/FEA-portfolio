"""
Stage 6 -- batch mesh generation across all 25 LHS design points.

Uses Stage 5's accepted L3 sizing values (hole_h=0.25, r3_h=0.25,
toe_h=0.10) as a FIXED sizing scheme applied to every design, per the
Stage 6 verification strategy (Option B): mesh sizing converged once
at the baseline (= Stage 5's geometry, confirmed exact reproduction),
held constant across the sweep. This establishes the SIZING SCHEME as
a controlled input -- it does NOT guarantee equivalent mesh quality
outcomes across all 25 (different) geometries. Every design's actual
element count, quality, and solver feasibility is logged individually
and must be inspected before any solving occurs.

NO CALCULIX INVOCATION IN THIS SCRIPT. Mesh generation and quality
logging only, per the explicit two-gate instruction: no solving until
the first mesh passes quality checks and the full batch is reviewed.
"""

import csv
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from mesh_stage6_design import build_mesh

INPUT_CSV = "results/stage6_parametric/lhs_design_points.csv"
OUTPUT_CSV = "results/stage6_parametric/mesh_quality_results.csv"
BREP_DIR = "mesh/stage6_parametric"
INP_DIR = "simulation/stage6_parametric"

# Fixed sizing scheme, reused directly from Stage 5's accepted L3 mesh
HOLE_H = 0.25
R3_H = 0.25
TOE_H = 0.10


def main():
    if not os.path.exists(INPUT_CSV):
        print(f"ERROR: {INPUT_CSV} not found.")
        sys.exit(1)

    with open(INPUT_CSV, newline="") as f:
        reader = csv.DictReader(f)
        design_points = [
            (row["design_id"], float(row["t"]), float(row["L_gusset"]), float(row["r_toe"]))
            for row in reader
        ]

    print(f"Loaded {len(design_points)} design points from {INPUT_CSV}")
    print(f"Fixed sizing scheme: hole={HOLE_H}, R3={R3_H}, toe={TOE_H} "
          f"(Stage 5 L3 values, held constant across the sweep)\n")

    results = []
    for design_id, t, L_gusset, r_toe in design_points:
        brep_path = os.path.join(BREP_DIR, f"stage6_design_{design_id}.brep")
        inp_out = os.path.join(INP_DIR, f"stage6_{design_id}_mesh.inp")

        diag = build_mesh(design_id, t, L_gusset, r_toe, HOLE_H, R3_H, TOE_H,
                           brep_path, inp_out)
        results.append(diag)

    os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)
    fieldnames = [
        "design_id", "t", "L_gusset", "r_toe", "hole_h", "r3_h", "toe_h",
        "success", "n_nodes", "n_elements", "n_linear_tets",
        "min_linear_volume", "mean_linear_volume", "n_bad_volume_elements",
        "n_root_nodes", "n_tip_nodes", "est_equations", "pct_of_ceiling",
        "inp_path", "failure_stage", "failure_reason",
    ]
    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for r in results:
            writer.writerow(r)

    n_pass = sum(1 for r in results if r["success"])
    n_fail = len(results) - n_pass

    print(f"\n{'='*100}")
    print(f"MESH GENERATION SUMMARY: {n_pass}/{len(results)} PASSED, {n_fail}/{len(results)} FAILED")
    print(f"{'='*100}")
    print(f"{'design_id':<10}{'nodes':<10}{'elements':<11}{'eq_%ceil':<11}{'bad_vol':<9}{'PASS':<7}{'failure'}")
    print("-" * 100)
    for r in results:
        pass_str = "PASS" if r["success"] else "FAIL"
        nodes_str = str(r.get("n_nodes", "-"))
        elems_str = str(r.get("n_elements", "-"))
        pct_str = f"{r.get('pct_of_ceiling', 0):.1f}" if r.get("pct_of_ceiling") is not None else "-"
        bad_str = str(r.get("n_bad_volume_elements", "-"))
        failure_str = f"{r.get('failure_stage','')}: {r.get('failure_reason','')}" if not r["success"] else ""
        print(f"{r['design_id']:<10}{nodes_str:<10}{elems_str:<11}{pct_str:<11}{bad_str:<9}{pass_str:<7}{failure_str}")

    print(f"\nSaved: {OUTPUT_CSV}")

    if n_fail > 0:
        print(f"\n{n_fail} design point(s) FAILED mesh generation.")
        print("Per the Stage 6 rule, these must be investigated and either resolved")
        print("or explicitly excluded -- NOT silently dropped -- before any solving.")
    else:
        print("\nAll designs passed mesh generation and quality checks.")
        print("NO SOLVING has been performed. Review this summary and the full CSV")
        print("before proceeding to the solve phase.")


if __name__ == "__main__":
    main()
