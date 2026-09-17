"""
Stage 6 -- geometry-only validation batch runner.

Reads the committed LHS design-point CSV and runs
build_bracket_geometry() (from build_parametric_bracket.py) on every
point, WITHOUT meshing or solving. Produces a pass/fail summary table
plus a results CSV recording mass, effective leg length, and (for any
failures) the exact failure stage and reason.

Per the Stage 6 rule: a design point that fails here is reported, not
silently dropped or repaired. No optimization or trade-space claim
may be made until this summary shows every sampled point either
passing or being explicitly excluded with a documented reason.
"""

import csv
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from build_parametric_bracket import build_bracket_geometry

INPUT_CSV = "results/stage6_parametric/lhs_design_points.csv"
OUTPUT_CSV = "results/stage6_parametric/geometry_validation_results.csv"
BREP_OUT_DIR = "mesh/stage6_parametric"


def main():
    if not os.path.exists(INPUT_CSV):
        print(f"ERROR: {INPUT_CSV} not found. Run generate_lhs_samples.py first.")
        sys.exit(1)

    with open(INPUT_CSV, newline="") as f:
        reader = csv.DictReader(f)
        design_points = [
            (row["design_id"], float(row["t"]), float(row["L_gusset"]), float(row["r_toe"]))
            for row in reader
        ]

    print(f"Loaded {len(design_points)} design points from {INPUT_CSV}\n")

    results = []
    for design_id, t, L_gusset, r_toe in design_points:
        success, brep_path, mass, diag = build_bracket_geometry(
            t, L_gusset, r_toe, design_id, BREP_OUT_DIR
        )
        results.append({
            "design_id": design_id,
            "t": t,
            "L_gusset": L_gusset,
            "r_toe": r_toe,
            "effective_leg": L_gusset - 3.0,  # R3 fixed at 3.0
            "pass": success,
            "mass": mass if success else "",
            "failure_stage": diag.get("failure_stage", ""),
            "failure_reason": diag.get("failure_reason", "; ".join(diag.get("prevalidity_reasons", []))),
        })

    os.makedirs(os.path.dirname(OUTPUT_CSV), exist_ok=True)
    with open(OUTPUT_CSV, "w", newline="") as f:
        fieldnames = ["design_id", "t", "L_gusset", "r_toe", "effective_leg",
                      "pass", "mass", "failure_stage", "failure_reason"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    n_pass = sum(1 for r in results if r["pass"])
    n_fail = len(results) - n_pass

    print(f"\n{'='*90}")
    print(f"GEOMETRY VALIDATION SUMMARY: {n_pass}/{len(results)} PASSED, {n_fail}/{len(results)} FAILED")
    print(f"{'='*90}")
    print(f"{'design_id':<10}{'t':<8}{'L_gusset':<10}{'r_toe':<8}{'eff_leg':<9}{'PASS':<7}{'mass':<14}{'failure'}")
    print("-" * 90)
    for r in results:
        mass_str = f"{r['mass']:.2f}" if r["pass"] else "-"
        pass_str = "PASS" if r["pass"] else "FAIL"
        failure_str = f"{r['failure_stage']}: {r['failure_reason']}" if not r["pass"] else ""
        print(f"{r['design_id']:<10}{r['t']:<8.3f}{r['L_gusset']:<10.3f}{r['r_toe']:<8.3f}"
              f"{r['effective_leg']:<9.3f}{pass_str:<7}{mass_str:<14}{failure_str}")

    print(f"\nSaved: {OUTPUT_CSV}")

    if n_fail > 0:
        print(f"\n{n_fail} design point(s) FAILED geometry validation.")
        print("These must be resolved or explicitly excluded before any")
        print("optimization or trade-space claim is made, per the Stage 6 rule.")


if __name__ == "__main__":
    main()
