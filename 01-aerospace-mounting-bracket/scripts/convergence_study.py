"""
Project 01 - Aerospace Mounting Bracket
Mesh convergence study driver.

Modes:
  --validate   Run ONE case (1.3 mm, same size as the hand-verified
               baseline) through the full automated pipeline and
               compare against the manually confirmed reference
               values. Must pass before --sweep is trusted.
  --sweep      Run all 5 refinement levels and write
               results/convergence_study.csv
"""

import argparse
import csv
import os
import sys

import convergence_lib as lib

SIM_ROOT = os.path.abspath("../simulation/convergence")
RESULTS_CSV = os.path.abspath("../results/convergence_study.csv")

# Manually confirmed baseline reference values (see chat log,
# hand-verified single case at 1.3 mm element size)
BASELINE_TIP_UZ = -1.048931   # mm
BASELINE_TOL_PCT = 5.0        # allow for Delaunay non-determinism


def run_case(elem_size, case_name):
    case_dir = os.path.join(SIM_ROOT, case_name)
    os.makedirs(case_dir, exist_ok=True)

    print(f"\n=== Case {case_name}: element size {elem_size} mm ===")

    geo_path = lib.write_geo(case_dir, elem_size)
    raw_inp, nset_inp, n_elements = lib.mesh_and_export(case_dir, geo_path)
    print(f"  Meshed: {n_elements} C3D10 elements")

    clean_inp = lib.clean_mesh(case_dir, raw_inp)
    nodes, elements = lib.parse_clean_mesh(clean_inp)

    _, n_faces = lib.build_facial_surface(case_dir, nodes, elements)
    print(f"  Tip facial surface: {n_faces} element faces")
    if n_faces < 50:
        raise RuntimeError(
            f"Facial surface has only {n_faces} faces - almost certainly "
            f"broken (expect ~100-1000 depending on mesh density). "
            f"Stopping before running CalculiX on a bad surface."
        )

    analysis_inp = lib.write_analysis_deck(case_dir)
    solve_time = lib.run_ccx(case_dir, analysis_inp)
    print(f"  Solved in {solve_time:.2f} s")

    tip_uz_avg, max_vm_root = lib.extract_results(case_dir, nset_inp, nodes)
    print(f"  Tip Uz avg: {tip_uz_avg:.6f} mm")
    print(f"  Max von Mises in root band [{lib.ROOT_X_MIN},{lib.ROOT_X_MAX}] mm: "
          f"{max_vm_root:.3f} MPa")

    return {
        "elem_size_mm": elem_size,
        "n_elements": n_elements,
        "tip_uz_avg_mm": tip_uz_avg,
        "max_vm_root_mpa": max_vm_root,
        "solve_time_s": solve_time,
    }


def validate():
    print("Running validation case (1.3 mm) through automated pipeline...")
    result = run_case(1.3, "validation_1.3mm")

    pct_diff = abs(result["tip_uz_avg_mm"] - BASELINE_TIP_UZ) / abs(BASELINE_TIP_UZ) * 100

    print("\n--- Validation against hand-verified baseline ---")
    print(f"  Hand-verified tip Uz avg: {BASELINE_TIP_UZ:.6f} mm")
    print(f"  Automated pipeline tip Uz avg: {result['tip_uz_avg_mm']:.6f} mm")
    print(f"  Percent difference: {pct_diff:.3f}%  (tolerance: {BASELINE_TOL_PCT}%)")

    if pct_diff > BASELINE_TOL_PCT:
        print("\nVALIDATION FAILED - automated pipeline does not reproduce "
              "the hand-verified baseline within tolerance. DO NOT proceed "
              "to --sweep. Investigate before continuing.")
        sys.exit(1)
    else:
        print("\nVALIDATION PASSED - automated pipeline reproduces the "
              "hand-verified baseline within tolerance. Safe to run --sweep.")


def sweep():
    sizes = [2.0, 1.3, 0.8, 0.5, 0.3]
    results = []
    for size in sizes:
        case_name = f"mesh_{size}mm"
        result = run_case(size, case_name)
        results.append(result)

        if len(results) >= 2:
            prev, curr = results[-2], results[-1]
            duz = abs(curr["tip_uz_avg_mm"] - prev["tip_uz_avg_mm"]) / \
                abs(prev["tip_uz_avg_mm"]) * 100
            dvm = abs(curr["max_vm_root_mpa"] - prev["max_vm_root_mpa"]) / \
                abs(prev["max_vm_root_mpa"]) * 100
            print(f"  Change vs previous level: Uz {duz:.2f}%, "
                  f"root max VM {dvm:.2f}%")

    os.makedirs(os.path.dirname(RESULTS_CSV), exist_ok=True)
    with open(RESULTS_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))
        writer.writeheader()
        writer.writerows(results)

    print(f"\nWrote: {RESULTS_CSV}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--validate", action="store_true")
    group.add_argument("--sweep", action="store_true")
    args = parser.parse_args()

    if args.validate:
        validate()
    elif args.sweep:
        sweep()
