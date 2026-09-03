"""
Project 01 - Aerospace Mounting Bracket
Structured hex convergence study driver (verification plate ONLY).

This documents VERIFICATION of the numerical methodology (does the
FE model solve the mathematical model correctly, and is the solution
mesh-independent) - NOT validation (does the model represent a real
physical bracket). No physical test data exists for this geometry.

The structured C3D20 hex meshing approach is scoped to this simple
plate. The final ribbed/filleted bracket geometry will require a
separate mesh methodology decision.
"""

import argparse
import csv
import os
import sys

import hex_lib as lib

SIM_ROOT = os.path.abspath("../simulation/convergence_hex")
RESULTS_CSV = os.path.abspath("../results/convergence_study_hex.csv")

# Hand-confirmed standalone reference (see chat log)
REF_TIP_UZ = -1.050090
REF_ROOT_VM = 121.639
REF_TOL_PCT = 2.0

LEVELS = [
    {"label": "2.0mm", "nx": 30,  "ny": 20,  "nz": 2},
    {"label": "1.3mm", "nx": 46,  "ny": 31,  "nz": 3},
    {"label": "0.8mm", "nx": 75,  "ny": 50,  "nz": 5},
    {"label": "0.5mm", "nx": 120, "ny": 80,  "nz": 8},
    {"label": "0.3mm", "nx": 200, "ny": 134, "nz": 14},
]


def run_case(nx, ny, nz, case_name):
    case_dir = os.path.join(SIM_ROOT, case_name)
    os.makedirs(case_dir, exist_ok=True)

    print(f"\n=== Case {case_name}: nx={nx} ny={ny} nz={nz} ===")

    geo_path = lib.write_geo(case_dir, nx, ny, nz)
    raw_inp, nset_inp, n_elements = lib.mesh_and_export(case_dir, geo_path)

    expected_elements = (nx-1) * (ny-1) * (nz-1)
    print(f"  Meshed: {n_elements} C3D20 elements "
          f"(expected {expected_elements})")
    if n_elements != expected_elements:
        raise RuntimeError(
            f"Element count mismatch: got {n_elements}, "
            f"expected {expected_elements} - stopping."
        )

    clean_inp = lib.clean_mesh(case_dir, raw_inp)
    nodes, elements = lib.parse_clean_mesh(clean_inp)

    _, n_faces = lib.build_facial_surface(case_dir, nodes, elements)
    expected_faces = (ny-1) * (nz-1)
    print(f"  Tip facial surface: {n_faces} faces (expected {expected_faces})")
    if n_faces != expected_faces:
        raise RuntimeError(
            f"Facial surface count mismatch: got {n_faces}, "
            f"expected {expected_faces} - stopping."
        )

    max_node_id = max(nodes.keys())
    ref_node_id = max_node_id + 1000
    print(f"  Max real node ID: {max_node_id}, using reference node "
          f"ID: {ref_node_id}")

    analysis_inp = lib.write_analysis_deck(case_dir, ref_node_id)
    solve_time = lib.run_ccx(case_dir, analysis_inp)
    print(f"  Solved in {solve_time:.2f} s")

    tip_uz_avg, max_vm_root = lib.extract_results(case_dir, nset_inp, nodes)
    print(f"  Tip Uz avg: {tip_uz_avg:.6f} mm")
    print(f"  Max von Mises in root band: {max_vm_root:.3f} MPa")

    dx = lib.L / (nx - 1)
    dy = lib.B / (ny - 1)
    dz = lib.H / (nz - 1)

    return {
        "nx": nx, "ny": ny, "nz": nz,
        "dx_mm": round(dx, 4), "dy_mm": round(dy, 4), "dz_mm": round(dz, 4),
        "n_elements": n_elements,
        "tip_uz_avg_mm": tip_uz_avg,
        "max_vm_root_mpa": max_vm_root,
        "solve_time_s": solve_time,
    }


def validate():
    print("Reproducing hand-confirmed standalone 0.8mm case "
          "(nx=75, ny=50, nz=5)...")
    result = run_case(75, 50, 5, "validation_0.8mm")

    pct_uz = abs(result["tip_uz_avg_mm"] - REF_TIP_UZ) / abs(REF_TIP_UZ) * 100
    pct_vm = abs(result["max_vm_root_mpa"] - REF_ROOT_VM) / REF_ROOT_VM * 100

    print("\n--- Validation against hand-confirmed standalone result ---")
    print(f"  Tip Uz avg: hand={REF_TIP_UZ:.6f}  automated="
          f"{result['tip_uz_avg_mm']:.6f}  diff={pct_uz:.3f}%")
    print(f"  Root max VM: hand={REF_ROOT_VM:.3f}  automated="
          f"{result['max_vm_root_mpa']:.3f}  diff={pct_vm:.3f}%")
    print(f"  Tolerance: {REF_TOL_PCT}%")

    if pct_uz > REF_TOL_PCT or pct_vm > REF_TOL_PCT:
        print("\nVALIDATION FAILED - automated pipeline does not reproduce "
              "the hand-confirmed standalone result. DO NOT proceed to "
              "--sweep. Investigate before continuing.")
        sys.exit(1)
    else:
        print("\nVALIDATION PASSED - safe to run --sweep.")


def sweep():
    results = []
    for level in LEVELS:
        case_name = f"mesh_{level['label']}"
        result = run_case(level["nx"], level["ny"], level["nz"], case_name)
        result["label"] = level["label"]
        results.append(result)

        if len(results) >= 2:
            prev, curr = results[-2], results[-1]
            duz = abs(curr["tip_uz_avg_mm"] - prev["tip_uz_avg_mm"]) / \
                abs(prev["tip_uz_avg_mm"]) * 100
            dvm = abs(curr["max_vm_root_mpa"] - prev["max_vm_root_mpa"]) / \
                abs(prev["max_vm_root_mpa"]) * 100
            curr["pct_change_uz"] = round(duz, 4)
            curr["pct_change_vm"] = round(dvm, 4)
            print(f"  Change vs previous level: Uz {duz:.2f}%, "
                  f"root max VM {dvm:.2f}%")
        else:
            result["pct_change_uz"] = ""
            result["pct_change_vm"] = ""

    os.makedirs(os.path.dirname(RESULTS_CSV), exist_ok=True)
    fieldnames = ["label", "nx", "ny", "nz", "dx_mm", "dy_mm", "dz_mm",
                  "n_elements", "tip_uz_avg_mm", "max_vm_root_mpa",
                  "pct_change_uz", "pct_change_vm", "solve_time_s"]
    with open(RESULTS_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
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
