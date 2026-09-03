"""
Project 01 - Aerospace Mounting Bracket
Collect results from the 5 already-completed structured hex
convergence cases directly from their saved files - no re-solving.

Each case's mesh dimensions (dx, dy, dz) are computed from the
locked nx/ny/nz table (see chat log). Tip Uz and root-band max von
Mises are extracted from analysis.frd using the same validated
fixed-width parsing as hex_lib.py's extract_results().
"""

import csv
import os
import sys

sys.path.insert(0, ".")
import hex_lib as lib

SIM_ROOT = os.path.abspath("../simulation/convergence_hex")
RESULTS_CSV = os.path.abspath("../results/convergence_study_hex.csv")

# Locked table (see chat log) - label, nx, ny, nz, case directory name
LEVELS = [
    {"label": "2.0mm", "nx": 30,  "ny": 20,  "nz": 2,  "dir": "mesh_2.0mm"},
    {"label": "1.3mm", "nx": 46,  "ny": 31,  "nz": 3,  "dir": "mesh_1.3mm"},
    {"label": "0.8mm", "nx": 75,  "ny": 50,  "nz": 5,  "dir": "mesh_0.8mm"},
    {"label": "0.5mm", "nx": 120, "ny": 80,  "nz": 8,  "dir": "mesh_0.5mm"},
    {"label": "0.3mm", "nx": 200, "ny": 134, "nz": 14, "dir": "mesh_0.3mm"},
]

results = []

for level in LEVELS:
    case_dir = os.path.join(SIM_ROOT, level["dir"])
    print(f"\n=== {level['label']} ({level['dir']}) ===")

    clean_inp = os.path.join(case_dir, "mesh_clean.inp")
    nset_inp = os.path.join(case_dir, "nsets.inp")

    nodes, elements = lib.parse_clean_mesh(clean_inp)
    n_elements = len(elements)
    print(f"  Elements: {n_elements}")

    tip_uz_avg, max_vm_root = lib.extract_results(case_dir, nset_inp, nodes)
    print(f"  Tip Uz avg: {tip_uz_avg:.6f} mm")
    print(f"  Max von Mises (root band): {max_vm_root:.3f} MPa")

    # Try to recover solve time from run_log.txt if present
    solve_time = None
    log_path = os.path.join(case_dir, "run_log.txt")
    if os.path.exists(log_path):
        with open(log_path) as f:
            content = f.read()
        import re
        m = re.search(r"Total CalculiX Time:\s*([\d.]+)", content)
        if m:
            solve_time = float(m.group(1))

    nx, ny, nz = level["nx"], level["ny"], level["nz"]
    dx = round(lib.L / (nx - 1), 4)
    dy = round(lib.B / (ny - 1), 4)
    dz = round(lib.H / (nz - 1), 4)

    results.append({
        "label": level["label"],
        "nx": nx, "ny": ny, "nz": nz,
        "dx_mm": dx, "dy_mm": dy, "dz_mm": dz,
        "n_elements": n_elements,
        "tip_uz_avg_mm": tip_uz_avg,
        "max_vm_root_mpa": max_vm_root,
        "solve_time_s": solve_time,
    })

# Percent change vs previous level (ordered coarse -> fine, as listed)
for i in range(1, len(results)):
    prev, curr = results[i-1], results[i]
    duz = abs(curr["tip_uz_avg_mm"] - prev["tip_uz_avg_mm"]) / \
        abs(prev["tip_uz_avg_mm"]) * 100
    dvm = abs(curr["max_vm_root_mpa"] - prev["max_vm_root_mpa"]) / \
        abs(prev["max_vm_root_mpa"]) * 100
    curr["pct_change_uz"] = round(duz, 4)
    curr["pct_change_vm"] = round(dvm, 4)
results[0]["pct_change_uz"] = ""
results[0]["pct_change_vm"] = ""

os.makedirs(os.path.dirname(RESULTS_CSV), exist_ok=True)
fieldnames = ["label", "nx", "ny", "nz", "dx_mm", "dy_mm", "dz_mm",
              "n_elements", "tip_uz_avg_mm", "max_vm_root_mpa",
              "pct_change_uz", "pct_change_vm", "solve_time_s"]
with open(RESULTS_CSV, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(results)

print(f"\nWrote: {RESULTS_CSV}")
print("\n=== Summary ===")
for r in results:
    pct = f"{r['pct_change_uz']}%" if r['pct_change_uz'] != "" else "—"
    print(f"{r['label']:>6}: {r['n_elements']:>7} elem  "
          f"Uz={r['tip_uz_avg_mm']:.6f}mm  VM={r['max_vm_root_mpa']:.3f}MPa  "
          f"Δ={pct}")
