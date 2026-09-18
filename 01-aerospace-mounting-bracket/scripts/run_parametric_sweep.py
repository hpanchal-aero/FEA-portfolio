"""
Stage 6 -- solve batch runner for the 22 mesh-feasible parametric designs.

v2: added after a full-batch run crashed the WSL host (11GB available,
no .wslconfig cap). Changes:
  - Hard per-process memory limit (RLIMIT_AS) so an over-budget ccx
    solve is killed cleanly by the OS instead of exhausting host RAM.
  - Incremental checkpointing: solve_results.csv is rewritten after
    EVERY design, and any design with an existing, successfully-
    parsed .frd is skipped (resumable across crashes/sessions).
  - Designs solved smallest-first (by mesh_quality_results.csv element
    count), so a crash on a large design doesn't cost already-viable
    smaller results.

Mesh sizing itself is UNCHANGED -- this only changes execution safety
and ordering, not the fixed L3 sizing scheme.

USAGE:
    python3 run_parametric_sweep.py --dry-run
    python3 run_parametric_sweep.py --design-id baseline
    python3 run_parametric_sweep.py --mem-limit-gb 7
    python3 run_parametric_sweep.py                     # full batch, resumable
"""

import argparse
import csv
import os
import threading
import subprocess
import sys
import time

DEFAULT_CCX_BIN = "/home/harsh/miniconda3/envs/ccx223/bin/ccx"
SIM_DIR = "simulation/stage6_parametric"
GEOM_CSV = "results/stage6_parametric/geometry_validation_results.csv"
MESH_CSV = "results/stage6_parametric/mesh_quality_results.csv"
OUT_CSV = "results/stage6_parametric/solve_results.csv"
LOG_FILE = "results/stage6_parametric/solve_batch.log"

APPLIED_LOAD_X = 235.44
EQUIL_ERROR_FLAG_PCT = 1.0

FIELDNAMES = [
    "design_id", "t", "L_gusset", "r_toe", "mass",
    "solver_status", "return_code", "solve_time_sec",
    "output_file_status", "parser_status", "equilibrium_status",
    "equil_error_pct", "max_ux_mm", "peak_vm_mpa",
    "failure_stage", "failure_reason",
]


def watch_rss(pid, mem_limit_bytes, killed_flag, poll_interval=0.2):
    """Background watchdog: polls actual RSS via /proc, SIGKILLs the
    process if physical memory (not virtual reservation) exceeds the
    cap. Runs until the process exits or is killed."""
    status_path = f"/proc/{pid}/status"
    while True:
        try:
            with open(status_path) as f:
                for line in f:
                    if line.startswith("VmRSS:"):
                        rss_kb = int(line.split()[1])
                        rss_bytes = rss_kb * 1024
                        if rss_bytes > mem_limit_bytes:
                            os.kill(pid, 9)
                            killed_flag["killed"] = True
                            killed_flag["rss_at_kill_gb"] = rss_bytes / (1024 ** 3)
                            return
                        break
        except (FileNotFoundError, ProcessLookupError):
            return  # process already exited
        time.sleep(poll_interval)


def parse_frd_disp_forc(frd_path):
    disp, rf = {}, {}
    with open(frd_path) as f:
        lines = f.readlines()
    mode = None
    for line in lines:
        if line.startswith(" -4") and "DISP" in line:
            mode = "disp"; continue
        if line.startswith(" -4") and "FORC" in line:
            mode = "forc"; continue
        if mode and line.startswith(" -3"):
            mode = None; continue
        if mode == "disp" and line.startswith(" -1"):
            node = int(line[3:13])
            disp[node] = [float(line[13:25]), float(line[25:37]), float(line[37:49])]
        if mode == "forc" and line.startswith(" -1"):
            node = int(line[3:13])
            rf[node] = [float(line[13:25]), float(line[25:37]), float(line[37:49])]
    return disp, rf


def parse_frd_stress(frd_path):
    stress = {}
    with open(frd_path) as f:
        lines = f.readlines()
    mode = False
    for line in lines:
        if line.startswith(" -4") and "STRESS" in line:
            mode = True; continue
        if mode and line.startswith(" -1"):
            node = int(line[3:13])
            stress[node] = [
                float(line[13:25]), float(line[25:37]), float(line[37:49]),
                float(line[49:61]), float(line[61:73]), float(line[73:85]),
            ]
            continue
        if mode and not line.startswith(" -1") and not line.startswith(" -5"):
            mode = False
    return stress


def von_mises(s):
    s11, s22, s33, s12, s13, s23 = s
    return (0.5 * ((s11 - s22) ** 2 + (s22 - s33) ** 2 + (s33 - s11) ** 2
                   + 6 * (s12 ** 2 + s13 ** 2 + s23 ** 2))) ** 0.5


def get_nset_nodes(inp_path, nset_name):
    nodes = []
    capture = False
    with open(inp_path) as f:
        for line in f:
            s = line.strip()
            if s.upper().startswith(f"*NSET, NSET={nset_name}"):
                capture = True; continue
            if capture:
                if s.startswith("*"):
                    break
                nodes.extend(int(x) for x in s.split(",") if x.strip())
    return nodes


def load_geometry_data():
    data = {}
    with open(GEOM_CSV, newline="") as f:
        for row in csv.DictReader(f):
            data[row["design_id"]] = row
    return data


def load_mesh_element_counts():
    counts = {}
    if not os.path.exists(MESH_CSV):
        return counts
    with open(MESH_CSV, newline="") as f:
        for row in csv.DictReader(f):
            if row.get("n_elements"):
                try:
                    counts[row["design_id"]] = int(row["n_elements"])
                except ValueError:
                    pass
    return counts


def discover_designs():
    designs = []
    for fname in sorted(os.listdir(SIM_DIR)):
        if fname.startswith("stage6_") and fname.endswith("_mesh.inp"):
            design_id = fname[len("stage6_"):-len("_mesh.inp")]
            designs.append(design_id)
    return designs


def load_existing_results():
    existing = {}
    if os.path.exists(OUT_CSV):
        with open(OUT_CSV, newline="") as f:
            for r in csv.DictReader(f):
                existing[r["design_id"]] = r
    return existing


def write_results_csv(existing):
    with open(OUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        for design_id in sorted(existing.keys()):
            writer.writerow(existing[design_id])


def already_done(design_id, existing):
    """A design counts as done if it previously solved AND parsed OK."""
    r = existing.get(design_id)
    if not r:
        return False
    return r.get("solver_status") == "ok" and r.get("parser_status") == "ok"


def solve_one(design_id, ccx_bin, timeout, mem_limit_bytes, geom_data, log_f):
    row = {"design_id": design_id}
    inp_name = f"stage6_{design_id}_mesh"
    inp_path = os.path.join(SIM_DIR, inp_name + ".inp")
    frd_path = os.path.join(SIM_DIR, inp_name + ".frd")

    geom_row = geom_data.get(design_id, {})
    row["t"] = geom_row.get("t", "")
    row["L_gusset"] = geom_row.get("L_gusset", "")
    row["r_toe"] = geom_row.get("r_toe", "")
    row["mass"] = geom_row.get("mass", "")

    for k in ["solver_status", "return_code", "solve_time_sec", "output_file_status",
              "parser_status", "equilibrium_status", "equil_error_pct",
              "max_ux_mm", "peak_vm_mpa", "failure_stage", "failure_reason"]:
        row[k] = ""

    if not os.path.exists(inp_path):
        row["solver_status"] = "skipped"
        row["failure_stage"] = "inp_missing"
        row["failure_reason"] = f"{inp_path} not found"
        log_f.write(f"[{design_id}] SKIPPED: {inp_path} not found\n")
        return row

    print(f"\n{'='*70}\n[{design_id}] Solving {inp_name} "
          f"(mem cap {mem_limit_bytes/1e9:.1f}GB) ...\n{'='*70}")
    log_f.write(f"\n[{design_id}] Solving {inp_name} (mem cap {mem_limit_bytes/1e9:.1f}GB) ...\n")

    t0 = time.time()
    t0 = time.time()
    killed_flag = {"killed": False, "rss_at_kill_gb": None}
    try:
        proc = subprocess.Popen(
            [ccx_bin, "-i", inp_name],
            cwd=SIM_DIR, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        )
        watchdog = threading.Thread(
            target=watch_rss, args=(proc.pid, mem_limit_bytes, killed_flag), daemon=True
        )
        watchdog.start()

        stdout, stderr = proc.communicate(timeout=timeout)
        elapsed = time.time() - t0
        row["solve_time_sec"] = round(elapsed, 2)
        row["return_code"] = proc.returncode
        log_f.write(stdout)
        log_f.write(stderr)

        if killed_flag["killed"]:
            row["solver_status"] = "oom_killed"
            row["failure_stage"] = "solve"
            row["failure_reason"] = (
                f"killed by watchdog: RSS reached {killed_flag['rss_at_kill_gb']:.2f}GB "
                f"(cap {mem_limit_bytes/1e9:.1f}GB)"
            )
            print(f"[{design_id}] KILLED by watchdog at "
                  f"{killed_flag['rss_at_kill_gb']:.2f}GB RSS (cap {mem_limit_bytes/1e9:.1f}GB)")
        elif proc.returncode == 0:
            row["solver_status"] = "ok"
        else:
            row["solver_status"] = "nonzero_exit"
            row["failure_stage"] = "solve"
            row["failure_reason"] = f"ccx exited with code {proc.returncode}"
            print(f"[{design_id}] ccx exited with code {proc.returncode}")

    except subprocess.TimeoutExpired:
        proc.kill()
        row["solve_time_sec"] = round(time.time() - t0, 2)
        row["solver_status"] = "timeout"
        row["failure_stage"] = "solve"
        row["failure_reason"] = f"exceeded {timeout}s timeout"
        log_f.write(f"[{design_id}] TIMEOUT after {timeout}s\n")
        print(f"[{design_id}] TIMEOUT after {timeout}s")
        return row

    except Exception as e:
        row["solve_time_sec"] = round(time.time() - t0, 2)
        row["solver_status"] = "exception"
        row["failure_stage"] = "solve"
        row["failure_reason"] = str(e)
        log_f.write(f"[{design_id}] EXCEPTION: {e}\n")
        print(f"[{design_id}] EXCEPTION: {e}")
        return row

    if row["solver_status"] != "ok":
        if not os.path.exists(frd_path):
            row["output_file_status"] = "missing"
            return row

    if not os.path.exists(frd_path):
        row["output_file_status"] = "missing"
        row["failure_stage"] = row["failure_stage"] or "output"
        row["failure_reason"] = row["failure_reason"] or f"{frd_path} not produced"
        print(f"[{design_id}] .frd not produced")
        return row
    row["output_file_status"] = "found"

    try:
        disp, rf = parse_frd_disp_forc(frd_path)
        stress = parse_frd_stress(frd_path)
        root_nodes = get_nset_nodes(inp_path, "NROOT")
        tip_nodes = get_nset_nodes(inp_path, "NTIP")

        if not disp or not rf or not stress:
            raise ValueError(f"empty parse (disp={len(disp)}, rf={len(rf)}, stress={len(stress)})")

        total_rf_x = sum(rf[n][0] for n in root_nodes if n in rf)
        equil_error_pct = abs(abs(total_rf_x) - APPLIED_LOAD_X) / APPLIED_LOAD_X * 100
        row["equil_error_pct"] = round(equil_error_pct, 6)
        row["equilibrium_status"] = "ok" if equil_error_pct < EQUIL_ERROR_FLAG_PCT else "large_error"

        ux_vals = [disp[n][0] for n in tip_nodes if n in disp]
        row["max_ux_mm"] = round(min(ux_vals), 6) if ux_vals else ""

        vm_vals = [von_mises(s) for s in stress.values()]
        row["peak_vm_mpa"] = round(max(vm_vals), 4) if vm_vals else ""

        row["parser_status"] = "ok"
        print(f"[{design_id}] Equil err={equil_error_pct:.4f}%, "
              f"MaxUx={row['max_ux_mm']}mm, PeakVM={row['peak_vm_mpa']}MPa")

    except Exception as e:
        row["parser_status"] = "parse_error"
        row["equilibrium_status"] = "na"
        row["failure_stage"] = row["failure_stage"] or "parse"
        row["failure_reason"] = row["failure_reason"] or str(e)
        print(f"[{design_id}] PARSE ERROR: {e}")
        log_f.write(f"[{design_id}] PARSE ERROR: {e}\n")

    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ccx-bin", default=DEFAULT_CCX_BIN)
    ap.add_argument("--timeout", type=int, default=1800)
    ap.add_argument("--mem-limit-gb", type=float, default=8.0,
                     help="Hard RLIMIT_AS cap per ccx process, in GB")
    ap.add_argument("--design-id", default=None, help="Solve only this one design (test mode)")
    ap.add_argument("--dry-run", action="store_true", help="List planned runs, solve nothing")
    ap.add_argument("--force", action="store_true",
                     help="Re-solve designs even if already successfully solved+parsed")
    args = ap.parse_args()

    if not os.path.exists(args.ccx_bin):
        print(f"ERROR: ccx binary not found at {args.ccx_bin}")
        sys.exit(1)

    mem_limit_bytes = int(args.mem_limit_gb * (1024 ** 3))
    geom_data = load_geometry_data()
    elem_counts = load_mesh_element_counts()
    all_designs = discover_designs()
    existing = load_existing_results()

    if args.design_id:
        if args.design_id not in all_designs:
            print(f"ERROR: design_id '{args.design_id}' has no .inp in {SIM_DIR}")
            sys.exit(1)
        designs = [args.design_id]
    else:
        designs = all_designs
        # smallest-first ordering; unknown element counts sort last
        designs = sorted(designs, key=lambda d: elem_counts.get(d, float("inf")))

    if not args.force:
        skipped_done = [d for d in designs if already_done(d, existing)]
        designs = [d for d in designs if not already_done(d, existing)]
        if skipped_done:
            print(f"Skipping {len(skipped_done)} already-solved design(s): {skipped_done}")

    print(f"Discovered {len(all_designs)} mesh .inp file(s) in {SIM_DIR}")
    print(f"Planned solves this run: {len(designs)} -> {designs}")
    print(f"ccx binary: {args.ccx_bin}")
    print(f"Timeout per design: {args.timeout}s")
    print(f"Memory cap per design: {args.mem_limit_gb}GB")

    if args.dry_run:
        print("\nDRY RUN -- no solves performed.")
        return

    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    with open(LOG_FILE, "a") as log_f:
        log_f.write(f"\n{'#'*70}\n# Batch run: {len(designs)} design(s), "
                     f"mem_cap={args.mem_limit_gb}GB, "
                     f"{time.strftime('%Y-%m-%d %H:%M:%S')}\n{'#'*70}\n")
        for design_id in designs:
            row = solve_one(design_id, args.ccx_bin, args.timeout, mem_limit_bytes,
                             geom_data, log_f)
            existing[design_id] = row
            write_results_csv(existing)  # checkpoint after EVERY design

    n_ok = sum(1 for d in designs
               if existing.get(d, {}).get("solver_status") == "ok"
               and existing.get(d, {}).get("parser_status") == "ok")
    print(f"\n{'='*70}\nBATCH SUMMARY (this run): {n_ok}/{len(designs)} fully OK "
          f"(solved + parsed)\n{'='*70}")
    for d in designs:
        r = existing.get(d, {})
        print(f"  {d:<10} solver={r.get('solver_status'):<12} "
              f"parser={r.get('parser_status'):<12} equil={r.get('equilibrium_status'):<12} "
              f"{r.get('failure_reason','')}")
    print(f"\nSaved: {OUT_CSV}")
    print(f"Log: {LOG_FILE}")


if __name__ == "__main__":
    main()
