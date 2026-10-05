#!/usr/bin/env python3
"""
Collect the Stage 4 results (distributed skin pressure load) into tables (results/*.csv, .json).

  results/stage4_static_summary.csv     tip deflection, stress statistics (full jobs), per job
  results/stage4_twist.csv              fitted nose-up rotation at the rib planes, per job
  results/stage4_buckling_summary.csv   buckling factors (6 modes) from the .dat files
  results/stage4_load_check.csv         independent re-read of the written *CLOAD blocks
  results/stage4_summary.json           all of the above

Definitions are those of Stage 3 (collect_stage3_results.stress_stats / dte): design skin stress = max VM
over expanded skin nodes with dTE >= 3 mm and z > 10 mm. The sensitivity rows (excluding z >= 995, or
z >= 990 together with dTE < 6) use cutoffs chosen AFTER inspecting the data. Skin stress = CalculiX
expanded-node values (interpretation verified near the root only). Yield 290 MPa is ASSUMED, uncited.
The applied lift is an ANALYTIC PLACEHOLDER (elliptical span x thin-airfoil chordwise, antisymmetric
upper/lower), not an aerodynamic prediction. Aircraft mass 20 kg is an assumption.

Run with the system Python (needs the local simulation/ files):
  /usr/bin/python3 scripts/collect_stage4_results.py
"""
import contextlib
import csv
import io
import json
import os
import re
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collect_stage3_results import stress_stats, dte, YIELD, REQ_BUCKLE, ROW_RE
from extract_stage4_twist import analyse, RIB_Z

JOBS = [
    dict(label="S3 spar-cap load 4/8", mesh="s4_k8", case="stage3_spar_cap", main=True,
         job="simulation/stage3_skin_s4_k8_t1_n3p8"),
    dict(label="S4 skin pressure 4/8", mesh="s4_k8", case="stage4_full", main=True,
         job="simulation/stage4_s4_k8_t1_n3p8"),
    dict(label="S4 lift only 4/8", mesh="s4_k8", case="stage4_lift", main=False,
         job="simulation/stage4_s4_k8_t1_n3p8_lift"),
    dict(label="S4 gravity only 4/8", mesh="s4_k8", case="stage4_grav", main=False,
         job="simulation/stage4_s4_k8_t1_n3p8_grav"),
    dict(label="S3 spar-cap load 4/6", mesh="s4_k6", case="stage3_spar_cap", main=True,
         job="simulation/stage3_skin_s4_k6_t1_n3p8"),
    dict(label="S4 skin pressure 4/6", mesh="s4_k6", case="stage4_full", main=True,
         job="simulation/stage4_s4_k6_t1_n3p8"),
]


def have(job):
    return os.path.isfile(job + ".frd") and os.path.isfile(job + ".inp")


def design_extras(vm, eid, ex):
    v = vm[eid]
    d = dte(ex[:, 2], ex[:, 0])
    des = (d >= 3.0) & (ex[:, 2] > 10.0)
    up = ex[:, 1] > 0.0
    return dict(
        design_max_upper=float(v[des & up].max()),
        design_max_lower=float(v[des & ~up].max()),
        design_max_excl_z_ge_995_posthoc=float(v[des & (ex[:, 2] < 995.0)].max()),
        design_max_excl_z_ge_990_dTE_lt_6_posthoc=float(v[des & (ex[:, 2] < 990.0) & (d >= 6.0)].max()),
    )


def buckle_factors(job):
    # same parsing as collect_stage3_results.buckling_table (copied: that logic has hardcoded jobs)
    lines = open(job + ".dat").readlines()
    hi = next(i for i, l in enumerate(lines) if "BUCKLINGFACTOROUTPUT" in "".join(l.upper().split()))
    fac = []
    for l in lines[hi:hi + 60]:
        m = ROW_RE.match(l)
        if m:
            fac.append(float(m.group(2).replace("D", "E").replace("d", "E")))
        elif fac and l.strip() == "":
            break
    return fac


def load_check(inp):
    ids, xyz, cl = [], [], []
    L = None
    sec = None
    with open(inp) as f:
        for line in f:
            if line.startswith("** LOADCASE"):
                L = float(re.search(r"LIFT_N=([0-9.eE+-]+)", line).group(1))
            if line.startswith("*") and not line.startswith("**"):
                sec = line.strip().upper().split(",")[0]
                continue
            if line.startswith("**") or not line.strip():
                continue
            p = line.split(",")
            if sec == "*NODE":
                ids.append(int(p[0]))
                xyz.append([float(v) for v in p[1:4]])
            elif sec == "*CLOAD":
                cl.append((int(p[0]), int(p[1]), float(p[2])))
    ids = np.array(ids)
    coords = np.zeros((ids.max() + 1, 3))
    coords[ids] = np.array(xyz)
    F = np.zeros((ids.max() + 1, 3))
    for nd, dof, val in cl:
        F[nd, dof - 1] += val
    sumF = F.sum(axis=0)
    z, x = coords[:, 2], coords[:, 0]
    mroot = float((F[:, 1] * z).sum())
    mexp = L * 1000.0 * 4.0 / (3.0 * np.pi)
    q0 = L / (1000.0 * np.pi / 4.0)
    zz = np.linspace(0.0, 1000.0, 200001)
    qq = q0 * np.sqrt(1.0 - (zz / 1000.0) ** 2)
    texp = float(np.trapz(qq * (-0.05 * (300.0 - 0.15 * zz)), zz))
    tfe = float((F[:, 1] * x).sum())
    c = 300.0 - 0.15 * z
    xi = (x + 0.3 * c) / c
    bins = []
    for b0 in range(0, 1000, 100):
        mk = (z >= b0) & (z < b0 + 100)
        sf = F[mk, 1].sum()
        bins.append(float((F[mk, 1] * xi[mk]).sum() / sf))
    return dict(L_header=L, n_cload_lines=len(cl), sum_Fx=float(sumF[0]), sum_Fy=float(sumF[1]),
                sum_Fz=float(sumF[2]), root_moment=mroot, root_moment_expected=float(mexp),
                root_moment_err_pct=100 * (mroot / mexp - 1.0), sum_FyX=tfe, sum_FyX_expected=texp,
                sum_FyX_err_pct=100 * (tfe / texp - 1.0), xi_cp_min=min(bins), xi_cp_max=max(bins))


def write_csv(path, rows):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"  wrote {path} ({len(rows)} rows)")


def main():
    static, twist, buck, loadc = [], [], [], []
    for J in JOBS:
        job = J["job"]
        if not have(job):
            print(f"  skip {job} (missing .frd/.inp)")
            continue
        jf = os.path.join("results", os.path.basename(job) + "_results.json")
        if not os.path.isfile(jf):
            print(f"  skip {job} (missing {jf})")
            continue
        j = json.load(open(jf))
        row = dict(label=J["label"], mesh=J["mesh"], case=J["case"], job=os.path.basename(job),
                   lift_N=j["lift_N"], n=j["n"], skin_t_mm=j["t_skin_mm"],
                   tip_uy_mm=j["uy_spar_axis_mm"], beam_B_hi_mm=j["bracket_hi_mm"],
                   tip_vs_B_hi_pct=100 * (j["uy_spar_axis_mm"] / j["bracket_hi_mm"] - 1.0))
        if J["main"]:
            st, (vm, eid, ex, so, coords) = stress_stats(job)
            row.update(st)
            row.update(design_extras(vm, eid, ex))
            row["FoS_design_skin_assumed_yield"] = YIELD / st["skin_max_dTE_ge3"]
        static.append(row)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            th = analyse(job)
        trow = dict(label=J["label"], mesh=J["mesh"], case=J["case"])
        for zc, t in zip(RIB_Z, th):
            trow[f"theta_z{int(zc)}_rad"] = float(t)
        trow["theta_tip_deg"] = float(np.degrees(th[-1]))
        twist.append(trow)
        if J["main"] and os.path.isfile(job + "_buckle.dat"):
            fac = buckle_factors(job + "_buckle")
            b = dict(label=J["label"], mesh=J["mesh"], case=J["case"])
            for k, f in enumerate(fac, 1):
                b[f"mode{k}"] = f
            b.update(critical=min(fac), requirement=REQ_BUCKLE, margin_pct=100 * (min(fac) / REQ_BUCKLE - 1.0),
                     result="PASS" if min(fac) >= REQ_BUCKLE else "FAIL")
            buck.append(b)
        if J["case"].startswith("stage4"):
            lc = load_check(job + ".inp") if J["case"] != "stage4_grav" else None
            if lc:
                lc.update(label=J["label"], mesh=J["mesh"], case=J["case"])
                loadc.append(lc)
        print(f"  done {J['label']}")

    os.makedirs("results", exist_ok=True)
    write_csv("results/stage4_static_summary.csv", static)
    write_csv("results/stage4_twist.csv", twist)
    write_csv("results/stage4_buckling_summary.csv", buck)
    write_csv("results/stage4_load_check.csv", loadc)
    notes = ("Stage 4 load = analytic placeholder (elliptical span x thin-airfoil chordwise, antisymmetric "
             "upper suction / lower pressure), NOT an aerodynamic prediction. 20 kg aircraft mass and 290 MPa "
             "yield are assumptions. Design-stress sensitivity cutoffs chosen after inspecting the data.")
    json.dump(dict(notes=notes, static=static, twist=twist, buckling=buck, load_check=loadc),
              open("results/stage4_summary.json", "w"), indent=1)
    print("  wrote results/stage4_summary.json")

    print("\nStage 4 vs Stage 3 (same mesh):")
    for mesh in ("s4_k8", "s4_k6"):
        a = next((r for r in static if r["mesh"] == mesh and r["case"] == "stage3_spar_cap"), None)
        b = next((r for r in static if r["mesh"] == mesh and r["case"] == "stage4_full"), None)
        if a and b:
            print(f"  {mesh}: tip {a['tip_uy_mm']:.6f} -> {b['tip_uy_mm']:.6f} mm "
                  f"({100*(b['tip_uy_mm']/a['tip_uy_mm']-1):+.3f}%) | design stress {a['skin_max_dTE_ge3']:.3f} -> "
                  f"{b['skin_max_dTE_ge3']:.3f} MPa ({100*(b['skin_max_dTE_ge3']/a['skin_max_dTE_ge3']-1):+.2f}%) | "
                  f"p99 {a['skin_p99']:.3f} -> {b['skin_p99']:.3f}")
        ba = next((r for r in buck if r["mesh"] == mesh and r["case"] == "stage3_spar_cap"), None)
        bb = next((r for r in buck if r["mesh"] == mesh and r["case"] == "stage4_full"), None)
        ta = next((r for r in twist if r["mesh"] == mesh and r["case"] == "stage3_spar_cap"), None)
        tb = next((r for r in twist if r["mesh"] == mesh and r["case"] == "stage4_full"), None)
        if ba and bb:
            print(f"        buckling {ba['critical']:.4f} -> {bb['critical']:.4f} "
                  f"({100*(bb['critical']/ba['critical']-1):+.2f}%)")
        if ta and tb:
            print(f"        tip twist {ta['theta_z1000_rad']:.3e} -> {tb['theta_z1000_rad']:.3e} rad")
    print("\nLoad check (re-read from written .inp):")
    for r in loadc:
        print(f"  {r['label']}: sumFy {r['sum_Fy']:.6f} (Fx {r['sum_Fx']:+.4f}, Fz {r['sum_Fz']:+.4f}) | "
              f"root moment err {r['root_moment_err_pct']:+.4f}% | sum(Fy*x) err {r['sum_FyX_err_pct']:+.3f}% | "
              f"xi_cp {r['xi_cp_min']:.4f}..{r['xi_cp_max']:.4f}")


if __name__ == "__main__":
    main()
