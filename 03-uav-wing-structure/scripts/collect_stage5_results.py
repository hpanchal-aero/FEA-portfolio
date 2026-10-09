#!/usr/bin/env python3
"""
Collect the Stage 5 results: rib spacing 200 mm vs 100 mm, Stage 4 skin-pressure load, t = 1.0 mm, solid 4 / skin 8 mm
mesh, into results/stage5_comparison.csv and results/stage5_summary.json.

Definitions as Stage 3/4 (collect_stage3_results.stress_stats; design stress = max VM over expanded skin nodes with
dTE >= 3 mm and z > 10 mm, which is MESH-SENSITIVE, see the sensitivity columns chosen after inspecting the data).
Yield 290 MPa is ASSUMED (uncited); aircraft mass 20 kg is assumed. The load is an analytic placeholder.
The 50 mm and 20 mm entries are NOT RUN (hardware limit); their counts are copied from earlier console output.

Run (system python, from the project directory):  /usr/bin/python3 scripts/collect_stage5_results.py
"""
import contextlib
import csv
import io
import json
import os
import re
import subprocess
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from collect_stage3_results import stress_stats, YIELD, REQ_BUCKLE
from collect_stage4_results import buckle_factors, design_extras
from extract_stage4_twist import analyse

G = 9.81
CASES = [
    dict(label="200 mm", spacing=200, job="simulation/stage4_s4_k8_t1_n3p8"),
    dict(label="100 mm", spacing=100, job="simulation/stage4_s4_k8_rs100_t1_n3p8"),
]
NOT_RUN = {
    "50 mm": dict(status="built and meshed, NOT solved (hardware limit)", ribs=21, volume_mm3=223469.4464,
                  tets=112518, skin_shells=65800, elements=178318, nodes=334992, equation_estimate=1401108,
                  equation_ceiling=1180000, equation_pct_of_ceiling=118.7,
                  source="build_stage3_skin_inp.py console output on raw_s4_k8_rs50_closeedgefacnone.npz (needs NETGEN_KW=closeedgefac=None)"),
    "20 mm": dict(status="geometry built, NOT meshed (hardware limit)", ribs=51, volume_mm3=413368.7873,
                  tets_extrapolated_not_measured=246000, source="linear extrapolation of tet counts for 6, 11, 21 ribs; "
                  "the imprint also contains four sliver pieces (two of 0.0037 mm^2)"),
}


def jload(job):
    return json.load(open(os.path.join("results", os.path.basename(job) + "_results.json")))


def inp_counts(inp):
    nodes = c3 = s6 = 0
    skin_nodes = set()
    sec = ""
    with open(inp) as f:
        for line in f:
            if line.startswith("*") and not line.startswith("**"):
                sec = line.strip().upper().replace(" ", "")
                continue
            if line.startswith("**") or not line.strip():
                continue
            if sec == "*NODE":
                nodes += 1
            elif sec.startswith("*ELEMENT,TYPE=C3D10"):
                c3 += 1
            elif sec.startswith("*ELEMENT,TYPE=S6"):
                s6 += 1
                skin_nodes.update(int(v) for v in line.split(",")[1:7])
    return dict(nodes=nodes, c3d10=c3, s6=s6, skin_nodes=len(skin_nodes),
                equation_estimate=3 * nodes + 3 * len(skin_nodes))


def beam_numbers(spacing):
    env = dict(os.environ)
    env["RIB_SPACING"] = str(spacing)
    out = subprocess.run([sys.executable, "scripts/beam_with_rib_inertia.py"], env=env,
                         capture_output=True, text=True, check=True).stdout
    res = {}
    pat = re.compile(r"^(lift only|gravity only|full)\s*: no rib inertia\s+(-?[\d.]+) mm \| "
                     r"with rib-section inertia\s+(-?[\d.]+) mm")
    for line in out.splitlines():
        m = pat.match(line)
        if m:
            res[m.group(1)] = (float(m.group(2)), float(m.group(3)))
    if len(res) != 3:
        raise SystemExit(f"could not parse beam_with_rib_inertia.py output:\n{out}")
    return res


def main():
    rows = []
    for C in CASES:
        job = C["job"]
        for ext in (".inp", ".frd"):
            if not os.path.isfile(job + ext):
                raise SystemExit(f"missing {job + ext}")
        j, jl, jg = jload(job), jload(job + "_lift"), jload(job + "_grav")
        st, (vm, eid, ex, so, coords) = stress_stats(job)
        extras = design_extras(vm, eid, ex)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            th = analyse(job)
        facs = buckle_factors(job + "_buckle")
        beam = beam_numbers(C["spacing"])
        cnt = inp_counts(job + ".inp")
        n_ribs = int(round(1000 / C["spacing"])) + 1
        m_solid = j["W_solid_N"] / (G * j["n"])
        m_skin = j["W_skin_N"] / (G * j["n"])
        row = dict(label=C["label"], spacing_mm=C["spacing"], n_ribs=n_ribs,
                   mass_solid_kg=m_solid, mass_skin_kg=m_skin, mass_total_kg=m_solid + m_skin,
                   tip_full_mm=j["uy_spar_axis_mm"], tip_lift_only_mm=jl["uy_spar_axis_mm"],
                   tip_gravity_only_mm=jg["uy_spar_axis_mm"],
                   beam_full_no_rib_inertia_mm=beam["full"][0], beam_full_with_rib_inertia_mm=beam["full"][1],
                   beam_lift_with_rib_inertia_mm=beam["lift only"][1], beam_grav_with_rib_inertia_mm=beam["gravity only"][1],
                   fea_vs_beam_with_rib_inertia_full_pct=100 * (j["uy_spar_axis_mm"] / beam["full"][1] - 1),
                   fea_vs_beam_with_rib_inertia_lift_pct=100 * (jl["uy_spar_axis_mm"] / beam["lift only"][1] - 1),
                   tip_twist_rad=float(th[-1]), tip_twist_deg=float(np.degrees(th[-1])),
                   skin_p99_MPa=st["skin_p99"], skin_max_TE_vertex_MPa=st["skin_max"],
                   solid_only_max_MPa=st["solid_only_max"],
                   design_stress_dTE_ge3_MPa=st["skin_max_dTE_ge3"], design_stress_dTE_ge15_MPa=st["skin_max_dTE_ge15"],
                   design_max_excl_z_ge_990_dTE_lt_6_posthoc_MPa=extras["design_max_excl_z_ge_990_dTE_lt_6_posthoc"],
                   FoS_official_assumed_yield=YIELD / st["skin_max_dTE_ge3"],
                   FoS_dTE_lt_6_variant_assumed_yield=YIELD / extras["design_max_excl_z_ge_990_dTE_lt_6_posthoc"],
                   buckling_critical=min(facs), buckling_requirement=REQ_BUCKLE,
                   buckling_margin_pct=100 * (min(facs) / REQ_BUCKLE - 1),
                   **{f"buckling_mode{k}": f for k, f in enumerate(facs, 1)}, **cnt)
        rows.append(row)
        print(f"  done {C['label']}")

    a, b = rows
    delta = {k: 100 * (b[k] / a[k] - 1) for k in (
        "mass_total_kg", "tip_full_mm", "tip_lift_only_mm", "skin_p99_MPa", "design_stress_dTE_ge3_MPa",
        "design_max_excl_z_ge_990_dTE_lt_6_posthoc_MPa", "buckling_critical", "tip_twist_rad")}
    delta["tip_gravity_only_mm_magnitude"] = 100 * (abs(b["tip_gravity_only_mm"]) / abs(a["tip_gravity_only_mm"]) - 1)
    os.makedirs("results", exist_ok=True)
    with open("results/stage5_comparison.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    notes = ("Stage 5 (reduced scope): 200 mm vs 100 mm rib spacing, t = 1.0 mm, solid 4 / skin 8 mm, Stage 4 analytic "
             "placeholder pressure (not aerodynamic), 20 kg mass and 290 MPa yield assumed. 50 mm and 20 mm not run "
             "(11 GB RAM). The design-stress metric (dTE >= 3 mm) is mesh-sensitive and not comparable across spacings.")
    json.dump(dict(notes=notes, cases=rows, delta_100_vs_200_pct=delta, not_run=NOT_RUN),
              open("results/stage5_summary.json", "w"), indent=1)
    print("  wrote results/stage5_comparison.csv and results/stage5_summary.json\n")
    for k in ("mass_total_kg", "tip_full_mm", "tip_lift_only_mm", "tip_gravity_only_mm", "beam_full_with_rib_inertia_mm",
              "fea_vs_beam_with_rib_inertia_full_pct", "tip_twist_deg", "skin_p99_MPa", "design_stress_dTE_ge3_MPa",
              "design_max_excl_z_ge_990_dTE_lt_6_posthoc_MPa", "buckling_critical", "buckling_mode2", "buckling_mode6",
              "nodes", "c3d10", "s6", "equation_estimate"):
        print(f"  {k:50s} {a[k]:14.6g} {b[k]:14.6g}")
    print("\n  100 mm vs 200 mm [%]:", {k: round(v, 3) for k, v in delta.items()})


if __name__ == "__main__":
    main()
