"""
Project 02, Configuration B -- cross-design comparison table and figures.

Reads:
  - results/config_b_<shape>_<level>pct_h<mesh>_configB_results.json
    (mesh level per shape chosen at the CONVERGED level established in
    this project's convergence study: 2.0mm for rect/circle, 1.8mm for
    cross/grid -- NOT simply "finest available")
  - mesh/config_b_study/build_assembly_<shape>_<level>pct.log
    (mass, grepped from the build script's stdout -- not stored elsewhere)
  - results/config_a_h2p0_results.json (Config A baseline, 2.0mm, the
    converged Config A level established earlier)

Config A's "net-section stress" is NOT computed by the same method as
Config B's (nodal-mean-with-corner-exclusion over a pocket throat) --
it is the theoretical/FEA-confirmed uniform applied stress (-0.1643 MPa),
since Config A is prismatic. This is flagged in the output, not silently
treated as equivalent precision.

Outputs:
  results/config_b_comparison.csv
  figures/config_b_mass_vs_level.png
  figures/config_b_stiffness_vs_level.png
  figures/config_b_netsection_vs_level.png
  figures/config_b_stiffness_per_mass.png
"""
import csv
import glob
import json
import os
import re

import matplotlib.pyplot as plt

SHAPES = ["rect", "circle", "cross", "grid"]
LEVELS = [10, 20, 30]
MESH_FOR_SHAPE = {"rect": "h2p0", "circle": "h2p0", "cross": "h1p8", "grid": "h1p8"}

CONFIG_A_VOLUME_MM3 = 95300.0
RHO = 2700e-9  # t/mm^3 (== kg/m^3 numerically in this unit system)


def get_mass_from_log(shape, level):
    path = f"mesh/config_b_study/build_assembly_{shape}_{level}pct.log"
    if not os.path.isfile(path):
        return None
    text = open(path).read()
    m = re.search(r"Mass at Al 6061-T6 density:\s*([0-9.]+)\s*g", text)
    return float(m.group(1)) if m else None


def get_config_b_result(shape, level):
    mesh = MESH_FOR_SHAPE[shape]
    path = f"results/config_b_{shape}_{level}pct_{mesh}_configB_results.json"
    if not os.path.isfile(path):
        print(f"WARNING: missing {path}")
        return None
    d = json.load(open(path))
    d["_mesh_level_used"] = mesh
    return d


def get_config_a_baseline():
    path = "results/config_a_h2p0_results.json"
    d = json.load(open(path))
    force = d["applied_N"]
    stiffness = force / abs(d["uz_top_mean_mm"])
    mass_g = CONFIG_A_VOLUME_MM3 * RHO * 1000.0  # recomputed from validated volume
    return {
        "mass_g": mass_g,
        "stiffness_N_per_mm": stiffness,
        "net_section_szz_MPa": -0.1643,  # theoretical/uniform, NOT the same method as Config B
        "net_section_method": "uniform_theoretical",
        "mesh_level_used": "h2p0",
    }


def main():
    rows = []
    for shape in SHAPES:
        for level in LEVELS:
            mass = get_mass_from_log(shape, level)
            res = get_config_b_result(shape, level)
            if mass is None or res is None:
                print(f"SKIP {shape} {level}%: missing data (mass={mass}, res={'ok' if res else None})")
                continue
            rows.append({
                "shape": shape, "level_pct": level,
                "mesh_level": res["_mesh_level_used"],
                "mass_g": mass,
                "stiffness_N_per_mm": res["stiffness_N_per_mm"],
                "net_section_szz_MPa": res["net_section_szz_full_MPa"],
                "net_section_method": "nodal_mean_corner_excl_proxy",
                "peak_vm_MPa": res["peak_vm_MPa"],
                "equilibrium_err_pct": res["equilibrium_err_pct"],
            })

    a = get_config_a_baseline()
    rows.append({
        "shape": "config_a_baseline", "level_pct": 0,
        "mesh_level": a["mesh_level_used"],
        "mass_g": a["mass_g"],
        "stiffness_N_per_mm": a["stiffness_N_per_mm"],
        "net_section_szz_MPa": a["net_section_szz_MPa"],
        "net_section_method": a["net_section_method"],
        "peak_vm_MPa": None,
        "equilibrium_err_pct": None,
    })

    if len(rows) < 13:
        print(f"WARNING: expected 13 rows (12 Config B + 1 Config A baseline), got {len(rows)}")

    csv_path = "results/config_b_comparison.csv"
    fieldnames = ["shape", "level_pct", "mesh_level", "mass_g", "stiffness_N_per_mm",
                  "net_section_szz_MPa", "net_section_method", "peak_vm_MPa",
                  "equilibrium_err_pct"]
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"Wrote {csv_path} ({len(rows)} rows)")

    b_rows = [r for r in rows if r["shape"] != "config_a_baseline"]
    colors = {"rect": "tab:blue", "circle": "tab:orange", "cross": "tab:red", "grid": "tab:green"}

    def by_shape(rows_, key):
        out = {}
        for r in rows_:
            out.setdefault(r["shape"], []).append((r["level_pct"], r[key]))
        for s in out:
            out[s].sort()
        return out

    os.makedirs("figures", exist_ok=True)

    # 1. Mass vs removal level
    plt.figure(figsize=(6, 4.5))
    md = by_shape(b_rows, "mass_g")
    for s, pts in md.items():
        xs, ys = zip(*pts)
        plt.plot(xs, ys, "o-", label=s, color=colors[s])
    plt.axhline(a["mass_g"], color="k", ls="--", lw=1, label="Config A baseline")
    plt.xlabel("Panel-area removal level (%)")
    plt.ylabel("Structure mass (g)")
    plt.title("Configuration B: mass vs. removal level\n(identical across shapes at a given level, by design)")
    plt.legend(fontsize=8)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("figures/config_b_mass_vs_level.png", dpi=150)
    plt.close()

    # 2. Stiffness vs removal level
    plt.figure(figsize=(6, 4.5))
    sd = by_shape(b_rows, "stiffness_N_per_mm")
    for s, pts in sd.items():
        xs, ys = zip(*pts)
        plt.plot(xs, ys, "o-", label=s, color=colors[s])
    plt.axhline(a["stiffness_N_per_mm"], color="k", ls="--", lw=1, label="Config A baseline")
    plt.xlabel("Panel-area removal level (%)")
    plt.ylabel("Axial stiffness F/|Uz| (N/mm)")
    plt.title("Configuration B: axial stiffness vs. removal level\n(mesh level per shape: converged level, see CSV)")
    plt.legend(fontsize=8)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("figures/config_b_stiffness_vs_level.png", dpi=150)
    plt.close()

    # 3. Net-section stress vs removal level
    plt.figure(figsize=(6, 4.5))
    nd = by_shape(b_rows, "net_section_szz_MPa")
    for s, pts in nd.items():
        xs, ys = zip(*pts)
        plt.plot(xs, [abs(y) for y in ys], "o-", label=s, color=colors[s])
    plt.axhline(abs(a["net_section_szz_MPa"]), color="k", ls="--", lw=1,
                label="Config A (uniform theoretical)")
    plt.xlabel("Panel-area removal level (%)")
    plt.ylabel("|Net-section mean SZZ| (MPa)")
    plt.title("Configuration B: net-section stress vs. removal level\n"
              "(nodal-mean corner-exclusion proxy; NOT an integrated force/area value;\n"
              "cross's small section-node count makes it least reliable)")
    plt.legend(fontsize=8)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("figures/config_b_netsection_vs_level.png", dpi=150)
    plt.close()

    # 4. Stiffness-per-mass (structural efficiency)
    plt.figure(figsize=(6, 4.5))
    eff = {r["shape"]: [] for r in b_rows}
    for r in b_rows:
        eff[r["shape"]].append((r["level_pct"], r["stiffness_N_per_mm"] / r["mass_g"]))
    for s in eff:
        eff[s].sort()
        xs, ys = zip(*eff[s])
        plt.plot(xs, ys, "o-", label=s, color=colors[s])
    plt.axhline(a["stiffness_N_per_mm"] / a["mass_g"], color="k", ls="--", lw=1,
                label="Config A baseline")
    plt.xlabel("Panel-area removal level (%)")
    plt.ylabel("Stiffness / mass (N/mm per g)")
    plt.title("Configuration B: structural efficiency vs. removal level")
    plt.legend(fontsize=8)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("figures/config_b_stiffness_per_mass.png", dpi=150)
    plt.close()

    print("\nWrote figures:")
    for f in ["config_b_mass_vs_level.png", "config_b_stiffness_vs_level.png",
              "config_b_netsection_vs_level.png", "config_b_stiffness_per_mass.png"]:
        print(f"  figures/{f}")

    print("\n--- Summary table ---")
    print(f"{'shape':<12}{'level':>6}{'mass_g':>9}{'stiff_N/mm':>13}{'net_MPa':>10}{'eff':>9}")
    for r in sorted(rows, key=lambda r: (r['shape'], r['level_pct'])):
        eff_v = r['stiffness_N_per_mm'] / r['mass_g']
        print(f"{r['shape']:<12}{r['level_pct']:>6}{r['mass_g']:>9.2f}"
              f"{r['stiffness_N_per_mm']:>13.1f}{r['net_section_szz_MPa']:>10.4f}{eff_v:>9.1f}")


if __name__ == "__main__":
    main()
