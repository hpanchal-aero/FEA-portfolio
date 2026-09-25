"""
Project 02, Configuration B random geometries -- full feasibility table
against the element/equation ceilings, for all 223 candidates at all 3
mesh sizes (1.8, 2.0, 3.0mm).

C3D10 node/equation counts here are ESTIMATES: linear-mesh tets * a
promotion ratio empirically measured on one actual matching Config B
geometry (cross_10pct at 3.0mm: 28,313 linear tets -> 58,000 C3D10
nodes after promotion, ratio 2.049). This is NOT an exact count --
exact counts only come from actually running the C3D10 promotion in
mesh_config_b.py. Flagged in the output header and every row.

Reads: results/config_b_random_mesh_feasibility.csv
Writes: results/config_b_random_feasibility_table.csv (full, all 223x3 rows)
Prints: per-size table to stdout (sorted by estimated equation % ascending)
"""
import csv

PROMOTION_RATIO = 58000 / 28313  # C3D10 nodes per linear tet, from cross_10pct@3.0mm ground truth
ELEM_CEILING = 230_000
EQ_CEILING = 1_180_000

IN_CSV = "results/config_b_random_mesh_feasibility.csv"
OUT_CSV = "results/config_b_random_feasibility_table.csv"


def main():
    rows = []
    with open(IN_CSV) as f:
        for r in csv.DictReader(f):
            if r["status"] != "ok":
                continue
            tets = int(r["tets"])
            linear_nodes = int(r["linear_nodes"])
            est_c3d10_nodes = round(tets * PROMOTION_RATIO)
            est_eq = 3 * est_c3d10_nodes
            rows.append({
                "tag": r["tag"], "maxh": r["maxh"],
                "linear_nodes": linear_nodes, "tets": tets,
                "est_c3d10_nodes": est_c3d10_nodes,
                "est_elem_pct_of_ceiling": round(100 * tets / ELEM_CEILING, 2),
                "est_eq": est_eq,
                "est_eq_pct_of_ceiling": round(100 * est_eq / EQ_CEILING, 2),
                "gen_time_s": r["gen_time_s"],
            })

    fieldnames = ["tag", "maxh", "linear_nodes", "tets", "est_c3d10_nodes",
                  "est_elem_pct_of_ceiling", "est_eq", "est_eq_pct_of_ceiling", "gen_time_s"]
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"Wrote full table: {OUT_CSV} ({len(rows)} rows)")
    print(f"Promotion ratio used (C3D10 nodes / linear tets): {PROMOTION_RATIO:.4f} "
          f"(estimate, grounded on cross_10pct@3.0mm actual data)\n")

    for size in ("1.8", "2.0", "3.0"):
        size_rows = sorted((r for r in rows if r["maxh"] == size),
                            key=lambda r: -r["est_eq_pct_of_ceiling"])
        print("=" * 100)
        print(f"MESH SIZE {size}mm -- {len(size_rows)} geometries, sorted by est. equation "
              f"% of ceiling (descending)")
        print("=" * 100)
        print(f"{'tag':<22}{'tets':>8}{'elem%':>8}{'est_nodes':>11}{'est_eq':>9}{'eq%':>8}")
        for r in size_rows:
            print(f"{r['tag']:<22}{r['tets']:>8}{r['est_elem_pct_of_ceiling']:>7.1f}%"
                  f"{r['est_c3d10_nodes']:>11}{r['est_eq']:>9}{r['est_eq_pct_of_ceiling']:>7.1f}%")
        elem_pcts = [r["est_elem_pct_of_ceiling"] for r in size_rows]
        eq_pcts = [r["est_eq_pct_of_ceiling"] for r in size_rows]
        print(f"\n  Summary: elem% min={min(elem_pcts):.1f} max={max(elem_pcts):.1f} "
              f"mean={sum(elem_pcts)/len(elem_pcts):.1f}  |  "
              f"eq% min={min(eq_pcts):.1f} max={max(eq_pcts):.1f} "
              f"mean={sum(eq_pcts)/len(eq_pcts):.1f}\n")


if __name__ == "__main__":
    main()
