#!/usr/bin/env python3
"""Read-only facts sheet for the Stage 6 write-up: counts and extremes per results CSV. Run: /usr/bin/python3 scripts/stage6_facts.py"""
import csv
import glob
import os


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


for path in ["results/stage6_results.csv"] + sorted(glob.glob("results/stage6_results_rt*.csv")):
    rows = list(csv.DictReader(open(path)))
    comp = [r for r in rows if r.get("status") == "complete" and f(r.get("buckling_crit")) is not None]
    stat = {}
    for r in rows:
        k = r.get("status", "?").split(":")[0]
        stat[k] = stat.get(k, 0) + 1
    print("\n=== %s : %d rows, status %s" % (path, len(rows), stat))
    if not comp:
        continue
    show = lambda tag, r: print("  %-14s %-5s buck %.4f tip %.4f mass %.4f rib_t %s stations %s" % (
        tag, r["id"], f(r["buckling_crit"]), f(r["tip_mm"]), f(r["mass_kg"]), r.get("rib_t"), r["stations"]))
    show("best buckling", max(comp, key=lambda r: f(r["buckling_crit"])))
    show("lowest tip", min(comp, key=lambda r: f(r["tip_mm"])))
    show("lowest mass", min(comp, key=lambda r: f(r["mass_kg"])))
    show("best buck/kg", max(comp, key=lambda r: f(r["buckling_crit"]) / f(r["mass_kg"])))
    for r in comp:
        if r["id"] == "L000":
            show("L000 baseline", r)

for p in ["results/stage6_scores.csv", "results/stage6_tier2_selection.csv"]:
    if os.path.exists(p):
        rows = list(csv.reader(open(p)))
        print("\n=== %s : %d rows ; header + first 6:" % (p, len(rows) - 1))
        for r in rows[:7]:
            print("  ", ",".join(r))
