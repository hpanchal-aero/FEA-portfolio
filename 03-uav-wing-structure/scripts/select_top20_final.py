#!/usr/bin/env python3
"""
Stage 6: final tier 2 selection. Front 1 (strict non-dominated sorting, no tolerance) in full, then the
remainder of the 20 from front 2 by NSGA-II crowding distance computed within front 2 on objectives
normalised to 0-1 over front 2 (ties: lower layout id). L000 is added as a reference if not selected.
Reads results/stage6_fronts.csv (from select_top20.py --first-pass with default tolerances).
Writes results/stage6_tier2_selection.csv.
"""
import csv

rows = list(csv.DictReader(open("results/stage6_fronts.csv")))
for r in rows:
    r["front"] = int(r["front"])
    for k in ("buckling", "tip_mm", "mass_kg"):
        r[k] = float(r[k])
f1 = [r for r in rows if r["front"] == 1]
f2 = [r for r in rows if r["front"] == 2]
need = 20 - len(f1)
dist = {r["id"]: 0.0 for r in f2}
for key in ("buckling", "tip_mm", "mass_kg"):
    s = sorted(f2, key=lambda r: (r[key], r["id"]))
    lo, hi = s[0][key], s[-1][key]
    dist[s[0]["id"]] = dist[s[-1]["id"]] = float("inf")
    for i in range(1, len(s) - 1):
        if dist[s[i]["id"]] != float("inf"):
            dist[s[i]["id"]] += (s[i + 1][key] - s[i - 1][key]) / (hi - lo)
rank = sorted(f2, key=lambda r: (-dist[r["id"]], r["id"]))
print("front 1: %d layouts; need %d from front 2 (%d layouts)" % (len(f1), need, len(f2)))
print("front 2 ranked by crowding distance:")
for i, r in enumerate(rank, 1):
    d = dist[r["id"]]
    print("  %2d %s %-34s buck %.4f tip %.4f mass %.4f  crowd %s %s" % (
        i, r["id"], r["stations"], r["buckling"], r["tip_mm"], r["mass_kg"],
        "inf" if d == float("inf") else "%.3f" % d, "<- selected" if i <= need else ""))
sel = f1 + rank[:need]
ids = {r["id"] for r in sel}
with open("results/stage6_tier2_selection.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["id", "stations", "source", "buckling", "tip_mm", "mass_kg"])
    for r in sel:
        w.writerow([r["id"], r["stations"], "front%d" % r["front"], r["buckling"], r["tip_mm"], r["mass_kg"]])
    if "L000" not in ids:
        b = next(r for r in rows if r["id"] == "L000")
        w.writerow([b["id"], b["stations"], "reference (front %d)" % b["front"], b["buckling"], b["tip_mm"], b["mass_kg"]])
        print("\nL000 not selected (front %d): added as reference" % b["front"])
print("\nselected %d + reference -> results/stage6_tier2_selection.csv" % len(sel))
