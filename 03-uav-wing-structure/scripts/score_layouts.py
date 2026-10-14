#!/usr/bin/env python3
"""
Stage 6: score all completed layout runs (69 layouts at 1.5 mm rib + 21 selected layouts at 1.2/1.0/0.8 mm).
Feasibility filter: critical buckling factor >= 1.5 (limit load).
Headline score : (buckling / mass) relative to L000 at 1.5 mm (buckling margin per kg; no weights).
Composite score: min-max normalised over the pool (buckling up, mass down, tip down),
                 weights mass 0.45, buckling 0.40, tip 0.15 (judgement for a small UAV wing, not derived).
Sensitivity    : weights swept (mass 0.30-0.60, tip 0.05-0.25, buckling = remainder, >= 0.15);
                 counts how often each run ranks first / in the top 5.
Differences below the discretization scatter (buckling about 0.3%, tip about 0.15%) are ties.
Writes results/stage6_scores.csv. Usage: /usr/bin/python3 scripts/score_layouts.py
"""
import csv
import itertools

FILES = ["results/stage6_results.csv", "results/stage6_results_rt1p2.csv",
         "results/stage6_results_rt1.csv", "results/stage6_results_rt0p8.csv"]
W = dict(mass=0.45, buck=0.40, tip=0.15)
REQ = 1.5

runs = []
for p in FILES:
    for r in csv.DictReader(open(p)):
        if r["status"] != "complete":
            continue
        runs.append(dict(id=r["id"], rt=float(r["rib_t"]), stations=r["stations"], family=r["family"],
                         buck=float(r["buckling_crit"]), tip=float(r["tip_mm"]), mass=float(r["mass_kg"])))
n0 = len(runs)
runs = [r for r in runs if r["buck"] >= REQ]
base = next(r for r in runs if r["id"] == "L000" and r["rt"] == 1.5)
print("runs: %d complete, %d pass buckling >= %.1f ; baseline L000 1.5 mm: buck %.4f mass %.4f tip %.4f"
      % (n0, len(runs), REQ, base["buck"], base["mass"], base["tip"]))

lo = {k: min(r[k] for r in runs) for k in ("buck", "tip", "mass")}
hi = {k: max(r[k] for r in runs) for k in ("buck", "tip", "mass")}
for r in runs:
    r["name"] = "%s@%g" % (r["id"], r["rt"])
    r["headline"] = (r["buck"] / r["mass"]) / (base["buck"] / base["mass"])
    r["nb"] = (r["buck"] - lo["buck"]) / (hi["buck"] - lo["buck"])
    r["nm"] = (hi["mass"] - r["mass"]) / (hi["mass"] - lo["mass"])
    r["nt"] = (hi["tip"] - r["tip"]) / (hi["tip"] - lo["tip"])


def comp(r, w):
    return w["buck"] * r["nb"] + w["mass"] * r["nm"] + w["tip"] * r["nt"]


for r in runs:
    r["composite"] = comp(r, W)

# sensitivity
first, top5, ncomb = {}, {}, 0
for wm in [x / 100 for x in range(30, 61, 5)]:
    for wt in [x / 100 for x in range(5, 26, 5)]:
        wb = 1.0 - wm - wt
        if wb < 0.15 - 1e-9:
            continue
        ncomb += 1
        rk = sorted(runs, key=lambda r: -comp(r, dict(mass=wm, tip=wt, buck=wb)))
        first[rk[0]["name"]] = first.get(rk[0]["name"], 0) + 1
        for r in rk[:5]:
            top5[r["name"]] = top5.get(r["name"], 0) + 1
for r in runs:
    r["first"] = first.get(r["name"], 0)
    r["top5"] = top5.get(r["name"], 0)

by_h = sorted(runs, key=lambda r: -r["headline"])
by_c = sorted(runs, key=lambda r: -r["composite"])
for i, r in enumerate(by_h, 1):
    r["rank_headline"] = i
for i, r in enumerate(by_c, 1):
    r["rank_composite"] = i

print("\nHEADLINE (buckling per kg relative to L000@1.5), top 12:")
for r in by_h[:12]:
    print("  %2d %-10s %-34s score %.3f  buck %.3f mass %.4f tip %.4f" % (
        r["rank_headline"], r["name"], r["stations"], r["headline"], r["buck"], r["mass"], r["tip"]))
print("\nCOMPOSITE (mass %.2f, buck %.2f, tip %.2f), top 12:" % (W["mass"], W["buck"], W["tip"]))
for r in by_c[:12]:
    print("  %2d %-10s %-34s score %.3f  buck %.3f mass %.4f tip %.4f" % (
        r["rank_composite"], r["name"], r["stations"], r["composite"], r["buck"], r["mass"], r["tip"]))
print("\nSENSITIVITY over %d weight combinations: runs ranked first / in top 5" % ncomb)
for r in sorted(runs, key=lambda r: (-r["first"], -r["top5"]))[:8]:
    print("  %-10s first %2d/%d  top5 %2d/%d" % (r["name"], r["first"], ncomb, r["top5"], ncomb))
print("\nL000@1.5 baseline: headline rank %d, composite rank %d of %d" % (
    base["rank_headline"], base["rank_composite"], len(runs)))

with open("results/stage6_scores.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["run", "id", "rib_t", "stations", "family", "buckling", "tip_mm", "mass_kg", "headline",
                "rank_headline", "composite", "rank_composite", "first_of_%d" % ncomb, "top5_of_%d" % ncomb])
    for r in by_c:
        w.writerow([r["name"], r["id"], r["rt"], r["stations"], r["family"], r["buck"], r["tip"], r["mass"],
                    "%.4f" % r["headline"], r["rank_headline"], "%.4f" % r["composite"], r["rank_composite"],
                    r["first"], r["top5"]])
print("\nwrote results/stage6_scores.csv")
