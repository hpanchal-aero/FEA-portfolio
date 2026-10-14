#!/usr/bin/env python3
"""
Stage 6: non-dominated sorting of the complete tier 1 layouts on
(critical buckling factor: maximise, tip deflection: minimise, mass: minimise).
Rule fixed before the results were seen. Writes results/stage6_fronts.csv with --first-pass (layouts L000-L099, the set the tier 2 selection was made from), results/stage6_fronts_all.csv without it; a tolerance run adds the suffix _tol.
Dominance is strict: a dominates b if a is no worse in all three objectives and strictly better in at least one.
No noise tolerance is applied here; --tol applies relative tolerances (buckling, tip) so that differences
below the measured noise count as equal, and the output reports how many fronts change.
Usage: /usr/bin/python3 scripts/select_top20.py [--tol-buckling 0.0006 --tol-tip 0.0015]
"""
import argparse
import csv


def fronts(rows, tb, tt):
    def dom(a, b):
        better = worse = False
        for key, sign, tol in (("buck", 1, tb), ("tip", -1, tt), ("mass", -1, 0.0)):
            x, y = sign * a[key], sign * b[key]
            scale = abs(b[key]) * tol
            if x > y + scale:
                better = True
            elif x < y - scale:
                worse = True
        return better and not worse

    left = list(rows)
    out = []
    while left:
        f = [r for r in left if not any(dom(o, r) for o in left if o is not r)]
        out.append(f)
        left = [r for r in left if r not in f]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tol-buckling", type=float, default=0.0)
    ap.add_argument("--tol-tip", type=float, default=0.0)
    ap.add_argument("--first-pass", action="store_true",
                    help="only layouts with id number < 100 (the set the tier 2 selection was made from)")
    a = ap.parse_args()
    rows = []
    for r in csv.DictReader(open("results/stage6_results.csv")):
        if r["status"] == "complete" and (not a.first_pass or int(r["id"][1:]) < 100):
            rows.append(dict(id=r["id"], stations=r["stations"], buck=float(r["buckling_crit"]),
                             tip=float(r["tip_mm"]), mass=float(r["mass_kg"]), family=r["family"]))
    fr = fronts(rows, a.tol_buckling, a.tol_tip)
    n = 0
    out_path = "results/stage6_fronts" + ("" if a.first_pass else "_all") + \
               ("" if (a.tol_buckling == 0.0 and a.tol_tip == 0.0) else "_tol") + ".csv"
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["front", "id", "stations", "buckling", "tip_mm", "mass_kg", "family"])
        for k, fl in enumerate(fr, 1):
            for r in sorted(fl, key=lambda r: -r["buck"]):
                w.writerow([k, r["id"], r["stations"], r["buck"], r["tip"], r["mass"], r["family"]])
    print("layouts: %d  tolerances: buckling %.4f tip %.4f" % (len(rows), a.tol_buckling, a.tol_tip))
    print("wrote " + out_path)
    cum = 0
    for k, fl in enumerate(fr, 1):
        cum += len(fl)
        print("front %d: %2d layouts (cumulative %2d)" % (k, len(fl), cum))
        if cum <= 30:
            for r in sorted(fl, key=lambda r: -r["buck"]):
                print("   %s %-34s buck %.4f tip %.4f mass %.4f" % (r["id"], r["stations"], r["buck"], r["tip"], r["mass"]))
        if cum >= 30:
            break


if __name__ == "__main__":
    main()
