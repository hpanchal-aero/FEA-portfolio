#!/usr/bin/env python3
"""
Stage 6C: skin-thickness sensitivity of four rib layouts (L100, L126, L006, L151).
Reads results/stage6/records/<id>.json (1.0 mm) and <id>_sk0p8.json, <id>_sk0p6.json (+ any extra skin runs, which
are listed as verification points and NOT used in the fit).
Fit: ln(buckling) and ln(tip) quadratic in ln(t) through the three base thicknesses; mass linear in t.
Outputs: results/stage6c_skin_sensitivity.csv, figures/stage6c_skin_sensitivity.png
Run: /usr/bin/python3 scripts/analyse_skin_sensitivity.py
"""
import csv
import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

LAYOUTS = {"L100": 2, "L126": 5, "L006": 6, "L151": 10}
BASE_T = (1.0, 0.8, 0.6)
REQ = 1.5
RECS = "results/stage6/records"

data = {}
for L in LAYOUTS:
    pts = []
    for p in sorted(glob.glob("%s/%s.json" % (RECS, L)) + glob.glob("%s/%s_sk*.json" % (RECS, L))):
        d = json.load(open(p))
        if not d.get("complete"):
            continue
        m = d["metrics"]
        pts.append(dict(t=float(d.get("skin_t", 1.0)), buck=m["buckling_crit"], tip=m["tip_mm"], mass=m["mass_kg"],
                        p99=m.get("skin_p99", float("nan")), solid=m.get("solid_only_max", float("nan")),
                        b2b1=m["buckling_factors"][1] / m["buckling_factors"][0]))
    pts.sort(key=lambda x: -x["t"])
    data[L] = pts

print("%-5s %-5s %-6s %-8s %-8s %-8s %-8s %-8s %-6s" % ("id", "ribs", "t", "buck", "tip", "mass", "p99", "solidmx", "b2/b1"))
for L, pts in data.items():
    for x in pts:
        print("%-5s %-5d %-6.3g %-8.4f %-8.4f %-8.4f %-8.2f %-8.2f %-6.3f" % (
            L, LAYOUTS[L], x["t"], x["buck"], x["tip"], x["mass"], x["p99"], x["solid"], x["b2b1"]))

res = []
fits = {}
tg = np.linspace(0.5, 1.0, 5001)
print("\nlocal scaling exponents (buckling ~ t^k, tip ~ t^-j) between base points, and estimates at buckling = %.1f:" % REQ)
for L, pts in data.items():
    base = [x for x in pts if any(abs(x["t"] - b) < 1e-9 for b in BASE_T)]
    if len(base) != 3:
        raise SystemExit("ABORT: %s has %d base-thickness records (need 3)" % (L, len(base)))
    lt = np.log([x["t"] for x in base])
    pb = np.polyfit(lt, np.log([x["buck"] for x in base]), 2)
    pt = np.polyfit(lt, np.log([x["tip"] for x in base]), 2)
    pm = np.polyfit([x["t"] for x in base], [x["mass"] for x in base], 1)
    fits[L] = (pb, pt, pm)
    k1 = np.log(base[0]["buck"] / base[1]["buck"]) / np.log(base[0]["t"] / base[1]["t"])
    k2 = np.log(base[1]["buck"] / base[2]["buck"]) / np.log(base[1]["t"] / base[2]["t"])
    j1 = np.log(base[1]["tip"] / base[0]["tip"]) / np.log(base[0]["t"] / base[1]["t"])
    j2 = np.log(base[2]["tip"] / base[1]["tip"]) / np.log(base[1]["t"] / base[2]["t"])
    B = np.exp(np.polyval(pb, np.log(tg)))
    idx = np.where(B >= REQ)[0]
    tstar = float(tg[idx[0]]) if len(idx) else float("nan")
    # refine by linear interpolation in the grid
    mass = float(np.polyval(pm, tstar))
    tip = float(np.exp(np.polyval(pt, np.log(tstar))))
    extrap = tstar < 0.6
    print("%-5s ribs %-2d  k(1.0-0.8)=%.2f k(0.8-0.6)=%.2f | j(1.0-0.8)=%.2f j(0.8-0.6)=%.2f | t*=%.3f mm%s  mass*=%.4f kg  tip*=%.3f mm"
          % (L, LAYOUTS[L], k1, k2, j1, j2, tstar, " (EXTRAPOLATED below 0.6)" if extrap else "", mass, tip))
    res.append(dict(id=L, ribs=LAYOUTS[L], k_10_08=k1, k_08_06=k2, t_star=tstar, mass_star=mass, tip_star=tip, extrapolated=extrap))

os.makedirs("figures", exist_ok=True)
with open("results/stage6c_skin_sensitivity.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(res[0])); w.writeheader(); w.writerows(res)

colors = {"L100": "#4c72b0", "L126": "#dd8452", "L006": "#55a868", "L151": "#c44e52"}
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.8), constrained_layout=True)
for L, pts in data.items():
    pb, pt, pm = fits[L]
    B = np.exp(np.polyval(pb, np.log(tg)))
    ax[0].plot(tg, B, color=colors[L], label="%s (%d ribs)" % (L, LAYOUTS[L]))
    ax[0].plot([x["t"] for x in pts], [x["buck"] for x in pts], "o", color=colors[L])
    ax[1].plot(B, np.polyval(pm, tg), color=colors[L], label="%s (%d ribs)" % (L, LAYOUTS[L]))
    ax[1].plot([x["buck"] for x in pts], [x["mass"] for x in pts], "o", color=colors[L])
ax[0].axhline(REQ, color="k", ls="--", lw=1); ax[0].text(0.505, REQ * 1.03, "requirement 1.5", fontsize=8)
ax[0].set_xscale("log"); ax[0].set_yscale("log")
ax[0].set_xlabel("skin thickness [mm]"); ax[0].set_ylabel("critical buckling factor"); ax[0].grid(alpha=0.3, which="both")
ax[0].legend(fontsize=8)
ax[0].set_title("Buckling vs skin thickness (points = solver, lines = fit)", fontsize=10)
ax[1].axvline(REQ, color="k", ls="--", lw=1)
ax[1].set_xlim(1.0, 5.0); ax[1].set_ylim(1.0, 1.8)
ax[1].set_xlabel("critical buckling factor"); ax[1].set_ylabel("structural mass [kg]"); ax[1].grid(alpha=0.3)
ax[1].legend(fontsize=8)
ax[1].set_title("Mass vs buckling: lowest curve at 1.5 is lightest", fontsize=10)
fig.savefig("figures/stage6c_skin_sensitivity.png", dpi=170)
print("\nwrote figures/stage6c_skin_sensitivity.png, results/stage6c_skin_sensitivity.csv")
