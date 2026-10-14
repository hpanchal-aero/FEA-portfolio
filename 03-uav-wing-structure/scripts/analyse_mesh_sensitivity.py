#!/usr/bin/env python3
"""
Stage 6 mesh sensitivity: skin-mesh size (solid mesh fixed at 4 mm) against the 8 mm baseline, per layout and skin gauge.
Reads results/stage6/records/<id>[_skXpY]_s4_k<K>.json (conv runs) and the matching base record <id>[_skXpY].json.
Writes results/stage6d_mesh_sensitivity.csv and figures/stage6d_mesh_sensitivity.png
Run: /usr/bin/python3 scripts/analyse_mesh_sensitivity.py
"""
import csv
import glob
import json
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = "results/stage6/records"
rows = {}     # (layout, skin) -> {k: metrics}
for p in sorted(glob.glob(R + "/L*_s4_k*.json")):
    key = p.split("/")[-1][:-5]
    m = re.match(r"(L\d+)(?:_sk([0-9p]+))?_s4_k([0-9p]+)$", key)
    if not m:
        continue
    lay, sk, k = m.group(1), m.group(2), float(m.group(3).replace("p", "."))
    skin = float(sk.replace("p", ".")) if sk else 1.0
    d = json.load(open(p))
    if not d["complete"]:
        raise SystemExit("ABORT: %s is not complete (failed step %s)" % (key, d.get("failed_step")))
    base_p = "%s/%s%s.json" % (R, lay, ("_sk" + sk) if sk else "")
    base = json.load(open(base_p))
    if not base["complete"]:
        raise SystemExit("ABORT: base record %s is not complete" % base_p)
    case = rows.setdefault((lay, skin), {})
    case[8.0] = dict(m=base["metrics"], eq=base["steps"]["inp"].get("eq_pct"), s6=base["steps"]["inp"].get("s6"))
    case[k] = dict(m=d["metrics"], eq=d["steps"]["inp"].get("eq_pct"), s6=d["steps"]["inp"].get("s6"))

out = []
print("%-5s %-5s %-5s %-9s %-9s %-8s %-9s %-9s %-7s %-6s" % ("lay", "skin", "k", "buck", "d_buck%", "tip", "d_tip%", "d_p99%", "eq%", "S6"))
for (lay, skin), case in sorted(rows.items()):
    b = case[8.0]["m"]
    for k in sorted(case, reverse=True):
        m = case[k]["m"]
        dp = 100 * (m["skin_p99"] / b["skin_p99"] - 1) if b.get("skin_p99") and m.get("skin_p99") else None
        rec = dict(layout=lay, skin_t=skin, skin_mesh=k, buckling=m["buckling_crit"],
                   d_buck_pct=100 * (m["buckling_crit"] / b["buckling_crit"] - 1), tip_mm=m["tip_mm"],
                   d_tip_pct=100 * (m["tip_mm"] / b["tip_mm"] - 1), mass_kg=m["mass_kg"],
                   d_p99_pct=dp, eq_pct=case[k]["eq"], s6=case[k]["s6"])
        out.append(rec)
        print("%-5s %-5g %-5g %-9.4f %+9.3f %-8.4f %+9.4f %-9s %-7s %-6s" % (
            lay, skin, k, rec["buckling"], rec["d_buck_pct"], rec["tip_mm"], rec["d_tip_pct"],
            "n/a" if dp is None else "%+.2f" % dp, rec["eq_pct"], rec["s6"]))
with open("results/stage6d_mesh_sensitivity.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)

print("\nmax |tip change| over all runs: %.4f %%" % max(abs(r["d_tip_pct"]) for r in out))
print("max |mass change| over all runs: %.4f %%" % max(abs(100 * (r["mass_kg"] / [o for o in out if o["layout"] == r["layout"] and o["skin_t"] == r["skin_t"] and o["skin_mesh"] == 8.0][0]["mass_kg"] - 1)) for r in out))
ps = [r["d_p99_pct"] for r in out if r["d_p99_pct"] is not None]
print("p99 change range: %+.2f %% to %+.2f %%" % (min(ps), max(ps)))

fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.6), constrained_layout=True)
col = {"L100": "#4c72b0", "L126": "#dd8452", "L006": "#55a868"}
for (lay, skin), case in sorted(rows.items()):
    ks = sorted(case)
    ax[0].plot(ks, [100 * (case[k]["m"]["buckling_crit"] / case[8.0]["m"]["buckling_crit"] - 1) for k in ks], "o",
               color=col[lay], ls="-" if skin == 1.0 else "--", label="%s, skin %g mm" % (lay, skin))
ax[0].axhline(0, color="k", lw=0.8)
ax[0].set_xlabel("skin mesh size [mm]  (finer to the left)"); ax[0].set_ylabel("buckling change vs 8 mm mesh [%]")
ax[0].invert_xaxis(); ax[0].grid(alpha=0.3); ax[0].legend(fontsize=8)
ax[0].set_title("Buckling vs skin mesh size (solid mesh 4 mm)", fontsize=10)
c = rows[("L100", 1.0)]
ks = sorted(c)
ax[1].plot(ks, [c[k]["m"]["buckling_crit"] for k in ks], "o-", color=col["L100"])
for k in ks:
    ax[1].annotate("%.4f" % c[k]["m"]["buckling_crit"], (k, c[k]["m"]["buckling_crit"]), textcoords="offset points", xytext=(0, 7), ha="center", fontsize=8)
ax[1].set_xlabel("skin mesh size [mm]  (finer to the left)"); ax[1].set_ylabel("critical buckling factor")
ax[1].invert_xaxis(); ax[1].grid(alpha=0.3)
ax[1].set_title("L100, 1.0 mm skin: change flattens from 6 to 5 mm", fontsize=10)
fig.savefig("figures/stage6d_mesh_sensitivity.png", dpi=170)
print("\nwrote results/stage6d_mesh_sensitivity.csv, figures/stage6d_mesh_sensitivity.png")
