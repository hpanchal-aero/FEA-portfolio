#!/usr/bin/python3
"""Project 05: aggregate the study cases in results/cases/ into summary tables.

Reads results/cases/NAME.txt and NAME_history.csv for every case; makes no solver run.
Writes results/study_summary.csv, results/study_summary.md and results/mesh_sensitivity.csv and prints the tables.
FoS_yield_bound = sigma_y / (1.267 x maximum von Mises over the pulse): bounded elastic value (Amendment 11), not an allowable.
Usage: /usr/bin/python3 scripts/aggregate_study.py
"""
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from material import MATERIAL

ROOT = os.path.dirname(HERE)
CASES = os.path.join(ROOT, "results", "cases")
L = W = 0.5
SY = MATERIAL["sigma_y"] / 1e6
RHO = MATERIAL["rho"]
TIMES = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0]
BCORD = {"rest": 0, "free": 1, "clamp": 2}


def parse_name(name):
    m = re.match(r"^h([0-9p]+)_(rest|free|clamp)_a(m?)([0-9p]+)_n(\d+)$", name)
    if not m:
        raise SystemExit("unexpected case name: " + name)
    a = float(m.group(4).replace("p", "."))
    if m.group(3) == "m":
        a = -a
    return float(m.group(1).replace("p", ".")), m.group(2), a, int(m.group(5))


def load(name):
    h, bc, a, n = parse_name(name)
    d = np.genfromtxt(os.path.join(CASES, name + "_history.csv"), delimiter=",", names=True)
    txt = open(os.path.join(CASES, name + ".txt")).read()
    t = d["time_s"]
    if len(t) < 60 or abs(t[-1] - 30.0) > 1e-6:
        raise SystemExit("%s: history has %d rows, last time %g" % (name, len(t), t[-1]))

    def at(col, time):
        k = int(np.argmin(np.abs(t - time)))
        if abs(t[k] - time) > 1e-6:
            raise SystemExit("%s: no output at t = %g s" % (name, time))
        return float(d[col][k])

    tp = float(re.search(r"peak temperature at 30 s ([0-9.]+) C", txt).group(1))
    cb = re.search(r"cutbacks at t = \[(.*?)\] s", txt).group(1)
    wt = np.array([at("W_centre_over_h", x) for x in TIMES])
    vm = float(np.max(d["max_vM_MPa"]))
    bound = float(np.max(d["bound_MPa"]))
    return {"name": name, "h": h, "bc": bc, "a": a, "n": n, "mass": RHO * L * W * h * 1e-3, "tpeak": tp,
            "ncut": len([x for x in cb.split(",") if x.strip()]), "passed": "RESULT: PASS" in txt,
            "Wt": wt, "W5": wt[0], "W30": wt[-1], "W30mm": wt[-1] * h, "vm": vm, "bound": bound, "fos": SY / bound,
            "Rx30": abs(at("Rx_face_x0_N", 30.0)), "Rxmax": float(np.max(np.abs(d["Rx_face_x0_N"])))}


def main():
    names = sorted(f[:-len("_history.csv")] for f in os.listdir(CASES) if f.endswith("_history.csv"))
    recs = [load(n) for n in names]
    print("cases found: %d (expected 14)" % len(recs))
    bad = [r["name"] for r in recs if not r["passed"]]
    print("cases without RESULT: PASS in their report: %s" % (bad if bad else "none"))
    recs.sort(key=lambda r: (r["h"], BCORD[r["bc"]], r["a"], r["n"]))
    idx = dict(((r["h"], r["bc"], r["a"], r["n"]), r) for r in recs)

    cols = ["case", "h_mm", "edge_condition", "A_over_h", "n_inplane", "mass_kg", "T_peak_30s_C", "cutbacks",
            "W_c_over_h_5s", "W_c_over_h_30s", "W_c_mm_30s", "max_vM_MPa", "vM_bound_MPa", "FoS_yield_bound",
            "Rx_x0_30s_N", "Rx_x0_max_N"]
    rows = [[r["name"], "%g" % r["h"], r["bc"], "%g" % r["a"], "%d" % r["n"], "%.4f" % r["mass"], "%.2f" % r["tpeak"],
             "%d" % r["ncut"], "%.4f" % r["W5"], "%.4f" % r["W30"], "%.3f" % r["W30mm"], "%.2f" % r["vm"],
             "%.2f" % r["bound"], "%.2f" % r["fos"], "%.1f" % r["Rx30"], "%.1f" % r["Rxmax"]] for r in recs]
    with open(os.path.join(ROOT, "results", "study_summary.csv"), "w") as f:
        f.write(",".join(cols) + "\n")
        f.writelines(",".join(r) + "\n" for r in rows)
    with open(os.path.join(ROOT, "results", "study_summary.md"), "w") as f:
        f.write("| " + " | ".join(cols) + " |\n|" + "---|" * len(cols) + "\n")
        f.writelines("| " + " | ".join(r) + " |\n" for r in rows)
        f.write("\nFoS_yield_bound = 880 MPa / (1.267 x maximum von Mises over the pulse); bounded elastic value, not an allowable.\n")
    print("\n=== all cases (h in mm; stresses in MPa; Rx in N) ===")
    print("%-22s %4s %-6s %6s %3s %6s %7s %3s %7s %7s %7s %7s %7s %6s %8s" % (
        "case", "h", "edge", "A/h", "n", "mass", "Tpeak", "cb", "W5", "W30", "vM", "bound", "FoS", "", "Rx30"))
    for r in recs:
        print("%-22s %4g %-6s %6g %3d %6.3f %7.2f %3d %7.4f %7.4f %7.2f %7.2f %7.2f %6s %8.1f" % (
            r["name"], r["h"], r["bc"], r["a"], r["n"], r["mass"], r["tpeak"], r["ncut"], r["W5"], r["W30"], r["vm"],
            r["bound"], r["fos"], "", r["Rx30"]))
    fos = [(r["fos"], r["name"]) for r in recs]
    print("lowest FoS_yield_bound: %.2f (%s); highest: %.2f (%s)" % (min(fos)[0], min(fos)[1], max(fos)[0], max(fos)[1]))

    print("\n=== edge-condition sensitivity relative to simply supported restrained (A/h = +0.1, 20 x 20) ===")
    print("h[mm]  edge   dW_c/h[%]  dvM[%]  Rx30[N]")
    for h in (1.0, 2.0, 3.0):
        ref = idx[(h, "rest", 0.1, 20)]
        for bc in ("free", "clamp"):
            r = idx[(h, bc, 0.1, 20)]
            print("%4g   %-6s %8.1f %8.1f %9.1f" % (h, bc, 100 * (r["W30"] / ref["W30"] - 1), 100 * (r["vm"] / ref["vm"] - 1), r["Rx30"]))

    print("\n=== imperfection sensitivity (2 mm, restrained, 20 x 20) ===")
    print("A/h     W5     W30     vM[MPa]  Rx30[N]  cutbacks")
    for a in sorted(k[2] for k in idx if k[0] == 2.0 and k[1] == "rest" and k[3] == 20):
        r = idx[(2.0, "rest", a, 20)]
        print("%5g %7.4f %7.4f %8.2f %8.1f %5d" % (a, r["W5"], r["W30"], r["vm"], r["Rx30"], r["ncut"]))

    print("\n=== mesh sensitivity, 20 x 20 -> 40 x 40 (restrained, A/h = +0.1) ===")
    out = ["h_mm,max_rel_dW_5_to_30s,rel_dvM,rel_dRx30,dTpeak_K,dTpeak_over_rise"]
    for h in (1.0, 2.0):
        a, b = idx[(h, "rest", 0.1, 20)], idx[(h, "rest", 0.1, 40)]
        dw = float(np.max(np.abs(b["Wt"] - a["Wt"]) / np.abs(a["Wt"])))
        dv = abs(b["vm"] - a["vm"]) / a["vm"]
        dr = abs(b["Rx30"] - a["Rx30"]) / a["Rx30"]
        dt = b["tpeak"] - a["tpeak"]
        dtr = dt / (a["tpeak"] - 20.0)
        print("h = %g mm: max dW_c/h %.3f %%, d(peak vM) %.2f %%, d(Rx30) %.2f %%, d(T_peak) %+.2f K (%.2f %% of the rise)" % (h, 100 * dw, 100 * dv, 100 * dr, dt, 100 * dtr))
        out.append("%g,%.6f,%.6f,%.6f,%.3f,%.6f" % (h, dw, dv, dr, dt, dtr))
    with open(os.path.join(ROOT, "results", "mesh_sensitivity.csv"), "w") as f:
        f.write("\n".join(out) + "\n")
    print("\nwrote results/study_summary.csv, results/study_summary.md, results/mesh_sensitivity.csv")


if __name__ == "__main__":
    main()
