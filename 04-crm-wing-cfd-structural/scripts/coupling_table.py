#!/usr/bin/python3
"""Consolidated coupling results table for Project 04.

usage: coupling_table.py

For each structural case it runs post_static.py (tip uz) and check_equilibrium.py (aerodynamic load A, gravity W, closure),
reads the CFD lift coefficient at the load-source iteration from forceCoeffs.dat and the elastic twist of the outermost 1 m
section from results/section_kinematics_<case>.csv.  Nothing is typed in by hand except which job belongs to which CFD run.
Prints the table and writes results/coupling_summary.csv and results/coupling_summary.md.
"""
import csv
import os
import re
import subprocess
import sys

RUNS = "simulation/cfd_v1/runs"
QINF = 10955.9      # Pa, dynamic pressure at the flight condition (Mach 0.85, 216.65 K, 21,662.7 Pa)
AREF = 191.87       # m2, half-wing reference area used by the force coefficients
STC = "simulation/static_cfd"
CASES = [
    dict(name="dataset loads", src="dataset reference loads (wing-body-tail configuration)",
         job="simulation/static_v3/coarse_tr_f2/uCRM9_tr_f2", cl=None, kin="dataset_trf2"),
    dict(name="one-way, L9", src="undeflected wing, alpha 0.36 deg, iteration 1500",
         job=f"{STC}/p0p36_full/uCRM9_cfd_p0p36_full", cl=("u_p0p36", 1500), kin="p0p36_full"),
    dict(name="one-way, L10", src="undeflected wing, alpha 0.36 deg, L10 mesh, iteration 1148",
         job=f"{STC}/L10_p0p36/uCRM9_cfd_L10_p0p36", cl=("u_L10_p0p36", 1148), kin="L10_p0p36"),
    dict(name="corrected, same alpha", src="deflected wing, alpha 0.36 deg, iteration 1500",
         job=f"{STC}/defl_p0p36_t1500_pb/uCRM9_cfd_defl_p0p36_t1500_pb", cl=("u_defl_p0p36", 1500),
         kin="defl_p0p36_t1500_pb"),
    dict(name="corrected, alpha 2.07", src="deflected wing, alpha 2.07 deg, iteration 1300",
         job=f"{STC}/defl_a2p07_t1300_pb/uCRM9_cfd_defl_a2p07_t1300_pb", cl=("u_defl_a2p07", 1300),
         kin="defl_a2p07_t1300_pb"),
    dict(name="corrected, CL 0.5", src="deflected wing, interpolated linearly in lift between the two deflected runs",
         job=f"{STC}/defl_blend_cl0p5_pb/uCRM9_cfd_defl_blend_cl0p5_pb", cl=("fixed", 0.5),
         kin="defl_blend_cl0p5_pb"),
]


def cfd_cl(spec):
    if spec is None:
        return None
    if spec[0] == "fixed":
        return spec[1]
    tag, it = spec
    path = f"{RUNS}/{tag}/postProcessing/forceCoeffsCompressible/0/forceCoeffs.dat"
    if not os.path.exists(path):
        sys.exit(f"ABORT: {path} not found")
    for line in open(path):
        t = line.split()
        if line.startswith("#") or len(t) != 6:
            continue
        try:
            if float(t[0]) == it:
                return float(t[3])
        except ValueError:
            continue
    sys.exit(f"ABORT: iteration {it} not found in {path}")


def call_script(script, job):
    r = subprocess.run([sys.executable, f"scripts/{script}", job], capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"ABORT: {script} {job} failed (exit {r.returncode}):\n{r.stderr[-800:]}")
    return r.stdout


def grab(pattern, text, what, job):
    m = re.search(pattern, text, re.M)
    if not m:
        sys.exit(f"ABORT: could not find {what} in the output for {job}")
    return float(m.group(1))


def twist(kin):
    path = f"results/section_kinematics_{kin}.csv"
    if not os.path.exists(path):
        sys.exit(f"ABORT: {path} not found")
    best = None
    with open(path) as f:
        for r in csv.DictReader(f):
            if best is None or float(r["y_c"]) > float(best["y_c"]):
                best = r
    return float(best["y_c"]), float(best["twist_deg"])


FZ = r"Fz=([+-][\d.]+e[+-]\d+) N"
rows = []
for c in CASES:
    job = c["job"]
    missing = [job + e for e in (".inp", ".dat", ".frd") if not os.path.exists(job + e)]
    if missing:
        sys.exit(f"ABORT: missing files: {missing}")
    ps = call_script("post_static.py", job)
    ce = call_script("check_equilibrium.py", job)
    tip = grab(r"tip section \(y >= ymax-0\.001\):\s+\d+ nodes\s+uz mean\s+(-?[\d.]+)", ps, "tip uz", job)
    a_fz = grab(r"CLOAD sum A\s*:.*?" + FZ, ce, "CLOAD sum A", job)
    w_fz = grab(r"gravity sum W[^:]*:.*?" + FZ, ce, "gravity sum W", job)
    clo = grab(r"closure R \+ A \+ W\s*:.*?" + FZ, ce, "closure", job)
    yc, tw = twist(c["kin"])
    rows.append(dict(case=c["name"], src=c["src"], cl=cfd_cl(c["cl"]), a=a_fz / 1e3, w=w_fz / 1e3,
                     ceq=a_fz / (QINF * AREF), tip=tip, tw=tw, yc=yc, clo=clo))

ref = rows[0]["tip"]
head = ["case", "load source", "CFD CL", "aero Fz (kN)", "Fz/(q*Aref)", "gravity Fz (kN)", "tip uz (m)",
        "tip uz vs dataset (%)", "twist, outer section (deg)", "Fz closure (N)"]
table = []
for r in rows:
    table.append([r["case"], r["src"], "-" if r["cl"] is None else f"{r['cl']:.4f}", f"{r['a']:.1f}", f"{r['ceq']:.4f}",
                  f"{r['w']:.1f}", f"{r['tip']:.4f}", "-" if r is rows[0] else f"{100.0 * (r['tip'] / ref - 1.0):+.1f}",
                  f"{r['tw']:+.3f}", f"{r['clo']:+.2f}"])
print(" | ".join(head))
for t in table:
    print(" | ".join(t))
notes = [
    f"Outer section = the section centred at y = {rows[1]['yc']:.2f} m (outermost 1 m bin); twist = -b of the plane fit uz = a + b (x - xbar) + c (y - ybar), "
    "negative = nose-down washout; not validated against the NASTRAN reference.",
    f"Fz/(q*Aref) uses q = {QINF} Pa and Aref = {AREF} m2; for the CFD cases it should equal the CFD CL within the small drag/inclination contribution.",
    "Gravity and point masses are the same in all cases. The dataset row uses the dataset's own wing-body-tail reference loads, the other rows wing-alone CFD loads "
    "applied to the same coarse CalculiX deck (static tip -7.41% vs the NASTRAN reference, documented): differences to the dataset row mix the load difference with that bias.",
]
for n in notes:
    print("note:", n)
with open("results/coupling_summary.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(head)
    w.writerows(table)
with open("results/coupling_summary.md", "w") as f:
    f.write("| " + " | ".join(head) + " |\n|" + "---|" * len(head) + "\n")
    for t in table:
        f.write("| " + " | ".join(t) + " |\n")
    f.write("\n" + "\n".join("- " + n for n in notes) + "\n")
print("saved results/coupling_summary.csv and results/coupling_summary.md")
