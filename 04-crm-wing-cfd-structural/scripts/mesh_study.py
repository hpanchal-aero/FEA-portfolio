#!/usr/bin/python3
"""Three-mesh study of the first-order upwind CFD load source at alpha = 0.36 deg (Project 04).

usage: mesh_study.py

For L8 / L9 / L10 (coarse -> fine): cells, wing faces, Cl/Cd/Cm at the common iteration, mean over the common window,
final value, and CalculiX tip-section mean uz / global max uz (from post_static.py).
Convergence analysis after Celik et al. (2008): observed order p, extrapolated value, GCI (safety factor 1.25),
for two refinement-ratio definitions, plus GCI with the formal order p = 1 (first-order upwind).
Non-monotonic or apparently diverging sequences are reported as such.
"""
import math
import os
import re
import subprocess
import sys

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
COMMON, WIN = 1148, 200
RUNS = [
    dict(name="L8", cells=636768, wing=109730, run="simulation/cfd_v1/runs/u_L8_p0p36",
         job="simulation/static_cfd/L8_p0p36/uCRM9_cfd_L8_p0p36"),
    dict(name="L9", cells=1103898, wing=198337, run="simulation/cfd_v1/runs/u_p0p36",
         job="simulation/static_cfd/p0p36_full/uCRM9_cfd_p0p36_full"),
    dict(name="L10", cells=2018605, wing=372363, run="simulation/cfd_v1/runs/u_L10_p0p36",
         job="simulation/static_cfd/L10_p0p36/uCRM9_cfd_L10_p0p36"),
]
FC = "postProcessing/forceCoeffsCompressible/0/forceCoeffs.dat"


def read_fc(path):
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            rows.append([float(p[0]), float(p[3]), float(p[2]), float(p[1])])  # iteration, Cl, Cd, Cm
    return np.array(rows)


def tip_values(job):
    out = subprocess.run([sys.executable, os.path.join(ROOT, "scripts/post_static.py"), os.path.join(ROOT, job)],
                         capture_output=True, text=True, stdin=subprocess.DEVNULL)
    if out.returncode != 0:
        sys.exit(f"ABORT: post_static failed for {job}\n{out.stderr[-500:]}")
    m1 = re.search(r"tip section \(y >= ymax-0.001\):\s+\d+ nodes\s+uz mean\s+([-\d.]+)", out.stdout)
    m2 = re.search(r"max uz\s*=\s*([-\d.]+)", out.stdout)
    if not (m1 and m2):
        sys.exit(f"ABORT: could not parse post_static output for {job}\n{out.stdout[-800:]}")
    return float(m1.group(1)), float(m2.group(1))


def gci(phi, r21, r32, order=None):
    """phi = (fine, medium, coarse).  Returns dict."""
    p1, p2, p3 = phi
    e21, e32 = p2 - p1, p3 - p2
    res = {}
    if e21 == 0 or e32 == 0:
        return {"status": "zero difference"}
    ratio = e32 / e21
    if ratio < 0:
        return {"status": f"non-monotonic (e32/e21 = {ratio:+.3f})"}
    if order is None:
        if abs(e32) < abs(e21):
            return {"status": f"apparent divergence (|e32/e21| = {ratio:.3f} < 1)"}
        p = 2.0
        for _ in range(500):
            q = math.log((r21 ** p - 1.0) / (r32 ** p - 1.0))
            pn = abs(math.log(abs(ratio)) + q) / math.log(r21)
            if abs(pn - p) < 1e-10:
                break
            p = pn
    else:
        p = order
    ext = (r21 ** p * p1 - p2) / (r21 ** p - 1.0)
    ea = abs((p1 - p2) / p1)
    res.update(status="ok", p=p, ext=ext, ea21=100 * ea,
               gci=100 * 1.25 * ea / (r21 ** p - 1.0), eext=100 * abs((ext - p1) / ext))
    return res


data = {}
for r in RUNS:
    path = os.path.join(ROOT, r["run"], FC)
    if not os.path.isfile(path):
        sys.exit(f"ABORT: missing {path}")
    a = read_fc(path)
    t = a[:, 0]
    ic = np.where(np.abs(t - COMMON) < 0.5)[0]
    if len(ic) == 0:
        sys.exit(f"ABORT: iteration {COMMON} not in {path}")
    at = a[ic[0], 1:]
    w = (t > COMMON - WIN) & (t <= COMMON)
    win = a[w, 1:].mean(axis=0)
    last = a[-1, 1:]
    tip, gmax = tip_values(r["job"])
    data[r["name"]] = dict(at=at, win=win, last=last, lastit=int(t[-1]), tip=tip, gmax=gmax, **r)

print(f"Common iteration {COMMON}; common window ({COMMON - WIN}, {COMMON}]; alpha = 0.36 deg, first-order upwind")
print(f"\n{'':4s} {'cells':>10s} {'wingfaces':>10s} {'lastit':>7s} | {'Cl@1148':>8s} {'Cd@1148':>8s} {'Cm@1148':>8s} | "
      f"{'Cl win':>8s} {'Cd win':>8s} {'Cm win':>8s} | {'Cl last':>8s} {'Cd last':>8s} {'Cm last':>8s} | {'tip uz':>7s} {'max uz':>7s}")
for n in ("L8", "L9", "L10"):
    d = data[n]
    print(f"{n:4s} {d['cells']:10d} {d['wing']:10d} {d['lastit']:7d} | "
          + " ".join(f"{v:8.5f}" for v in d["at"]) + " | " + " ".join(f"{v:8.5f}" for v in d["win"]) + " | "
          + " ".join(f"{v:8.5f}" for v in d["last"]) + f" | {d['tip']:7.4f} {d['gmax']:7.4f}")

N = [data[n]["cells"] for n in ("L10", "L9", "L8")]
W = [data[n]["wing"] for n in ("L10", "L9", "L8")]
ratios = {
    "volume  r=(N1/N2)^(1/3)": (( N[0] / N[1]) ** (1 / 3), (N[1] / N[2]) ** (1 / 3)),
    "surface r=(F1/F2)^(1/2)": ((W[0] / W[1]) ** 0.5, (W[1] / W[2]) ** 0.5),
}
for k, (r21, r32) in ratios.items():
    print(f"\nrefinement ratios [{k}]: r21 (L10/L9) = {r21:.4f}   r32 (L9/L8) = {r32:.4f}")

quantities = [("Cl @1148", lambda d: d["at"][0]), ("Cd @1148", lambda d: d["at"][1]), ("Cm @1148", lambda d: d["at"][2]),
              ("Cl window", lambda d: d["win"][0]), ("Cd window", lambda d: d["win"][1]), ("Cm window", lambda d: d["win"][2]),
              ("tip uz", lambda d: d["tip"]), ("max uz", lambda d: d["gmax"])]
for k, (r21, r32) in ratios.items():
    print(f"\n=== GCI, {k} ===")
    print(f"{'quantity':10s} {'fine':>9s} {'medium':>9s} {'coarse':>9s} | {'p_obs':>6s} {'extrap':>9s} {'e_a21%':>7s} {'GCI21%':>7s} | {'GCI21% (p=1)':>12s}")
    for lab, f in quantities:
        phi = (f(data["L10"]), f(data["L9"]), f(data["L8"]))
        g = gci(phi, r21, r32)
        g1 = gci(phi, r21, r32, order=1.0)
        if g["status"] == "ok":
            print(f"{lab:10s} {phi[0]:9.5f} {phi[1]:9.5f} {phi[2]:9.5f} | {g['p']:6.2f} {g['ext']:9.5f} {g['ea21']:7.2f} {g['gci']:7.2f} | "
                  f"{g1['gci']:12.2f}")
        else:
            print(f"{lab:10s} {phi[0]:9.5f} {phi[1]:9.5f} {phi[2]:9.5f} | {g['status']}   (p=1 GCI21 {g1.get('gci', float('nan')):.2f} %)")
