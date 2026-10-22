#!/usr/bin/env python3
"""Feasibility test (no solver): can a symmetric-geometry stack of N equal-thickness orthotropic plies reproduce
NASTRAN's smeared A, B, D of every PSHELL property of the uCRM-9 wingbox?

Per property (NASTRAN PSHELL convention, as in section_stiffness.py), in the property's material axes:
  A = T G(MID1),  B = T^2 G(MID4),  D = (12I/T^3) T^3/12 G(MID2)
Laminate with N plies, thickness h = H/N, centres z_k = -H/2 + (k-1/2) h, ply stiffness Q_k (3x3 orthotropic, angle 0):
  A = sum h Q_k,  B = sum h z_k Q_k,  D = sum h (z_k^2 + h^2/12) Q_k      (linear in Q_k)
Components 11, 22, 66: bounded linear least squares with Q_k >= FLOOR * (A_c/H).
Component 12 (variants): |Q12_k| <= sqrt(Q11_k Q22_k) (ply positive definite)
  strict : all three rows (A12, B12, D12) within TOL
  R1     : row weights (1, 0.01, 0.01): A12 within TOL, B12 and D12 only approximately fitted
  R2     : component 12 ignored (the three diagonal components must still pass)
H is scanned from T to HMAX*T (geometric grid of NH points); the smallest H meeting each rule is reported.
Relative error per row = |residual| / max(|target|/scale, 0.01) with scales (A, A H/2, A H^2/12).
Mass is preserved by a uniform ply density rho T / H (not part of the fit).

Usage: /usr/bin/python3 scripts/layered_fit.py [coarse|medium|fine] [N list, default 25] [FLOOR, default 1e-3]
"""
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.optimize import lsq_linear

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from ucrm_bdf import ival, num, read_cards  # noqa: E402

root = HERE.parent
level = sys.argv[1] if len(sys.argv) > 1 else "coarse"
n_list = [int(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else [25]
FLOOR = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-3
TOL, NH, HMAX = 1e-3, 50, 25.0
fem = root / "data" / "ucrm9" / "uCRM9" / "Mesh" / level.capitalize() / "FEM"


def fnum(s, default=0.0):
    return num(s) if s.strip() else default


pshell, mat2, count = {}, {}, Counter()
for name, f in read_cards(fem / f"uCRM-9_wingbox_{level}.bdf"):
    if name == "CQUADR":
        count[ival(f[1])] += 1
    elif name == "PSHELL":
        pshell[ival(f[0])] = dict(mid1=ival(f[1]), t=num(f[2]), mid2=ival(f[3]), r=fnum(f[4], 1.0), mid4=ival(f[10]))
    elif name == "MAT2":
        mat2[ival(f[0])] = np.array([[fnum(f[1]), fnum(f[2]), fnum(f[3])],
                                     [fnum(f[2]), fnum(f[4]), fnum(f[5])],
                                     [fnum(f[3]), fnum(f[5]), fnum(f[6])]])

targets = {}
for pid, p in pshell.items():
    t = p["t"]
    A = t * mat2[p["mid1"]]
    B = t ** 2 * mat2[p["mid4"]]
    D = p["r"] * t ** 3 / 12.0 * mat2[p["mid2"]]
    targets[pid] = (t, A, B, D)
print(f"level {level}: {len(pshell)} properties, {sum(count.values())} elements, FLOOR {FLOOR:g}, TOL {TOL:g}")


def solve(At, Bt, Dt, sA, H, N, lb, ub, w=(1.0, 1.0, 1.0)):
    """Return ply stiffnesses [Pa], relative errors of rows (A, B, D) and the fitted physical (A, B, D)."""
    h = H / N
    z = -H / 2.0 + (np.arange(N) + 0.5) * h
    q0 = sA / H
    sB, sD = sA * H / 2.0, sA * H * H / 12.0
    M = np.vstack([h * q0 * np.ones(N) / sA, h * q0 * z / sB, h * q0 * (z * z + h * h / 12.0) / sD])
    tgt = np.array([At / sA, Bt / sB, Dt / sD])
    lbv = np.full(N, lb) if np.isscalar(lb) else lb / q0
    ubv = np.full(N, ub) if np.isscalar(ub) else ub / q0
    wv = np.asarray(w, float)
    res = lsq_linear(M * wv[:, None], tgt * wv, bounds=(lbv, ubv), method="bvls")
    r = M @ res.x - tgt
    rel = np.abs(r) / np.maximum(np.abs(tgt), 0.01)
    return res.x * q0, rel, (M @ res.x) * np.array([sA, sB, sD])


def fit_property(pid, N):
    t, A, B, D = targets[pid]
    s12 = float(np.sqrt(A[0, 0] * A[1, 1]))
    dscale = float(np.sqrt(D[0, 0] * D[1, 1]))
    out = dict(R2=None, strict=None, R1=None, passing={"11": [], "22": [], "66": []}, rows_at_Hlo=None)
    for H in t * np.geomspace(1.0, HMAX, NH):
        Q, err = {}, {}
        for key, (i, j) in (("11", (0, 0)), ("22", (1, 1)), ("66", (2, 2))):
            q, rel, _ = solve(A[i, j], B[i, j], D[i, j], A[i, j], H, N, FLOOR, np.inf)
            Q[key], err[key] = q, float(rel.max())
            if err[key] <= TOL:
                out["passing"][key].append(H / t)
        if max(err.values()) > TOL:
            continue
        info = dict(ratio=H / t, contrast=float(Q["11"].max() / Q["11"].min()))
        if out["R2"] is None:
            out["R2"] = info
        ub = (1.0 - 1e-6) * np.sqrt(Q["11"] * Q["22"])
        _, rel_s, _ = solve(A[0, 1], B[0, 1], D[0, 1], s12, H, N, -ub, ub)
        if out["rows_at_Hlo"] is None:
            out["rows_at_Hlo"] = (H / t, rel_s)
        if out["strict"] is None and rel_s.max() <= TOL:
            out["strict"] = dict(info)
        _, rel_r, fit_r = solve(A[0, 1], B[0, 1], D[0, 1], s12, H, N, -ub, ub, w=(1.0, 0.01, 0.01))
        if out["R1"] is None and rel_r[0] <= TOL:
            out["R1"] = dict(info, dD=abs(fit_r[2] - D[0, 1]) / dscale, dB=abs(fit_r[1] - B[0, 1]) / (s12 * H / 2.0))
    return out


def span(vals):
    return "none" if not vals else f"[{min(vals):.1f},{max(vals):.1f}]"


total_el = sum(count.values())
for N in n_list:
    res = {pid: fit_property(pid, N) for pid in sorted(pshell)}
    print(f"\n=== N = {N} plies ===")
    for rule in ("strict", "R1", "R2"):
        ok = [p for p in res if res[p][rule]]
        el = sum(count[p] for p in ok)
        line = f"  {rule:6s} feasible {len(ok):3d}/{len(res)} properties, {100.0 * el / total_el:5.1f} % of elements"
        if ok:
            r = np.array([res[p][rule]["ratio"] for p in ok])
            line += f"   H/T median {np.median(r):.2f} p90 {np.percentile(r, 90):.2f} max {r.max():.2f}"
        print(line)
    r1 = [p for p in res if res[p]["R1"]]
    if r1:
        dD = np.array([res[p]["R1"]["dD"] for p in r1])
        dB = np.array([res[p]["R1"]["dB"] for p in r1])
        print(f"  R1 residuals  D12 error/sqrt(D11 D22): median {np.median(dD):.3g} max {dD.max():.3g};  "
              f"B12 error/(sqrt(A11 A22) H/2): median {np.median(dB):.3g} max {dB.max():.3g}")
    rows = [(p, res[p]["rows_at_Hlo"]) for p in res if res[p]["rows_at_Hlo"] is not None]
    cnt = Counter()
    for p, (hlo, rel) in rows:
        if rel.max() > TOL:
            cnt["ABD"[int(np.argmax(rel))]] += count[p]
    print(f"  at the thinnest H where 11, 22, 66 pass: component-12 row that fails worst (elements): {dict(cnt)}")
    for p in (1, 94, 95):
        q = res[p]
        tag = " ".join(f"{k}={'H/T %.2f' % q[k]['ratio'] if q[k] else 'infeasible'}" for k in ("strict", "R1", "R2"))
        extra = ""
        if q["rows_at_Hlo"] is not None:
            hlo, rel = q["rows_at_Hlo"]
            extra = f"  12-row errors at H/T {hlo:.2f}: A {rel[0]:.3g} B {rel[1]:.3g} D {rel[2]:.3g}"
        print(f"  pid {p:3d}: {tag}{extra}")
    bad = sorted((p for p in res if not res[p]["R2"]), key=lambda p: -count[p])[:6]
    print(f"  properties infeasible even under R2: {sum(1 for p in res if not res[p]['R2'])}")
    for p in bad:
        pa = res[p]["passing"]
        print(f"   pid {p:4d} ({count[p]:5d} el): 11:{span(pa['11'])} 22:{span(pa['22'])} 66:{span(pa['66'])}")
