#!/usr/bin/env python3
"""Export the best-effort layered-shell fit of the uCRM-9 smeared A, B, D and report its membrane errors.

Per property (symmetric stack of N equal-thickness plies at 0 degrees, see layered_fit.py):
  components 11, 22, 66: fitted at the smallest H where all three pass TOL (floor FLOOR on ply stiffness)
  component 12: strict (all three rows within TOL) if feasible at some H; otherwise at the smallest diagonal-feasible H
                with row weights (1, 0.01, 0.01) (A12 prioritised), |Q12_k| <= sqrt(Q11_k Q22_k)
Writes an .npz with pid, H, strict flag, per-ply Q (Q11, Q22, Q12, Q66) and the fitted 3x3 A, B, D.
Reports element-weighted errors of the fitted membrane law (A12 and effective E1, E2 from the 2x2 plane-stress inverse).

Usage: /usr/bin/python3 scripts/layered_abd_export.py [coarse|medium|fine] [N, default 25] [FLOOR, default 1e-3] [out.npz]
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
N = int(sys.argv[2]) if len(sys.argv) > 2 else 25
FLOOR = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-3
out_path = Path(sys.argv[4]) if len(sys.argv) > 4 else root / "simulation" / "layered_v1" / f"layered_abd_{level}.npz"
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
    targets[pid] = (t, t * mat2[p["mid1"]], t ** 2 * mat2[p["mid4"]], p["r"] * t ** 3 / 12.0 * mat2[p["mid2"]])


def solve(At, Bt, Dt, sA, H, lb, ub, w=(1.0, 1.0, 1.0)):
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
    rel = np.abs(M @ res.x - tgt) / np.maximum(np.abs(tgt), 0.01)
    return res.x * q0, rel


def laminate(q, H):
    h = H / N
    z = -H / 2.0 + (np.arange(N) + 0.5) * h
    a, b, d = h * q.sum(), h * (z * q).sum(), h * ((z * z + h * h / 12.0) * q).sum()
    return a, b, d


def fit_property(pid):
    t, A, B, D = targets[pid]
    s12 = float(np.sqrt(A[0, 0] * A[1, 1]))
    fallback = None
    for H in t * np.geomspace(1.0, HMAX, NH):
        Q, ok = {}, True
        for key, (i, j) in (("11", (0, 0)), ("22", (1, 1)), ("66", (2, 2))):
            q, rel = solve(A[i, j], B[i, j], D[i, j], A[i, j], H, FLOOR, np.inf)
            Q[key] = q
            ok = ok and rel.max() <= TOL
        if not ok:
            continue
        ub = (1.0 - 1e-6) * np.sqrt(Q["11"] * Q["22"])
        q12, rel = solve(A[0, 1], B[0, 1], D[0, 1], s12, H, -ub, ub)
        if rel.max() <= TOL:
            return H, True, Q["11"], Q["22"], q12, Q["66"]
        if fallback is None:
            q12w, _ = solve(A[0, 1], B[0, 1], D[0, 1], s12, H, -ub, ub, w=(1.0, 0.01, 0.01))
            fallback = (H, False, Q["11"], Q["22"], q12w, Q["66"])
    return fallback


pids = sorted(pshell)
H_arr, strict, qply, A_fit, B_fit, D_fit = [], [], [], [], [], []
for pid in pids:
    r = fit_property(pid)
    if r is None:
        sys.exit(f"ABORT: no diagonal-feasible H for pid {pid}")
    H, st, q11, q22, q12, q66 = r
    H_arr.append(H)
    strict.append(st)
    qply.append(np.stack([q11, q22, q12, q66], axis=1))
    Af, Bf, Df = np.zeros((3, 3)), np.zeros((3, 3)), np.zeros((3, 3))
    for M_, idx in ((Af, 0), (Bf, 1), (Df, 2)):
        pass
    for (i, j), q in (((0, 0), q11), ((1, 1), q22), ((0, 1), q12), ((2, 2), q66)):
        a, b, d = laminate(q, H)
        Af[i, j] = Af[j, i] = a
        Bf[i, j] = Bf[j, i] = b
        Df[i, j] = Df[j, i] = d
    A_fit.append(Af)
    B_fit.append(Bf)
    D_fit.append(Df)

H_arr, strict = np.array(H_arr), np.array(strict)
out_path.parent.mkdir(parents=True, exist_ok=True)
np.savez(out_path, pid=np.array(pids), H=H_arr, strict=strict, qply=np.array(qply), A=np.array(A_fit),
         B=np.array(B_fit), D=np.array(D_fit), N=N, floor=FLOOR, t=np.array([targets[p][0] for p in pids]))
print(f"wrote {out_path}  ({len(pids)} properties, N = {N}, FLOOR {FLOOR:g})")

w = np.array([count[p] for p in pids], float)


def wstat(name, v):
    o = np.argsort(v)
    vs, ws = v[o], w[o]
    c = np.cumsum(ws) / ws.sum()
    q = lambda x: vs[min(np.searchsorted(c, x), len(vs) - 1)]
    print(f"  {name:44s} median {q(0.5):9.3g}  p90 {q(0.9):9.3g}  max {vs.max():9.3g}")


print(f"strict (all rows of component 12 within TOL): {strict.sum()}/{len(pids)} properties, "
      f"{100.0 * w[strict].sum() / w.sum():.1f} % of elements")
ratio = H_arr / np.array([targets[p][0] for p in pids])
wstat("H/T used", ratio)
cont = np.array([q[:, 0].max() / q[:, 0].min() for q in qply])
pdm = np.array([np.min(1.0 - q[:, 2] ** 2 / (q[:, 0] * q[:, 1])) for q in qply])
wstat("ply contrast Q11 max/min", cont)
print(f"  worst ply (1 - Q12^2/(Q11 Q22)): {pdm.min():.3g}")
e1, e2, a12 = [], [], []
for k, pid in enumerate(pids):
    t, A, B, D = targets[pid]
    S_t = np.linalg.inv(A[:2, :2])
    S_f = np.linalg.inv(A_fit[k][:2, :2])
    e1.append(abs(S_t[0, 0] / S_f[0, 0] - 1.0))
    e2.append(abs(S_t[1, 1] / S_f[1, 1] - 1.0))
    a12.append(abs(A_fit[k][0, 1] - A[0, 1]) / np.sqrt(A[0, 0] * A[1, 1]))
print("element-weighted errors of the fitted membrane law (fit vs NASTRAN):")
wstat("A12 error / sqrt(A11 A22)", np.array(a12))
wstat("effective E1 (1/S11 of plane-stress inverse)", np.array(e1))
wstat("effective E2 (1/S22)", np.array(e2))
sel = [pids.index(p) for p in (1, 94, 95)]
for k in sel:
    print(f"  pid {pids[k]:3d}: strict {strict[k]}  H/T {ratio[k]:.2f}  A12 err {a12[k]:.3g}  E1 err {e1[k]:.3g}  E2 err {e2[k]:.3g}")
