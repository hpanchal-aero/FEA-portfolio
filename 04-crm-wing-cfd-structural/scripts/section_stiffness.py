#!/usr/bin/env python3
"""H2 test: what do the neglected membrane-bending coupling (MID4) and bending (MID2) matrices do to the
bending stiffness of the uCRM-9 wingbox?  Beam-theory estimate from the BDF data, no solver involved.

Per property (NASTRAN PSHELL convention): A = t G(MID1), B = t^2 G(MID4), D = (12I/T^3) t^3/12 G(MID2).
Per span station y0 the wingbox is cut by the plane y = y0; every cut panel segment (length l, |u| = length of
the projection of the y axis on the panel) contributes
    EI = sum  l |u|^3 [ a (z - zn)^2 + 2 b (z - zn) n_z + d n_z^2 ]
with (a, b, d) the reduced stiffness of the pair (membrane strain, curvature) along u after eliminating the
other components of the panel in the frame (u, v) that is aligned with the beam axis:
  free       : all other strains and curvatures free (their resultants are zero)
  restrained : transverse curvature and twist held at zero, other strains free
Model M (the CalculiX deck): membrane stiffness only (b = 0) plus the thin shell's own d = a t^2/12.
Cumulative stages per variant: 1 = a only, 2 = a and b, 3 = a, b and d (full ABD).
Tip deflection by the unit-load method under the nodal z forces of the v2 deck; only the ratio
delta_M / delta_X is interpreted.

Usage: /usr/bin/python3 scripts/section_stiffness.py [coarse|medium|fine] [n_stations]
"""
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from ucrm_bdf import ival, num, read_cards  # noqa: E402

level = sys.argv[1] if len(sys.argv) > 1 else "coarse"
NST = int(sys.argv[2]) if len(sys.argv) > 2 else 80
if level not in ("coarse", "medium", "fine"):
    sys.exit("level must be coarse, medium or fine")
root = HERE.parent
fem = root / "data" / "ucrm9" / "uCRM9" / "Mesh" / level.capitalize() / "FEM"
deck = root / "simulation" / "static_v2" / level / f"uCRM9_{level}_v2.inp"
if not deck.exists():
    sys.exit(f"ABORT: {deck} not found (run ucrm_to_ccx.py {level} --all-loads first)")


def fnum(s, default=0.0):
    return num(s) if s.strip() else default


grids, quads, pshell, mat2, cord = {}, [], {}, {}, {}
comps = defaultdict(set)
for name, f in read_cards(fem / f"uCRM-9_wingbox_{level}.bdf"):
    if name == "GRID":
        grids[ival(f[0])] = np.array([num(f[2]), num(f[3]), num(f[4])])
    elif name == "CQUADR":
        quads.append((ival(f[0]), ival(f[1]), [ival(f[i]) for i in range(2, 6)], ival(f[6])))
    elif name == "PSHELL":
        pshell[ival(f[0])] = dict(mid1=ival(f[1]), t=num(f[2]), mid2=ival(f[3]), r=fnum(f[4], 1.0), mid4=ival(f[10]))
    elif name == "MAT2":
        mat2[ival(f[0])] = np.array([[fnum(f[1]), fnum(f[2]), fnum(f[3])],
                                     [fnum(f[2]), fnum(f[4]), fnum(f[5])],
                                     [fnum(f[3]), fnum(f[5]), fnum(f[6])]])
    elif name == "CORD2R":
        cord[ival(f[0])] = [np.array([num(f[i]) for i in range(j, j + 3)]) for j in (2, 5, 8)]
    elif name == "SPC":
        i = 1
        while i + 1 < len(f):
            if f[i].strip() and ival(f[i]) > 0:
                comps[ival(f[i])].update(c for c in f[i + 1].strip() if c in "123456")
            i += 3

for pid, p in pshell.items():
    for key in ("mid1", "mid2", "mid4"):
        if p[key] not in mat2:
            sys.exit(f"ABORT: property {pid} references missing material {p[key]} ({key})")


def mat_x(cid):
    A, B, C = cord[cid]
    z = (B - A) / np.linalg.norm(B - A)
    x = (C - A) - np.dot(C - A, z) * z
    return x / np.linalg.norm(x)


def schur(K, g, f):
    Kgg = K[np.ix_(g, g)]
    if not f:
        return Kgg
    Kgf = K[np.ix_(g, f)]
    Kff = K[np.ix_(f, f)]
    return Kgg - Kgf @ np.linalg.solve(Kff, Kgf.T)


# ---- per-property ABD (scaled: strains (eps, t kappa), resultants (N, M/t)) -----------------------
prop = {}
n_psd_bad = n_ill = 0
for pid, p in pshell.items():
    t = p["t"]
    A = t * mat2[p["mid1"]]
    B = t ** 2 * mat2[p["mid4"]]
    D = p["r"] * t ** 3 / 12.0 * mat2[p["mid2"]]
    Ms = np.block([[A, B / t], [B / t, D / t ** 2]])
    Ms = 0.5 * (Ms + Ms.T)
    w = np.linalg.eigvalsh(Ms)
    if w[0] < -1e-9 * w[-1]:
        n_psd_bad += 1
    if np.linalg.cond(Ms) > 1e12:
        n_ill += 1
    prop[pid] = dict(t=t, Ms=Ms, Ainv=np.linalg.inv(A), Cs=np.linalg.pinv(Ms, rcond=1e-12))
print(f"properties {len(pshell)}: ABD not positive semi-definite for {n_psd_bad}; ill-conditioned {n_ill}")

# ---- per-element beam data -----------------------------------------------------------------------
Ps, nvs, nzs, uus = [], [], [], []
aM, dM, aF, bF, dF, aR, bR, dR = ([] for _ in range(8))
skip = Counter()
max_rel_diff = 0.0
for eid, pid, ns, mcid in quads:
    P = np.array([grids[n] for n in ns])
    n = np.cross(P[2] - P[0], P[3] - P[1])
    nn = np.linalg.norm(n)
    if nn < 1e-12:
        skip["degenerate"] += 1
        continue
    n = n / nn
    if abs(n[1]) >= 0.5:
        skip["rib (normal mostly spanwise)"] += 1
        continue
    xc = mat_x(mcid)
    x1 = xc - np.dot(xc, n) * n
    if np.linalg.norm(x1) < 1e-6:
        skip["material axis normal to panel"] += 1
        continue
    x1 /= np.linalg.norm(x1)
    x2 = np.cross(n, x1)
    u = np.array([0.0, 1.0, 0.0]) - n[1] * n
    uu = np.linalg.norm(u)
    if uu < 1e-6:
        skip["u degenerate"] += 1
        continue
    uh = u / uu
    c, s = float(np.dot(uh, x1)), float(np.dot(uh, x2))
    pr = prop[pid]
    t = pr["t"]
    m = np.array([c * c, s * s, c * s])
    a_m = 1.0 / float(m @ pr["Ainv"] @ m)
    # frame (u, v): rotate every 3x3 block with T (stress-like Voigt, same for resultants)
    T3 = np.array([[c * c, s * s, 2 * c * s], [s * s, c * c, -2 * c * s], [-c * s, c * s, c * c - s * s]])
    T6 = np.zeros((6, 6))
    T6[:3, :3] = T3
    T6[3:, 3:] = T3
    Mr = T6 @ pr["Ms"] @ T6.T
    Kf = schur(Mr, [0, 3], [1, 2, 4, 5])
    Kr = schur(Mr, [0, 3], [1, 2])
    # consistency with the v1 formulation (inverse of the projected compliance)
    Pm = np.zeros((6, 2))
    Pm[:3, 0] = m
    Pm[3:, 1] = m
    Kold = np.linalg.inv(Pm.T @ pr["Cs"] @ Pm)
    max_rel_diff = max(max_rel_diff, float(np.max(np.abs(Kold - Kf)) / np.max(np.abs(Kold))))
    Ps.append(P)
    nvs.append(n)
    nzs.append(n[2])
    uus.append(uu)
    aM.append(a_m)
    dM.append(a_m * t * t / 12.0)
    aF.append(Kf[0, 0])
    bF.append(t * Kf[0, 1])
    dF.append(t * t * Kf[1, 1])
    aR.append(Kr[0, 0])
    bR.append(t * Kr[0, 1])
    dR.append(t * t * Kr[1, 1])

Ps, nvs = np.array(Ps), np.array(nvs)
nzs, uus, aM, dM, aF, bF, dF, aR, bR, dR = (np.array(v) for v in (nzs, uus, aM, dM, aF, bF, dF, aR, bR, dR))
ymin, ymax = Ps[:, :, 1].min(axis=1), Ps[:, :, 1].max(axis=1)
print(f"elements used {len(Ps)} of {len(quads)}; skipped {dict(skip)}")
print(f"free variant vs v1 formulation: max relative difference {max_rel_diff:.2e}")
print(f"a_F/a_M median {np.median(aF / aM):.4f} [{np.min(aF / aM):.4f}, {np.max(aF / aM):.4f}]   "
      f"a_R/a_M median {np.median(aR / aM):.4f} [{np.min(aR / aM):.4f}, {np.max(aR / aM):.4f}]")

# ---- are the element normals pointing into the box? ----------------------------------------------
cen = Ps.mean(axis=1)
bins = np.floor(cen[:, 1] / 0.25).astype(int)
inward = counted = 0
for bnum in np.unique(bins):
    msk = bins == bnum
    mx, mz = cen[msk, 0].mean(), cen[msk, 2].mean()
    dx, dz = mx - cen[msk, 0], mz - cen[msk, 2]
    far = np.hypot(dx, dz) > 0.02
    dot = nvs[msk, 0] * dx + nvs[msk, 2] * dz
    inward += int(np.sum((dot > 0) & far))
    counted += int(np.sum(far))
print(f"element normals pointing into the box: {inward} of {counted} ({100.0 * inward / max(counted, 1):.1f} %)")

# ---- supports and loads --------------------------------------------------------------------------
clamped = [n for n, c in comps.items() if set("123456") <= c]
yc = np.array([grids[n][1] for n in clamped])
y_root = float(np.median(yc))
y_tip = max(g[1] for g in grids.values())
cl, kw = [], None
for ln in deck.read_text().splitlines():
    s_ = ln.strip()
    if s_.startswith("*"):
        kw = s_.split(",")[0].upper()
        continue
    if kw == "*CLOAD" and s_:
        t_ = [x.strip() for x in s_.split(",")]
        if int(t_[1]) == 3:
            cl.append((int(t_[0]), float(t_[2])))
ycl = np.array([grids[n][1] for n, _ in cl])
Fz = np.array([v for _, v in cl])
print(f"loads: sum Fz = {Fz.sum():.4e} N; beam root y = {y_root:.3f} m, tip y = {y_tip:.3f} m")

# ---- sections ------------------------------------------------------------------------------------
dy = (y_tip - y_root) / NST
ys = y_root + (np.arange(NST) + 0.5) * dy + 0.0137
cutskip = 0


def ei(w, zr, nz, a, b, d):
    zn = np.sum(w * (a * zr + b * nz)) / np.sum(w * a)
    return np.sum(w * (a * (zr - zn) ** 2 + 2 * b * (zr - zn) * nz + d * nz ** 2))


NAMES = ["M", "F1", "F2", "F3", "R1", "R2", "R3"]


def section(y0):
    global cutskip
    ws, zs, ks = [], [], []
    for k in np.nonzero((ymin < y0) & (ymax > y0))[0]:
        P = Ps[k]
        pts = []
        for i in range(4):
            p, q = P[i], P[(i + 1) % 4]
            if (p[1] - y0) * (q[1] - y0) < 0:
                pts.append(p + (y0 - p[1]) / (q[1] - p[1]) * (q - p))
        if len(pts) != 2:
            cutskip += 1
            continue
        ws.append(np.linalg.norm(pts[0] - pts[1]) * uus[k] ** 3)
        zs.append(0.5 * (pts[0][2] + pts[1][2]))
        ks.append(k)
    if not ks:
        return [np.nan] * 7
    w, zr, k = np.array(ws), np.array(zs), np.array(ks)
    nz, zero = nzs[k], np.zeros(len(k))
    return [ei(w, zr, nz, aM[k], zero, dM[k]),
            ei(w, zr, nz, aF[k], zero, dM[k]), ei(w, zr, nz, aF[k], bF[k], dM[k]), ei(w, zr, nz, aF[k], bF[k], dF[k]),
            ei(w, zr, nz, aR[k], zero, dM[k]), ei(w, zr, nz, aR[k], bR[k], dM[k]), ei(w, zr, nz, aR[k], bR[k], dR[k])]


EI = np.array([section(y0) for y0 in ys])           # (NST, 7)
ok = np.all(np.isfinite(EI) & (EI > 0), axis=1)
print(f"\nstations {NST}, usable {int(ok.sum())}, cuts with an unexpected crossing count {cutskip}")
Mk = np.array([np.sum(Fz[ycl > y] * (ycl[ycl > y] - y)) for y in ys])
g = Mk * (y_tip - ys)
delta = {nm: float(np.sum(g[ok] / EI[ok, j]) * dy) for j, nm in enumerate(NAMES)}

labels = {"F1": "free        1: a only", "F2": "free        2: + eccentricity b", "F3": "free        3: + intrinsic d (= v1)",
          "R1": "restrained  1: a only", "R2": "restrained  2: + eccentricity b", "R3": "restrained  3: + intrinsic d"}
print(f"\nbeam tip deflection, model M (CalculiX deck): {delta['M']:.4f} m  (CalculiX v2 coarse: 2.3368 m)")
print(f"{'variant':36s} {'EI_X/EI_M median [min, max]':>32s} {'delta_M/delta_X':>16s}")
for j, nm in enumerate(NAMES[1:], start=1):
    r = EI[ok, j] / EI[ok, 0]
    print(f"{labels[nm]:36s} {np.median(r):8.4f} [{r.min():.4f}, {r.max():.4f}] {delta['M'] / delta[nm]:16.4f}")
print("\nprediction (written before the run): free stage 3 reproduces 0.9089 (v1); normals > 95 % inward; "
      "restrained stage 3 about 0.95 to 0.98")
print("observed CalculiX / NASTRAN tip ratio: 0.926 (tip-section mean) to 0.937 (tip max)")
print("caveat: plane sections, no shear lag or shear deformation, uniaxial stress, root at the clamped nodes; "
      "the absolute beam deflection is not trusted")
