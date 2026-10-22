#!/usr/bin/env python3
"""Verification strip: CalculiX layered / homogeneous shell vs classical laminate theory for one PSHELL property.

  strip_test.py gen <npz> <pid> <layered|homog> <axis 1|2> <outdir>
                [--elem S4|S8R] [--load 100] [--level coarse] [--stack fit|sym|uniform] [--n plies]
  strip_test.py check <outdir> [tol_percent, default 3]

Cantilever strip L = 3 m, b = 0.3 m (60 x 6 elements), clamped at x = 0, tip load P (z) at the tip edge. Beam axis is
material direction <axis>. Expected tip deflection:
  bending  P L^3 / (3 b K),  K = Schur complement of the 6x6 [[A,B],[B,D]] onto the beam-axis curvature, all other strains
           and curvatures free (narrow-strip assumption); ABD of a layered deck is computed from its actual plies.
  shear    P L / (b T G3_<axis>)       (no shear correction factor)
layered stacks (total thickness H of the fit):
  fit      the fitted plies of the npz
  sym      fitted plies symmetrised, q_k <- (q_k + q_(N-1-k))/2   (B = 0)
  uniform  --n identical plies, Q = A_NASTRAN / H   (homogeneous thick plate, B = 0)
homog: membrane-only plate (A of NASTRAN, B = 0, D = A T^2/12), i.e. the representation of the current CalculiX deck.
CalculiX 2.23 allows *SHELL SECTION COMPOSITE only for S8R or S6 elements: layered decks need --elem S8R.
Units: m, N, kg, Pa.
"""
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from ucrm_bdf import ival, num, read_cards  # noqa: E402

ROOT = HERE.parent
L, BW, NX, NY = 3.0, 0.3, 60, 6


def fnum(s, default=0.0):
    return num(s) if s.strip() else default


def load_target(pid, level):
    fem = ROOT / "data" / "ucrm9" / "uCRM9" / "Mesh" / level.capitalize() / "FEM"
    psh, mats, rho = None, {}, {}
    for name, f in read_cards(fem / f"uCRM-9_wingbox_{level}.bdf"):
        if name == "PSHELL" and ival(f[0]) == pid:
            psh = dict(mid1=ival(f[1]), t=num(f[2]), mid2=ival(f[3]), r=fnum(f[4], 1.0), mid3=ival(f[5]),
                       ts=fnum(f[6], 0.8333333), mid4=ival(f[10]))
        elif name == "MAT2":
            m = ival(f[0])
            mats[m] = np.array([[fnum(f[1]), fnum(f[2]), fnum(f[3])],
                                [fnum(f[2]), fnum(f[4]), fnum(f[5])],
                                [fnum(f[3]), fnum(f[5]), fnum(f[6])]])
            rho[m] = fnum(f[7])
    if psh is None:
        sys.exit(f"ABORT: PSHELL {pid} not found")
    t = psh["t"]
    A = t * mats[psh["mid1"]]
    B = t ** 2 * mats[psh["mid4"]]
    D = psh["r"] * t ** 3 / 12.0 * mats[psh["mid2"]]
    G3 = mats[psh["mid3"]]
    return t, A, B, D, rho[psh["mid1"]], G3, psh["ts"]


def k_eff(A, B, D, axis):
    M6 = np.block([[A, B], [B, D]])
    g = [3] if axis == 1 else [4]
    f = [k for k in range(6) if k not in g]
    return float((M6[np.ix_(g, g)] - M6[np.ix_(g, f)] @ np.linalg.solve(M6[np.ix_(f, f)], M6[np.ix_(f, g)]))[0, 0])


def laminate_abd(qpl, H):
    """ABD (3x3 each) of equal-thickness plies; qpl columns Q11, Q22, Q12, Q66."""
    n = len(qpl)
    h = H / n
    z = -H / 2.0 + (np.arange(n) + 0.5) * h
    A, B, D = np.zeros((3, 3)), np.zeros((3, 3)), np.zeros((3, 3))
    for (i, j), c in (((0, 0), 0), ((1, 1), 1), ((0, 1), 2), ((2, 2), 3)):
        a = h * qpl[:, c].sum()
        b = h * (z * qpl[:, c]).sum()
        d = h * ((z * z + h * h / 12.0) * qpl[:, c]).sum()
        A[i, j] = A[j, i] = a
        B[i, j] = B[j, i] = b
        D[i, j] = D[j, i] = d
    return A, B, D


def mesh(elem, P):
    """Return (node lines, element lines, root nodes, tip nodes, load dict {node: force})."""
    if elem == "S4":
        nid = lambda i, j: j * (NX + 1) + i + 1
        nodes = [f"{nid(i, j)}, {i * L / NX:.9e}, {j * BW / NY:.9e}, 0." for j in range(NY + 1) for i in range(NX + 1)]
        els, e = [], 0
        for j in range(NY):
            for i in range(NX):
                e += 1
                els.append(f"{e}, {nid(i, j)}, {nid(i + 1, j)}, {nid(i + 1, j + 1)}, {nid(i, j + 1)}")
        root = [nid(0, j) for j in range(NY + 1)]
        tip = [nid(NX, j) for j in range(NY + 1)]
        loads = {n: P / (NY + 1) for n in tip}
        return nodes, els, root, tip, loads
    ids = {}

    def get(i, j):
        if (i, j) not in ids:
            ids[(i, j)] = len(ids) + 1
        return ids[(i, j)]

    elist = []
    for b in range(NY):
        for a in range(NX):
            c = [(2 * a, 2 * b), (2 * a + 2, 2 * b), (2 * a + 2, 2 * b + 2), (2 * a, 2 * b + 2),
                 (2 * a + 1, 2 * b), (2 * a + 2, 2 * b + 1), (2 * a + 1, 2 * b + 2), (2 * a, 2 * b + 1)]
            elist.append([get(i, j) for i, j in c])
    nodes = [f"{n}, {i * L / (2 * NX):.9e}, {j * BW / (2 * NY):.9e}, 0."
             for (i, j), n in sorted(ids.items(), key=lambda kv: kv[1])]
    els = [f"{k + 1}, " + ", ".join(str(n) for n in el) for k, el in enumerate(elist)]
    root = [ids[(0, j)] for j in range(2 * NY + 1)]
    tip = [ids[(2 * NX, j)] for j in range(2 * NY + 1)]
    loads = {n: 0.0 for n in tip}
    for b in range(NY):
        seg = P / NY
        loads[ids[(2 * NX, 2 * b)]] += seg / 6.0
        loads[ids[(2 * NX, 2 * b + 1)]] += seg * 4.0 / 6.0
        loads[ids[(2 * NX, 2 * b + 2)]] += seg / 6.0
    return nodes, els, root, tip, loads


def gen(argv):
    npz, pid, mode, axis, outdir = argv[0], int(argv[1]), argv[2], int(argv[3]), Path(argv[4])
    opt = lambda key, default: argv[argv.index(key) + 1] if key in argv else default
    elem, P, level = opt("--elem", "S4"), float(opt("--load", 100.0)), opt("--level", "coarse")
    stack, nply = opt("--stack", "fit"), opt("--n", None)
    if mode not in ("layered", "homog") or axis not in (1, 2) or elem not in ("S4", "S8R") \
            or stack not in ("fit", "sym", "uniform"):
        sys.exit("ABORT: mode layered|homog, axis 1|2, elem S4|S8R, stack fit|sym|uniform")
    if mode == "layered" and elem != "S8R":
        sys.exit("ABORT: CalculiX 2.23 allows COMPOSITE only with S8R or S6 elements")
    t, A, B, D, rho, G3, ts = load_target(pid, level)
    z = np.load(npz)
    pids = [int(p) for p in z["pid"]]
    if pid not in pids:
        sys.exit(f"ABORT: pid {pid} not in {npz}")
    k = pids.index(pid)
    H = float(z["H"][k])
    qpl = None
    if mode == "layered":
        qpl = np.array(z["qply"][k], float)
        if stack == "sym":
            qpl = 0.5 * (qpl + qpl[::-1])
        elif stack == "uniform":
            n = int(nply) if nply else len(qpl)
            qpl = np.tile(np.array([A[0, 0], A[1, 1], A[0, 1], A[2, 2]]) / H, (n, 1))
        Af, Bf, Df = laminate_abd(qpl, H)
    else:
        Af, Bf, Df = A, np.zeros((3, 3)), A * t * t / 12.0
    N = len(qpl) if qpl is not None else 0
    K = k_eff(Af, Bf, Df, axis)
    g3 = G3[0, 0] if axis == 1 else G3[1, 1]
    d_b = P * L ** 3 / (3.0 * BW * K)
    d_s = P * L / (BW * t * ts * g3) if ts > 0 else 0.0
    outdir.mkdir(parents=True, exist_ok=True)
    nodes, els, root, tip, loads = mesh(elem, P)
    X = ["*NODE"] + nodes + [f"*ELEMENT, TYPE={elem}, ELSET=EALL"] + els
    X += ["*NSET, NSET=ROOT"] + [str(n) for n in root] + ["*NSET, NSET=TIP"] + [str(n) for n in tip]
    X += ["*ORIENTATION, NAME=OR1, SYSTEM=RECTANGULAR"]
    X += ["1., 0., 0., 0., 1., 0." if axis == 1 else "0., 1., 0., -1., 0., 0."]
    if mode == "layered":
        h = H / N
        for p in range(N):
            q11, q22, q12, q66 = qpl[p]
            det = q11 * q22 - q12 ** 2
            if det <= 0:
                sys.exit(f"ABORT: ply {p} not positive definite")
            e1, e2, nu12 = det / q22, det / q11, q12 / q22
            X += [f"*MATERIAL, NAME=PLY{p + 1}", "*ELASTIC, TYPE=ENGINEERING CONSTANTS",
                  f"{e1:.9e}, {e2:.9e}, {e2:.9e}, {nu12:.9e}, 0.3, 0.3, {q66:.9e}, {G3[0, 0] * t / H:.9e}",
                  f"{G3[1, 1] * t / H:.9e}, 0.", "*DENSITY", f"{rho * t / H:.6f}"]
        X += ["*SHELL SECTION, ELSET=EALL, COMPOSITE"]
        X += [f"{h:.9e}, , PLY{p + 1}, OR1" for p in range(N)]
    else:
        S = np.linalg.inv(A[:2, :2])
        e1, e2, nu12 = 1.0 / S[0, 0] / t, 1.0 / S[1, 1] / t, -S[0, 1] / S[0, 0]
        X += ["*MATERIAL, NAME=HOM", "*ELASTIC, TYPE=ENGINEERING CONSTANTS",
              f"{e1:.9e}, {e2:.9e}, {e2:.9e}, {nu12:.9e}, 0.3, 0.3, {A[2, 2] / t:.9e}, {G3[0, 0]:.9e}",
              f"{G3[1, 1]:.9e}, 0.", "*DENSITY", f"{rho:.6f}",
              "*SHELL SECTION, ELSET=EALL, MATERIAL=HOM, ORIENTATION=OR1", f"{t:.9e}"]
    X += ["*STEP", "*STATIC", "*BOUNDARY", "ROOT, 1, 6", "*CLOAD"]
    X += [f"{n}, 3, {v:.9e}" for n, v in loads.items()]
    X += ["*NODE PRINT, NSET=TIP", "U", "*END STEP"]
    (outdir / "strip.inp").write_text("\n".join(X) + "\n")
    exp = dict(pid=pid, mode=mode, stack=stack if mode == "layered" else "-", axis=axis, elem=elem, N=N, H=H, T=t,
               P=P, K=K, B11=float(Bf[0, 0]), delta_bending=d_b, delta_shear=d_s, delta_total=d_b + d_s,
               load_sum=sum(loads.values()))
    (outdir / "expected.json").write_text(json.dumps(exp, indent=1))
    desc = f"{N} plies, stack {stack}" if mode == "layered" else "homogeneous"
    print(f"wrote {outdir}/strip.inp ({len(els)} {elem} elements, {len(nodes)} nodes, {desc}, "
          f"shell thickness {H if mode == 'layered' else t:.5f} m, load sum {sum(loads.values()):.4f} N)")
    print(f"expected: K = {K:.5e} N m per m width, B11 = {Bf[0, 0]:.4e}, bending {d_b:.5e} m, "
          f"shear {d_s:.5e} m, total {d_b + d_s:.5e} m")


def check(argv):
    d = Path(argv[0])
    tol = float(argv[1]) if len(argv) > 1 else 3.0
    exp = json.loads((d / "expected.json").read_text())
    txt = (d / "strip.dat").read_text().split("\n")
    vals, on = [], False
    for ln in txt:
        if "displacements" in ln:
            on = True
            continue
        if on:
            tok = ln.split()
            if len(tok) == 4 and re.fullmatch(r"\d+", tok[0]):
                vals.append(float(tok[3]))
            elif vals and not tok:
                break
    if not vals:
        sys.exit("ABORT: no displacements found in strip.dat")
    uz = float(np.mean(vals))
    ratio = uz / exp["delta_total"]
    verdict = "PASS" if abs(ratio - 1.0) * 100.0 <= tol else "FAIL"
    print(f"{exp['mode']}/{exp.get('stack', '-')} {exp.get('elem', 'S4')} pid {exp['pid']} axis {exp['axis']} "
          f"N {exp.get('N', 0)}: tip uz (mean of {len(vals)} nodes) = {uz:.5e} m, expected {exp['delta_total']:.5e} m, "
          f"ratio {ratio:.4f} -> {verdict} (tolerance {tol} %)")


if __name__ == "__main__":
    if len(sys.argv) < 3 or sys.argv[1] not in ("gen", "check"):
        sys.exit(__doc__)
    (gen if sys.argv[1] == "gen" else check)(sys.argv[2:])
