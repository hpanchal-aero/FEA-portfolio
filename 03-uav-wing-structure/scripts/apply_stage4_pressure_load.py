#!/usr/bin/env python3
"""
Stage 4: replace the spar-cap-face lift (*CLOAD) of a Stage 3 skin-model .inp by
consistent nodal loads from an ANALYTIC PLACEHOLDER pressure distribution on the skin.

NOT an aerodynamic prediction:
  spanwise  : elliptical, lift per unit span q(z) = q0*sqrt(1-(z/1000)^2), 0<=z<=1000
  chordwise : thin-airfoil flat-plate loading, f(xi) = (2/pi)*sqrt((1-xi)/xi), xi from LE
              (centre of pressure at 25% chord; spar axis is at 30% chord)
  split     : antisymmetric. Upper skin: suction (force along +outward normal),
              lower skin: pressure (force along -outward normal); each carries half of the
              net lift, i.e. normal pressure magnitude = p_net/2 on each surface, where
              p_net = q(z) f(xi) / c(z) is the net load per unit PLANFORM area.

Usage (system python, numpy only):
  /usr/bin/python3 scripts/apply_stage4_pressure_load.py <static.inp> <out.inp> <full|lift|grav>
    full : pressure load + existing gravity *DLOAD
    lift : pressure load only (no *DLOAD)     [superposition check]
    grav : gravity only (no *CLOAD)            [superposition check]
"""
import sys
import re
import numpy as np

SPAN = 1000.0
M_SUB = 8                 # sub-triangles per element = M_SUB**2 = 64
CHUNK = 1000
GATE_SCALE = 0.01
GATE_MOMENT = 0.005
GATE_XI = 0.01
GATE_TORQUE = 0.02


def chord(z):
    return 300.0 - 0.15 * z


def G(xi):
    """cumulative of f(xi) = (2/pi) sqrt((1-xi)/xi); G(1) = 1."""
    xi = np.clip(xi, 0.0, 1.0)
    th = np.arcsin(np.sqrt(xi))
    return (2.0 / np.pi) * (th + np.sin(th) * np.cos(th))


def f_pt(xi):
    xi = np.clip(xi, 1e-9, 1.0)
    return (2.0 / np.pi) * np.sqrt((1.0 - xi) / xi)


def Qs(s):
    """integral_0^s sqrt(1-t^2) dt, s clipped to [0,1]."""
    s = np.clip(s, 0.0, 1.0)
    return 0.5 * (s * np.sqrt(1.0 - s * s) + np.arcsin(s))


def is_kw(line):
    return line.startswith("*") and not line.startswith("**")


def block_end(lines, start):
    k = start + 1
    while k < len(lines) and not is_kw(lines[k]):
        k += 1
    return k


def find_kw(lines, pred):
    return [i for i, l in enumerate(lines) if pred(l.strip().upper())]


def data_rows(lines, a, b):
    return [l.split(",") for l in lines[a + 1:b] if l.strip() and not l.startswith("**")]


def subtri(m):
    vtx = []
    for i in range(m):
        for j in range(m - i):
            vtx.append([(i, j), (i + 1, j), (i, j + 1)])
    for i in range(m - 1):
        for j in range(m - 1 - i):
            vtx.append([(i + 1, j), (i, j + 1), (i + 1, j + 1)])
    V = np.array(vtx, float) / m          # (S,3,2)
    return V, V.mean(axis=1)              # centroids (S,2)


def main():
    if len(sys.argv) != 4 or sys.argv[3] not in ("full", "lift", "grav"):
        print(__doc__)
        sys.exit(2)
    src, dst, mode = sys.argv[1], sys.argv[2], sys.argv[3]
    lines = open(src).read().split("\n")

    # ---- load case header
    ilc = [i for i, l in enumerate(lines) if l.startswith("** LOADCASE")]
    assert len(ilc) == 1, "expected exactly one '** LOADCASE' line"
    m = re.search(r"LIFT_N=([0-9.eE+-]+)", lines[ilc[0]])
    assert m, "LIFT_N not found in LOADCASE header"
    L = float(m.group(1))
    print(f"LOADCASE line: {lines[ilc[0]]}")
    print(f"target half-wing lift L = {L:.6f} N")

    # ---- parse nodes / elements
    i_node = find_kw(lines, lambda s: s == "*NODE")
    assert len(i_node) == 1
    rows = data_rows(lines, i_node[0], block_end(lines, i_node[0]))
    nid = np.array([int(r[0]) for r in rows])
    xyz = np.array([[float(v) for v in r[1:4]] for r in rows])
    order = np.argsort(nid)
    sid = nid[order]
    assert len(np.unique(nid)) == len(nid)

    def to_idx(conn):
        k = np.searchsorted(sid, conn)
        assert np.all(sid[k] == conn), "connectivity refers to unknown node"
        return order[k]

    i_s6 = find_kw(lines, lambda s: s.startswith("*ELEMENT") and "TYPE=S6" in s.replace(" ", ""))
    i_c3 = find_kw(lines, lambda s: s.startswith("*ELEMENT") and "TYPE=C3D10" in s.replace(" ", ""))
    assert len(i_s6) == 1 and len(i_c3) == 1
    rs = data_rows(lines, i_s6[0], block_end(lines, i_s6[0]))
    assert all(len(r) == 7 for r in rs), "S6 rows must have 7 fields"
    s6 = to_idx(np.array([[int(v) for v in r[1:7]] for r in rs]))
    rc = data_rows(lines, i_c3[0], block_end(lines, i_c3[0]))
    assert all(len(r) == 11 for r in rc), "C3D10 rows must have 11 fields"
    c3 = to_idx(np.array([[int(v) for v in r[1:11]] for r in rc]))
    solid_mask = np.zeros(len(nid), bool)
    solid_mask[c3.ravel()] = True
    print(f"nodes {len(nid)}, S6 {len(s6)}, C3D10 {len(c3)}, nodes in solid {solid_mask.sum()}")

    # ---- locate CLOAD / DLOAD blocks
    i_cl = find_kw(lines, lambda s: s.startswith("*CLOAD"))
    i_dl = find_kw(lines, lambda s: s.startswith("*DLOAD"))
    assert len(i_cl) == 1 and len(i_dl) == 1, "expected one *CLOAD and one *DLOAD"
    cs, ce = i_cl[0], block_end(lines, i_cl[0])
    ds, de = i_dl[0], block_end(lines, i_dl[0])
    old = data_rows(lines, cs, ce)
    old_F = np.zeros(3)
    for r in old:
        old_F[int(r[1]) - 1] += float(r[2])
    print(f"existing *CLOAD block: {len(old)} lines, sum F = {old_F}")
    print(f"existing *DLOAD block: {de - ds - 1} data lines: {[l.strip() for l in lines[ds + 1:de]]}")

    # ---- consistent nodal loads
    V, Cn = subtri(M_SUB)
    S = V.shape[0]
    u, v = Cn[:, 0], Cn[:, 1]
    Lb = np.stack([1 - u - v, u, v], axis=1)
    Nw = np.empty((S, 6))
    for i in range(3):
        Nw[:, i] = Lb[:, i] * (2 * Lb[:, i] - 1)
    Nw[:, 3] = 4 * Lb[:, 0] * Lb[:, 1]
    Nw[:, 4] = 4 * Lb[:, 1] * Lb[:, 2]
    Nw[:, 5] = 4 * Lb[:, 2] * Lb[:, 0]
    u3, v3 = V[:, :, 0], V[:, :, 1]
    q0 = L / (SPAN * np.pi / 4.0)

    nodeF = np.zeros((len(nid), 3))
    n_up = n_lo = bad_orient = bad_sign = 0
    area_tot = 0.0
    for e0 in range(0, len(s6), CHUNK):
        idx = s6[e0:e0 + CHUNK]
        A, B, C = xyz[idx[:, 0]], xyz[idx[:, 1]], xyz[idx[:, 2]]
        cr = np.cross(B - A, C - A)
        ar = 0.5 * np.linalg.norm(cr, axis=1)
        n = cr / (2.0 * ar)[:, None]
        cen = (A + B + C) / 3.0
        rad = cen.copy()
        rad[:, 0] = 0.0
        rad[:, 1] = 0.0
        rad = cen - rad                                    # centroid - (0,0,zc)
        bad_orient += int(np.sum((n * rad).sum(1) <= 0.0))
        sgn = np.where(cen[:, 1] > 0.0, 1.0, -1.0)
        n_up += int(np.sum(sgn > 0))
        n_lo += int(np.sum(sgn < 0))
        bad_sign += int(np.sum(sgn * n[:, 1] <= 0.0))
        area_tot += ar.sum()

        P = (A[:, None, None, :]
             + u3[None, :, :, None] * (B - A)[:, None, None, :]
             + v3[None, :, :, None] * (C - A)[:, None, None, :])      # (e,S,3,3)
        X, Z = P[..., 0], P[..., 2]
        cz = chord(Z)
        xi = np.clip((X + 0.3 * cz) / cz, 0.0, 1.0)
        xmin, xmax = xi.min(2), xi.max(2)
        d = xmax - xmin
        dsafe = np.where(d > 1e-9, d, 1.0)
        fmean = np.where(d > 1e-9, (G(xmax) - G(xmin)) / dsafe, f_pt(0.5 * (xmin + xmax)))
        za, zb = Z.min(2), Z.max(2)
        dz = zb - za
        dzs = np.where(dz > 1e-9, dz, 1.0)
        zmid = 0.5 * (za + zb)
        q_pt = q0 * np.sqrt(np.clip(1.0 - (zmid / SPAN) ** 2, 0.0, 1.0)) * (zmid >= 0.0)
        qm = np.where(dz > 1e-9, q0 * SPAN * (Qs(zb / SPAN) - Qs(za / SPAN)) / dzs, q_pt)
        pnet = qm * fmean / chord(Z.mean(2))               # N/mm^2 net per planform area
        press = 0.5 * pnet
        Fsub = (sgn[:, None] * press * (ar / M_SUB ** 2)[:, None])[:, :, None] * n[:, None, :]
        Fn = np.einsum("sk,esd->ekd", Nw, Fsub)            # (e,6,3)
        np.add.at(nodeF, idx, Fn)

    print(f"\nskin area {area_tot:.4f} mm^2 | upper elems {n_up}, lower elems {n_lo}")
    print(f"orientation test n.(c-(0,0,zc))<=0 : {bad_orient} elements (expect 0)")
    print(f"elements with sgn*n_y<=0 (force would point down): {bad_sign}")

    Fy_raw = nodeF[:, 1].sum()
    scale = L / Fy_raw
    print(f"raw sum Fy = {Fy_raw:.6f} N -> scale to target = {scale:.6f}")
    ok = True
    if abs(scale - 1.0) > GATE_SCALE:
        print(f"GATE FAIL: scale deviates > {GATE_SCALE*100:.1f}%")
        ok = False
    nodeF *= scale

    Fx, Fy, Fz = nodeF.sum(axis=0)
    print(f"net force after scaling: Fx = {Fx:.6f}  Fy = {Fy:.6f}  Fz = {Fz:.6f} N")

    Mroot = float((nodeF[:, 1] * xyz[:, 2]).sum())
    Mexp = L * SPAN * 4.0 / (3.0 * np.pi)
    err_m = Mroot / Mexp - 1.0
    print(f"root moment sum(Fy*z) = {Mroot:.3f} N*mm | analytic {Mexp:.3f} | err {err_m*100:+.4f}%")
    if abs(err_m) > GATE_MOMENT:
        print("GATE FAIL: root moment")
        ok = False

    zz = np.linspace(0.0, SPAN, 200001)
    qq = q0 * np.sqrt(1.0 - (zz / SPAN) ** 2)
    Texp = float(np.trapz(qq * (-0.05 * chord(zz)), zz))
    Tfe = float((nodeF[:, 1] * xyz[:, 0]).sum())
    err_t = Tfe / Texp - 1.0
    print(f"sum(Fy*x) about spar axis = {Tfe:.3f} N*mm | analytic (CP at 25% chord) {Texp:.3f} | err {err_t*100:+.3f}%")
    print(f"  (for reference sum(-Fx*y) = {float((-nodeF[:,0]*xyz[:,1]).sum()):.3f} N*mm)")
    if abs(err_t) > GATE_TORQUE:
        print("GATE FAIL: torque about spar axis")
        ok = False

    cz_n = chord(xyz[:, 2])
    xi_n = (xyz[:, 0] + 0.3 * cz_n) / cz_n
    print("\nspan bin [mm]      sum Fy [N]   xi_cp (expect 0.250)")
    for b0 in range(0, 1000, 100):
        mk = (xyz[:, 2] >= b0) & (xyz[:, 2] < b0 + 100)
        sf = nodeF[mk, 1].sum()
        xcp = (nodeF[mk, 1] * xi_n[mk]).sum() / sf
        flag = "" if abs(xcp - 0.25) <= GATE_XI else "  <-- GATE FAIL"
        if flag:
            ok = False
        print(f"{b0:5d}-{b0+100:<5d}   {sf:12.4f}   {xcp:8.4f}{flag}")

    loaded = np.any(np.abs(nodeF) > 1e-12, axis=1)
    sh_only = loaded & ~solid_mask
    shared = loaded & solid_mask
    print(f"\nloaded nodes {loaded.sum()} | shell-only {sh_only.sum()} (sum Fy {nodeF[sh_only,1].sum():.4f} N)"
          f" | shared with solid {shared.sum()} (sum Fy {nodeF[shared,1].sum():.4f} N)")
    print(f"max |F| on a node: {np.abs(nodeF).max():.5f} N")

    if not ok:
        print("\nABORT: gate failed, nothing written.")
        sys.exit(1)

    # ---- write output .inp
    new_cl = ["*CLOAD"]
    for k in np.nonzero(loaded)[0]:
        for d in range(3):
            if abs(nodeF[k, d]) > 1e-12:
                new_cl.append(f"{nid[k]}, {d + 1}, {nodeF[k, d]:.8e}")
    ops = []
    if mode in ("full", "lift"):
        ops.append((cs, ce, new_cl))
    else:
        ops.append((cs, ce, []))
    if mode == "lift":
        ops.append((ds, de, []))
    out = lines[:]
    for a, b, new in sorted(ops, key=lambda t: -t[0]):
        out[a:b] = new
    out.insert(ilc[0] + 1,
               f"** STAGE4 MODE={mode}: lift as ANALYTIC PLACEHOLDER skin pressure (elliptical span, "
               f"thin-airfoil chord ~sqrt((1-xi)/xi), antisymmetric upper suction/lower pressure), "
               f"consistent nodal loads, sum Fy normalized to LIFT_N")
    open(dst, "w").write("\n".join(out))

    # ---- re-read check
    chk = open(dst).read().split("\n")
    jc = find_kw(chk, lambda s: s.startswith("*CLOAD"))
    jd = find_kw(chk, lambda s: s.startswith("*DLOAD"))
    if mode in ("full", "lift"):
        assert len(jc) == 1
        rr = data_rows(chk, jc[0], block_end(chk, jc[0]))
        S3 = np.zeros(3)
        for r in rr:
            S3[int(r[1]) - 1] += float(r[2])
        print(f"\nre-read CLOAD: {len(rr)} lines, sum = {S3}")
        assert abs(S3[1] - L) < 1e-6 * L
    else:
        assert len(jc) == 0
    assert len(jd) == (0 if mode == "lift" else 1)
    print(f"written {dst} (mode {mode}); *CLOAD blocks {len(jc)}, *DLOAD blocks {len(jd)}")


if __name__ == "__main__":
    main()
