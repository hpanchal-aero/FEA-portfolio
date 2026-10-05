import numpy as np
F = "simulation/stage3_skin_s4_k8_t1_n3p8.inp"
nd = np.loadtxt(F, delimiter=",", skiprows=5, max_rows=150322)
el = np.loadtxt(F, delimiter=",", skiprows=195228, max_rows=34698, dtype=int)
assert nd[0, 0] == 1 and nd[-1, 0] == 150322 and el[0, 0] == 44900
assert np.all(nd[:, 0] == np.arange(1, 150323))
xyz = nd[:, 1:4]
A, B, C = xyz[el[:, 1] - 1], xyz[el[:, 2] - 1], xyz[el[:, 3] - 1]
cr = np.cross(B - A, C - A); ar = 0.5 * np.linalg.norm(cr, axis=1); n = cr / (2 * ar)[:, None]
cen = (A + B + C) / 3
sgn = np.where(cen[:, 1] > 0, 1.0, -1.0)
bad = sgn * n[:, 1] <= 0
c = 300 - 0.15 * cen[:, 2]
xi = (cen[:, 0] + 0.3 * c) / c
print("count", bad.sum(), "| total area", ar[bad].sum(), "mm^2")
print("xi range", xi[bad].min(), xi[bad].max(), "| z range", cen[bad, 2].min(), cen[bad, 2].max())
print("xi<0.05:", int((xi[bad] < 0.05).sum()), " xi>0.90:", int((xi[bad] > 0.90).sum()), " other:", int(((xi[bad] >= 0.05) & (xi[bad] <= 0.90)).sum()))
print("signed planform area of these (sgn*n_y*A):", float((sgn * n[:, 1] * ar)[bad].sum()), "mm^2")
print("min sgn*n_y:", (sgn * n[:, 1])[bad].min())
for k in np.nonzero(bad)[0][:8]:
    print(f"  elem {el[k,0]} z={cen[k,2]:.1f} xi={xi[k]:.4f} y={cen[k,1]:.3f} n={np.round(n[k],3)} A={ar[k]:.3f}")
