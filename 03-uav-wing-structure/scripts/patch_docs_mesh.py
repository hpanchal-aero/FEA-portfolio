#!/usr/bin/env python3
"""
Guarded documentation patch: skin-mesh sensitivity (Stage 6, section 6).
Numbers are read from results/stage6d_mesh_sensitivity.csv.
Aborts (writing nothing) unless every anchor occurs exactly once.
Backups: /tmp/README.md.bak_docs, /tmp/development-log.md.bak_docs
Run: /usr/bin/python3 scripts/patch_docs_mesh.py
"""
import csv
import shutil

rows = list(csv.DictReader(open("results/stage6d_mesh_sensitivity.csv")))
D = {}
for r in rows:
    D[(r["layout"], float(r["skin_t"]), float(r["skin_mesh"]))] = r


def f(key, col):
    return float(D[key][col])


def pct(x):
    return "%+.2f%%" % x


# ---------- numbers from the CSV ----------
need = [("L100", 1.0, k) for k in (10, 8, 6, 5)] + \
       [(l, 1.0, k) for l in ("L126", "L006") for k in (10, 8, 6)] + \
       [(l, 0.6, k) for l in ("L126", "L006") for k in (10, 8, 6)]
for n in need:
    if n not in D:
        raise SystemExit("ABORT: CSV is missing %s" % (n,))

tip_max = max(abs(float(r["d_tip_pct"])) for r in rows)
p99 = [float(r["d_p99_pct"]) for r in rows if r["d_p99_pct"] not in ("", "None", "n/a")]
b1_10 = [f((l, 1.0, 10), "d_buck_pct") for l in ("L100", "L126", "L006")]
b1_6 = [f((l, 1.0, 6), "d_buck_pct") for l in ("L100", "L126", "L006")]
b06_10 = [f((l, 0.6, 10), "d_buck_pct") for l in ("L126", "L006")]
b06_6 = [f((l, 0.6, 6), "d_buck_pct") for l in ("L126", "L006")]
l100_5 = f(("L100", 1.0, 5), "d_buck_pct")
l100_6 = f(("L100", 1.0, 6), "d_buck_pct")
step65 = l100_5 - l100_6
gap = {k: 100 * (f(("L006", 1.0, k), "buckling") / f(("L126", 1.0, k), "buckling") - 1) for k in (10, 8, 6)}
mean06_6 = sum(b06_6) / 2.0
l151_est = 1.507 * (1 + mean06_6 / 100.0)   # 1.507 = L151 buckling at 0.6 mm, 8 mm skin mesh (README table)
nexp = 2.2
tstar = {}
for lay, t0 in (("L126", 0.605), ("L006", 0.603)):
    d6 = f((lay, 0.6, 6), "d_buck_pct")
    tstar[lay] = t0 * (1 / (1 + d6 / 100.0)) ** (1 / nexp)


def table():
    out = ["| Layout | Skin | 10 mm | 8 mm (baseline) | 6 mm | 5 mm | Max tip change |", "|---|---|---|---|---|---|---|"]
    for lay, skin in sorted({(r["layout"], float(r["skin_t"])) for r in rows}):
        cells = []
        for k in (10, 8, 6, 5):
            cells.append(pct(f((lay, skin, k), "d_buck_pct")) if (lay, skin, k) in D else "not run")
        tm = max(abs(f((lay, skin, k), "d_tip_pct")) for k in (10, 8, 6, 5) if (lay, skin, k) in D)
        out.append("| %s | %g mm | %s | %.4f%% |" % (lay, skin, " | ".join(cells), tm))
    return "\n".join(out)


README_SEC = """### 6. Skin mesh sensitivity (solid mesh 4 mm; skin mesh 10, 8, 6 and 5 mm)

The Stage 3–5 mesh (skin 8 mm) was used for all Stage 6 layouts. Because buckling is a local aft-skin phenomenon, a limited check was made by re-meshing selected cases with a different skin mesh size (solid mesh unchanged) and comparing with the 8 mm baseline of the same case. Change in critical buckling factor relative to 8 mm (solver results):

%(table)s

- **Tip deflection and mass are mesh-insensitive**: the largest tip change over all 16 runs is %(tip).4f%%.
- **Buckling is mesh-sensitive at the 0.5–1.5%% level and the 8 mm baseline is slightly unconservative.** Refining 8 → 6 mm lowers buckling by %(b1a).2f%% to %(b1b).2f%% at 1.0 mm skin and by %(b6a).2f%% to %(b6b).2f%% at 0.6 mm; coarsening to 10 mm raises it by %(c1a).2f%% to %(c1b).2f%% (1.0 mm) and about %(c6).2f%% (0.6 mm). The thinner skin is about twice as sensitive.
- **Convergence is not demonstrated.** For L126 and L006 only three levels (10, 8, 6 mm) could be run, and the change per step grows toward the finer mesh. L100 was also run at 5 mm: 6 → 5 mm changes buckling by only %(step).2f%% (against %(l6).2f%% for 8 → 6 mm), which suggests flattening on that layout, but it is a single 1 mm step on a single layout and is not proof of an asymptote. L126 at 5 mm exceeds the equation-count ceiling, and 4 mm skin meshes are infeasible for all layouts. No Richardson extrapolation or grid-convergence index is claimed.
- **p99 skin stress is not converged**: it changes by %(p99lo)+.2f%% to %(p99hi)+.2f%% across the runs and keeps drifting with refinement, so it is quoted to about ±4%%. Stress is not binding, so no conclusion depends on it.
- **What survives:** layout rankings (L006 exceeds L126 by %(g10).2f%%, %(g8).2f%% and %(g6).2f%% at 10, 8 and 6 mm, i.e. a constant gap), all tip comparisons, and the lightest-at-requirement ordering. Applying the 6 mm shift at 0.6 mm skin with a local exponent of 2.2 moves the skin gauge for buckling 1.5 from 0.605 to about %(tL126).3f mm (L126) and from 0.603 to about %(tL006).3f mm (L006), a shift of about 1%% (estimate from the shift, not a re-solve).
- **What does not survive as stated:** the pass of L151 at 0.6 mm (1.507, 0.5%% above the requirement at the 8 mm mesh). The mesh shift at 0.6 mm is about %(m06).1f%%, which would put it near %(l151).2f (estimate; L151 was not re-solved), so that pass is within the mesh uncertainty and should not be treated as established.
- **Earlier evidence is consistent in size but not in sign:** the Stage 3 skin-mesh check (8 → 6 mm, development log) gave buckling −0.5%% at t = 1.0 mm and +2.2%% at t = 0.5 mm, and the solid mesh change 2.5 → 4 mm gave a tip change of −0.02%%. The sign of the thin-skin effect therefore depends on the case, which is a further reason to treat about ±1.5%% as the uncertainty on buckling at 0.6 mm skin rather than as a correction.

![Skin mesh sensitivity](figures/stage6d_mesh_sensitivity.png)

""" % dict(table=table(), tip=tip_max, b1a=min(b1_6), b1b=max(b1_6), b6a=min(b06_6), b6b=max(b06_6),
           c1a=min(b1_10), c1b=max(b1_10), c6=sum(b06_10) / 2.0, step=step65, l6=l100_6,
           p99lo=min(p99), p99hi=max(p99), g10=gap[10], g8=gap[8], g6=gap[6],
           tL126=tstar["L126"], tL006=tstar["L006"], m06=mean06_6, l151=l151_est)

DEVLOG_SEC = """### Mesh-sensitivity check (added after the Stage 6 commit)

**Decision (user, "decision A"):** check the top layouts against a finer skin mesh. Solid mesh fixed at 4 mm; skin mesh 10, 8 (baseline), 6 and 5 mm; L100, L126 and L006 at 1.0 mm skin and L126 and L006 at 0.6 mm skin. L151 was not run (74.8%% of the equation ceiling already at 8 mm). L126 at 5 mm would exceed the equation gate and was not attempted.

**Pipeline change:** `scripts/patch_driver_mesh.py` (17 guarded edits) added `--mesh SOLID,SKIN` and `--stage conv` to `scripts/run_layout_study.py`; records are keyed with a `_s4_k<skin>` suffix and reuse the baseline imprint. My first conv guard was over-strict and was relaxed; the relaxed guard requires the base record complete and the step, json and vol files present. The default (no `--mesh`) path was regression-tested as unchanged earlier.

**Results** (`results/stage6d_mesh_sensitivity.csv`, `scripts/analyse_mesh_sensitivity.py`; buckling change against the 8 mm baseline of the same case):

%(table)s

**Predictions (written before the runs) and outcome:**

| Quantity | Predicted | Observed | Result |
|---|---|---|---|
| 0.6 mm skin, 8 → 6 mm buckling | +0.5%% to +3%% | %(m06).2f%% (mean of L126, L006) | wrong sign |
| 0.6 mm skin, 10 mm buckling | −1%% to −4%% | +%(c6).2f%% | wrong sign |
| L100 1.0 mm, 5 mm buckling | −1.0%% to −0.8%% (continued decline) | %(l5).3f%% | missed; flatter |
| L100 1.0 mm, 5 mm tip | within about ±0.01%% | +0.0100%% | hit (at the edge) |
| L100 1.0 mm, 5 mm p99 | −3%% to −4%% | −3.55%% | hit |

Only predictions whose numbers are recorded in this log are tabulated; the sign error at 0.6 mm is the important miss, since it means the skin-mesh effect on buckling cannot be predicted from the 1.0 mm behaviour.

**Retractions:** (1) my first reading of the L126/L006 sequences over 10/8/6 mm as "not converging" was revised after L100 at 5 mm; (2) the interim statement that L100 buckling is "converged to about 0.04%%" overstated a single 1 mm step and is not made in the README, which says "suggests flattening".

**Netgen:** two more intermittent failures (L126 at 10 mm, a 150 s timeout; L006 at 0.6 mm and 10 mm, a failed first attempt); both succeeded on retry, cause still unknown.

**Findings:** tip converged to %(tip).4f%% maximum change; buckling baseline about 0.5–1.5%% unconservative; rankings and tip conclusions robust; the L151 pass at 0.6 mm (+0.5%% at the 8 mm mesh) is not robust; p99 stress drifts by %(p99lo)+.2f%% to %(p99hi)+.2f%%. Limits: asymptotic convergence not shown (L126 and L006 stop at 6 mm), no GCI, thin-skin sign differs from the Stage 3 L000 check.

**Reproduce:** `scripts/run_layout_study.py --stage conv --mesh 4,10 --layouts L100,L126,L006` (likewise `4,6`; `4,5 --layouts L100`; add `--skin-t 0.6` for the thin skin), then `scripts/analyse_mesh_sensitivity.py`. Figure: `figures/stage6d_mesh_sensitivity.png`.

""" % dict(table=table(), m06=mean06_6, c6=sum(b06_10) / 2.0, l5=l100_5, tip=tip_max, p99lo=min(p99), p99hi=max(p99))

# ---------- guarded edits ----------
R = open("README.md").read()
V = open("development-log.md").read()


def edit(text, old, new, label):
    n = text.count(old)
    if n != 1:
        raise SystemExit("ABORT: anchor '%s' occurs %d times (need 1)" % (label, n))
    return text.replace(old, new)


R = edit(R, "no finer-mesh convergence check on\nthe Stage 6 layouts yet.",
         "only a limited skin-mesh sensitivity check on\nthree layouts (section 6 of the Stage 6 findings; not a convergence study).", "README status")
R = edit(R, "At 0.6 mm only the 10-rib layout meets the buckling requirement of 1.5, and then by 0.5%;",
         "At 0.6 mm and the 8 mm skin mesh only the 10-rib layout meets the buckling requirement of 1.5, and then by 0.5% (smaller than the skin-mesh effect, see section 6, so this pass is not robust);", "README L151")
R = edit(R, "### Engineering interpretation\n", README_SEC + "### Engineering interpretation\n", "README eng interp")
R = edit(R, "- Only verification was performed (baseline and repeat reproduction, fit check). There is no validation against experiment.",
         "- Skin-mesh sensitivity was checked on three layouts only (section 6): buckling at the 8 mm baseline is about 0.5–1.5% unconservative, asymptotic convergence is not shown, and the p99 stress is not converged. The solid mesh was checked only in Stage 3.\n- Only verification was performed (baseline and repeat reproduction, fit check, limited mesh sensitivity). There is no validation against experiment.", "README limitations")
R = edit(R, "`scripts/analyse_skin_sensitivity.py`, `scripts/modal_span_distribution.py`",
         "`scripts/analyse_skin_sensitivity.py`, `scripts/analyse_mesh_sensitivity.py` (after `--stage conv --mesh SOLID,SKIN`), `scripts/modal_span_distribution.py`", "README reproduce")

old_bullet = "- No Stage 6 mesh-convergence check has been run. A finer-mesh check on the top layouts was agreed (decision A) but is not done; the baseline mesh is the Stage 3–5 mesh (solid 4 mm, skin 8 mm, Netgen)."
V = edit(V, old_bullet,
         "- Mesh: only a limited skin-mesh sensitivity check was run (subsection \"Mesh-sensitivity check\" below): three layouts, skin mesh 10/8/6 mm (L100 also 5 mm); asymptotic convergence is not shown, the baseline mesh (solid 4 mm, skin 8 mm, Netgen) is about 0.5–1.5% unconservative for buckling, and p99 stress is not converged.", "devlog bullet")
V = edit(V, "finer-mesh check of the top layouts not run; committed locally (9684c45), nothing pushed",
         "committed locally (9684c45) before the limited skin-mesh sensitivity check, which is documented afterwards in the Stage 6 entry; nothing pushed", "devlog status")
V = edit(V, "\n## Status Log\n", "\n" + DEVLOG_SEC + "## Status Log\n", "devlog Status Log heading")

shutil.copy("README.md", "/tmp/README.md.bak_docs")
shutil.copy("development-log.md", "/tmp/development-log.md.bak_docs")
open("README.md", "w").write(R)
open("development-log.md", "w").write(V)
print("PATCH_OK  tip_max=%.4f  l100 step 6->5=%+.3f  L151 estimate=%.3f  t*: L126 %.3f L006 %.3f" % (tip_max, step65, l151_est, tstar["L126"], tstar["L006"]))
