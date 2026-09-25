# Project 02 — 1U CubeSat Primary Structure

A staged FEA study of a 1U CubeSat primary structure under
quasi-static launch loading, built with CalculiX, Netgen, and Python.
Full development history lives in
[`development-log.md`](development-log.md); this file is the summary.

## Research Question

How does CubeSat structural configuration affect mass, stiffness,
and stress margin under representative launch loading?

## Methodology

**Status: Configuration A and Configuration B both complete and
verified within the limits stated below.**

Two structural configurations of a 1U CubeSat frame are compared
under an assumed quasi-static axial launch load:

- **Configuration A** — solid-panel frame: four 8.5 x 8.5 x 100 mm
  corner rails and four 2 mm x 83 mm side panels inside a
  100 x 100 x 100 mm envelope, fused into a single solid. Top and
  bottom panels are omitted, and bolted joints are not modeled (both
  documented simplifications).
- **Configuration B** — same rail/panel envelope, with a through-cut
  pocket in each of the four panels, centered at z=50 with a 5 mm
  perimeter frame retained. Four pocket shapes (rectangle, circle,
  cross, 2x2 grid) at three removal levels (10%, 20%, 30% of panel
  area), with equal removed area across shapes at a given level.

Material: Al 6061-T6 (E = 68,900 MPa, nu = 0.33, rho = 2700 kg/m^3).
Load: 12 g quasi-static axial, a preliminary design assumption and not
a sourced launch-vehicle requirement. A 1.33 kg total mass assumption
gives 156.6 N (156.5779 N as applied). Base face (z = 0) fixed.
Constraint: von Mises <= yield / 1.5 (~184 MPa). Mass is checked
against the 1U CDS budget (1.33 kg) as a design target; that figure is
not verified against the current CDS text, and the yield value
(276 MPa) is a commonly quoted typical value with no source cited yet.

**Load application differs by configuration.** Configuration A
(prismatic) uses uniform 0.1643 MPa pressure on the top face — valid
because every member carries equal stress under equal strain.
Configuration B (non-prismatic, once pocketed) instead ties every
top-face node's axial displacement to a reference node carrying the
total 156.5779 N force (`*EQUATION` + `*CLOAD`), so the load splits
between rails and panels by actual stiffness rather than assuming
equal stress. This mechanism was verified against Configuration A's
own pressure-loaded results (both fixed-base and roller-base) before
use; see development log.

Modal and random-vibration response are out of scope and deferred to
Project 07 (CubeSat modal analysis).

**Toolchain**: geometry built with the Gmsh OCC kernel (geometry only,
no Gmsh meshing), meshed with Netgen (straight-sided C3D10), solved
with CalculiX 2.23, post-processed with Python (NumPy, PyVista,
Matplotlib).

**Verification strategy**: no single closed-form reference exists for
the full assembly. Individual members (an isolated panel and an
isolated rail) were verified against uniaxial (axial-bar) theory, with
closed-form buckling checks. Each assembly (A and B) is an
internal-consistency study: force equilibrium, mesh convergence
(2 levels for A and for Config B's rect/circle; 3 levels for Config
B's cross/grid, where the trend did not stabilize at 2 levels), and,
for Config A, a roller-base diagnostic. Neither is claimed as
closed-form-verified.

## Findings — Configuration A

Isolated members (1.5 mm mesh, fixed base, uniform top pressure):

| | Panel (2a) | Rail (2b) |
|---|---|---|
| Base reaction vs applied | +0.0000% | +0.0000% |
| Top displacement magnitude vs theory | 1.29% below | 0.33% below |
| Mean von Mises, 25 < z < 75 mm | 0.1638 MPa (-0.29%) | 0.1643 MPa (-0.001%) |
| Peak von Mises | 0.345 MPa, at base | 0.259 MPa, at base corners |
| Applied / critical load (closed form) | 0.0011 | 0.0016 |

Configuration A assembly (fixed base, C3D10):

| | 3.0 mm mesh | 2.0 mm mesh |
|---|---|---|
| Base reaction vs 156.578 N | +0.0000% | -0.0000% |
| Top displacement magnitude vs theory | 1.71% below | 1.65% below |
| Mean von Mises, panel region | 0.16424 MPa | 0.16424 MPa |
| Mean von Mises, rail region | 0.16047 MPa | 0.16061 MPa |
| Peak von Mises (base corners) | 0.368 MPa | 0.422 MPa |

- Mean stress and displacement are **stable across the two mesh
  levels**; two levels do not give a convergence order.
- The peak von Mises stress sits at the fixed-base corners, rose 15%
  on refinement, and is **not converged and sensitive to the base
  boundary condition**. A 3.0 mm diagnostic with a roller base gave a
  uniform 0.1643 MPa field and displacement matching theory to
  0.0002%, showing the peak and the ~1.7% displacement offset come
  from the lateral restraint of the fixed base. The peak is not
  reported as a design stress. Mean stress (0.164 MPa) is about 1,100
  times below the ~184 MPa allowable.

## Findings — Configuration B

All 12 designs (4 shapes x 3 levels) built, meshed, solved and
extracted; equilibrium closes to <1e-4% in every case. Mesh levels
used below are the converged level per shape: 2.0 mm for rect/circle
(stable already at 2.0 mm vs 3.0 mm), 1.8 mm for cross/grid (net-
section stress needed a third level to converge — see development
log).

| Shape | Level | Mass (g) | Stiffness (N/mm) | Net-section \|SZZ\| (MPa) | Stiffness/mass |
|---|---|---|---|---|---|
| rect | 10% | 239.38 | 557,187 | 0.227 | 2327.6 |
| rect | 20% | 221.45 | 482,042 | 0.256 | 2176.8 |
| rect | 30% | 203.53 | 426,935 | 0.281 | 2097.7 |
| circle | 10% | 239.38 | 546,304 | 0.265 | 2282.2 |
| circle | 20% | 221.45 | 460,685 | 0.329 | 2080.3 |
| circle | 30% | 203.53 | 395,045 | 0.407 | 1941.0 |
| cross | 10% | 239.38 | 317,161 | 0.718 | 1324.9 |
| cross | 20% | 221.45 | 306,379 | 0.619 | 1383.5 |
| cross | 30% | 203.53 | 297,091 | 0.561 | 1459.7 |
| grid | 10% | 239.38 | 555,670 | 0.158 | 2321.3 |
| grid | 20% | 221.45 | 475,237 | 0.141 | 2146.0 |
| grid | 30% | 203.53 | 412,189 | 0.120 | 2025.2 |
| Config A baseline | — | 257.31 | 667,635 | 0.164 (theoretical) | 2594.7 |

Figures: `figures/config_b_mass_vs_level.png`,
`config_b_stiffness_vs_level.png`, `config_b_netsection_vs_level.png`,
`config_b_stiffness_per_mass.png`.

**Cross is a qualitatively different design, not just a smaller
version of the others.** Its cutter's vertical arm always spans the
full interior panel height (H_INT, independent of removal level), so
even at the 10% level it creates a full-width structural throat with
local area reduced by about 88% at z=50. This shows up as: stiffness
44–48% of Config A's (versus 60–85% for the other three shapes at the
same nominal removal), net-section stress 3.4–4.4x the Config A
baseline (versus 0.7–2.5x for the others), and a peak von Mises stress
around 1.35–1.45 MPa (about 8x Config A's peak) that is dominated by
this throat, not a converged design-stress value. The design point was
kept deliberately, as a legitimate "equal removed area, very different
stress concentration" comparison.

**Rect and grid consistently lead on stiffness and stiffness/mass**,
tracking each other closely at every level; circle trails a modest but
consistent margin behind them. **Grid's net-section stress sits below
even Config A's theoretical baseline and decreases with removal
level** — consistent with the grid pattern's four separated blocks
preserving continuous load-bearing rib material, so removed material
comes disproportionately from low-stress regions. Rect and circle's
net-section stress rises with removal level, as expected for ordinary
stress concentration around a growing central cutout. Why rect and
grid track so closely has not been investigated and is reported as an
open observation, not an explained mechanism.

**Net-section stress is a pragmatic proxy, not an integrated
force/area value**: the mean axial stress (SZZ) over panel-region
corner nodes within z = 50 +/- 1.5 mm, reported both as a full-section
mean and with the highest/lowest decile excluded (a corner-exclusion
attempt). For cross, the excl-decile trim collapses to nearly the same
value as the full mean (gap under 2% in every case) because the entire
throat cross-section is elevated, not just the corners — the proxy
does not cleanly separate corner-singularity behavior from bulk
behavior there. For grid, the metric needed a third mesh level to
converge; the 3.0->2.0mm step alone showed a *growing* percentage
difference with level (10.8% to 20.3%) that looked like a possible
real or metric-flaw trend, but the 2.0->1.8mm step shrank sharply for
every case (to 0.03–2.5%), confirming convergence rather than
divergence. This reversal is recorded in the development log as a
caution against concluding from two mesh levels alone.

### Random-shape exploration (10% level)

Beyond the 4 named shapes, a random-polygon search (223 candidates
generated, 20 solved after a resource-based screen, top 5 converged
across all 3 mesh levels) found designs that beat the best named shape
(rect) on the combined score:

| Design | Stiffness (N/mm) | Net stress (MPa) | Score | vs. rect |
|---|---|---|---|---|
| attempt16_hull (winner) | 562,908 | 0.2210 | 0.6372 | +0.46% |
| rect (best named) | 557,187 | 0.2266 | 0.6343 | -- |

All top performers were convex-hull-derived (triangle/quad) polygons;
no irregular (star-order) polygon reached the top 6. The margin is
modest (0.2-0.5%) but mesh-converged and consistent across all 5
finalists. This screened 20 of 223 generated candidates, not the full
set or a proven global optimum -- see development-log.md for the full
selection funnel, convergence data, and known limitations of the
screening method.
## Limitations

- **Prismatic-vs-non-prismatic distinction matters for the axial-only
  load case.** For Configuration A (prismatic), stress is set by total
  cross-section area and stiffness by EA/L — the roller-base
  diagnostic confirmed this exactly, so for prismatic variants the
  comparison reduces to a mass-versus-area trade. Configuration B's
  pocketed panels are non-prismatic, so this load case retains real
  discriminating content between shapes (as the findings above show);
  lateral or combined quasi-static loading remains a proposed, not
  adopted, extension.
- Two mesh levels for Config A and Config B rect/circle; three for
  Config B cross/grid. The Config A roller-base run is a diagnostic at
  3.0 mm only.
- One element spans each 2 mm wall at Config A's 3.0 mm (checked in
  one panel region); local stress gradients at junctions are not fully
  resolved at the coarser mesh levels used throughout.
- Fully fused rails and panels in both configurations: no bolted-joint
  contact or local joint stress.
- The rail-region mean stress in Config A's assembly is about 2% below
  the panel-region mean; this is not explained.
- Net-section stress is a nodal-mean proxy with a corner-exclusion
  heuristic that does not work equally well for all four pocket
  shapes (see Configuration B findings above).
- Peak von Mises stress is not reported as a design value for either
  configuration; it is boundary-condition- or singularity-dominated
  and non-converged at the mesh resolutions used.

## Conclusion

Configuration A is a verified prismatic baseline. Configuration B
shows that, under this axial-only load case, rect and grid pocket
patterns retain the most stiffness and stiffness-per-mass among the
4 named shapes for a given mass reduction, circle trails them
modestly, and the cross pattern — despite removing the same area —
creates a severe structural throat that is a fundamentally different,
weaker design at every level tested. All mass reductions relative to
Configuration A are real (6.97% / 13.94% / 20.90% of structure mass at
10/20/30% panel-area removal); the corresponding stiffness losses vary
sharply by shape.

A subsequent random-polygon search at the 10% level found designs
(convex, triangle/quad pocket shapes) that beat rect, the best named
shape, by a modest but mesh-converged margin (0.2-0.5% on the combined
score). This is the primary engineering result of the full study:
hand-picked, symmetric/regular pocket shapes are a reasonable
starting point, but are not guaranteed optimal even within a simple
removal-level and fairness constraint — an unconstrained shape search
found better designs, and did so systematically (every top performer
was convex-hull-derived, none were irregular/star-derived), suggesting
"smoother, more convex pocket boundaries" is a real, generalizable
design heuristic for this load case rather than a one-off result.

## Reproducibility

```bash
cd 02-cubesat-structure
# Configuration A (see development log for full command sequence)
/usr/bin/python3 scripts/build_config_a_base_frame.py
/usr/bin/python3 scripts/mesh_config_a.py 3.0
/usr/bin/python3 scripts/run_ccx_watchdog.py simulation/config_a_h3p0 8
/usr/bin/python3 scripts/extract_config_a.py simulation/config_a_h3p0

# Configuration B: geometry (12 designs)
for s in rect circle cross grid; do
  for lvl in 10 20 30; do
    /usr/bin/python3 scripts/build_config_b_assembly.py $s $lvl
  done
done

# Configuration B: mesh + solve + extract (mesh size per shape: see README findings)
/usr/bin/python3 scripts/mesh_config_b.py <shape> <level_pct> <maxh_mm>
/usr/bin/python3 scripts/run_ccx_watchdog.py simulation/config_b_<shape>_<level>pct_h<mesh> 8
/usr/bin/python3 scripts/extract_config_b.py simulation/config_b_<shape>_<level>pct_h<mesh>

# Configuration B: cross-design comparison table and figures
/usr/bin/python3 scripts/compare_config_b.py
```

Netgen scripts run under the system Python (`/usr/bin/python3`); the
solver is CalculiX 2.23 from the `ccx223` conda environment.
