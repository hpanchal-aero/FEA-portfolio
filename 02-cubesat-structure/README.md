# Project 02 — 1U CubeSat Primary Structure

A staged FEA study of a 1U CubeSat primary structure under
quasi-static launch loading, built with CalculiX, Netgen, and Python.
Full development history lives in
[`development-log.md`](development-log.md); this file is the summary.

## Research Question

How does CubeSat structural configuration affect mass, stiffness,
and stress margin under representative launch loading?

## Approach

Two structural configurations of a 1U CubeSat frame are compared
under an assumed 12 g quasi-static axial launch load:

- **Configuration A** — solid-panel frame: four 8.5 x 8.5 x 100 mm
  corner rails and four 2 mm x 83 mm side panels in a 100 mm cube,
  fused into one solid. Top/bottom panels and bolted joints are not
  modeled (documented simplifications).
- **Configuration B** — same frame, with a through-cut pocket in each
  panel (centered at z = 50, 5 mm perimeter frame retained). Four
  named pocket shapes (rectangle, circle, cross, 2x2 grid) at three
  removal levels (10/20/30% of panel area), plus a random-polygon
  search at the best-performing level.

Material: Al 6061-T6. Total mass assumption 1.33 kg gives 156.6 N
applied axially; base fixed. Constraint: von Mises <= yield/1.5
(~184 MPa). The 276 MPa yield figure and the 1.33 kg CDS mass budget
are commonly quoted values, not yet verified against a primary source.
Modal/vibration response is out of scope (deferred to Project 07).

Configuration A's uniform load is valid only because it is prismatic
(every member carries equal stress under equal strain). Configuration
B instead ties every top-face node's displacement to a reference node
carrying the total load, so it splits between rails and panels by
actual stiffness rather than an assumed equal stress — verified
against Configuration A's own results before use.

No single closed-form reference exists for either full assembly, so
each is verified as an internal-consistency study: isolated-member
checks against uniaxial theory, force equilibrium, mesh convergence,
and (for Configuration A) a roller-base diagnostic. Full details,
including two real modeling bugs found and fixed, are in the
development log.

## Findings — Configuration A

Isolated panel and rail checks (against uniaxial theory) matched to
within 0.001–0.3% on stress and 0.3–1.3% on displacement, with the
Saint-Venant boundary peak correctly located and explained rather than
reported as a design value.

The full assembly closes equilibrium to <0.0001% at both mesh levels
tested. Mean stress and displacement are stable across refinement.
The peak stress (at the fixed-base corners) is **not a converged
design value** — a roller-base diagnostic showed it comes entirely
from the base's lateral restraint, not from geometry or mesh. Mean
stress sits about 1,100x below the allowable.

## Findings — Configuration B: named shapes

All 12 designs (4 shapes x 3 levels) built, meshed, solved, and
verified; equilibrium closes to <0.0001% throughout.

| Shape | Level | Mass (g) | Stiffness (N/mm) | Net stress (MPa) | Stiffness/mass |
|---|---|---|---|---|---|
| rect | 10% | 239.4 | 557,187 | 0.227 | 2327.6 |
| rect | 20% | 221.5 | 482,042 | 0.256 | 2176.8 |
| rect | 30% | 203.5 | 426,935 | 0.281 | 2097.7 |
| circle | 10% | 239.4 | 546,304 | 0.265 | 2282.2 |
| circle | 20% | 221.5 | 460,685 | 0.329 | 2080.3 |
| circle | 30% | 203.5 | 395,045 | 0.407 | 1941.0 |
| cross | 10% | 239.4 | 317,161 | 0.718 | 1324.9 |
| cross | 20% | 221.5 | 306,379 | 0.619 | 1383.5 |
| cross | 30% | 203.5 | 297,091 | 0.561 | 1459.7 |
| grid | 10% | 239.4 | 555,670 | 0.158 | 2321.3 |
| grid | 20% | 221.5 | 475,237 | 0.141 | 2146.0 |
| grid | 30% | 203.5 | 412,189 | 0.120 | 2025.2 |
| Config A baseline | — | 257.3 | 667,635 | 0.164* | 2594.7 |

\* theoretical, not the same nodal-mean method used for Config B.

**Cross is a fundamentally different design, not a scaled-down
version of the others.** Its cutter always spans the panel's full
interior height regardless of removal %, so even at 10% it creates a
near-total structural throat: stiffness stays at 44-48% of Config A
at every level (versus 60-85% for the other three), and net stress
runs 3.4-4.4x the baseline. Kept deliberately as an "equal area, very
different concentration" comparison point.

**Rect and grid lead on stiffness and stiffness/mass**, tracking each
other closely; circle trails modestly. **Grid's net stress sits below
even the Config A baseline and falls with more removal** — consistent
with its four separated blocks preserving continuous load paths, so
material comes disproportionately from low-stress regions. Rect and
circle's net stress rises with removal, as expected for an ordinary
growing cutout.

Net-section stress is a pragmatic nodal-mean proxy (not an integrated
force/area value) and needed a third mesh level to converge for
cross/grid — the 3.0->2.0mm step alone showed a misleadingly *growing*
trend that reversed once a third level was added. See the development
log for the full convergence data and the proxy's known limits (it
does not cleanly separate corner-singularity behavior from bulk
behavior for cross).

## Findings — Configuration B: random-shape search

**Why 10%, not 20% or 30%:** a combined score — 1/3 stiffness ratio to
Config A, 1/3 mass-removed fraction, 1/3 stress margin to the 184 MPa
allowable — averages highest at 10% removal across the three
well-behaved named shapes (0.632 vs. 0.615 at 20%, 0.608 at 30%), so
10% was carried into the random search.

**Generation:** 223 candidate polygons (5-12 vertices, random points
connected by convex hull or angular order, exact-area-scaled, 3 mm
minimum feature size, fixed seed for reproducibility) were generated
and validated as clean single-solid geometry — all 223 also passed a
mesh-feasibility check at every one of the three converged mesh sizes.

**Why 20 of 223, not all 223:** solving all 223 at full convergence
was judged disproportionate to the question being asked. The 20 taken
forward were the candidates with the **lowest worst-case equation-
ceiling usage** across the three mesh sizes — a resource/safety
screen, not a performance filter, since no candidate had been solved
yet at that point. In practice this was closer to a tie-breaker than
a meaningful cut: all 223 were comfortably solvable, and the worst-
case usage across all of them spanned only 48.5-52.9% of the ceiling.

**All 20, solved at 3.0 mm** (score computed the same way as the named
shapes; rect's 3.0 mm score, 0.6343, is the bar to beat):

| Rank | Design | Stiffness (N/mm) | Net stress (MPa) | Score | vs. rect |
|---|---|---|---|---|---|
| 1 | attempt16_hull | 562,908\*\* | 0.2210\*\* | 0.6372\*\* | +0.46% |
| 2 | attempt37_hull | 561,457\*\* | 0.2155\*\* | 0.6365\*\* | +0.35% |
| 3 | attempt54_hull | 560,472\*\* | 0.2348\*\* | 0.6360\*\* | +0.27% |
| 4 | attempt15_hull | 559,901\*\* | 0.2240\*\* | 0.6357\*\* | +0.22% |
| 5 | attempt40_hull | 559,867\*\* | 0.2167\*\* | 0.6357\*\* | +0.22% |
| 6 | attempt114_hull | 558,168 | 0.2355 | 0.6348 | +0.08% |
| 7 | attempt94_hull | 555,852 | 0.2459 | 0.6336 | -0.11% |
| 8 | attempt112_hull | 551,094 | 0.2075 | 0.6313 | -0.47% |
| 9 | attempt124_hull | 550,712 | 0.2365 | 0.6311 | -0.50% |
| 10 | attempt134_hull | 547,680 | 0.2398 | 0.6296 | -0.74% |
| 11 | attempt29_hull | 545,874 | 0.2256 | 0.6287 | -0.88% |
| 12 | attempt13_hull | 545,741 | 0.2332 | 0.6286 | -0.90% |
| 13 | attempt108_hull | 525,584 | 0.2391 | 0.6185 | -2.49% |
| 14 | attempt77_hull | 524,776 | 0.2706 | 0.6181 | -2.55% |
| 15 | attempt102_star | 520,417 | 0.1811 | 0.6161 | -2.87% |
| 16 | attempt26_hull | 519,778 | 0.2470 | 0.6156 | -2.95% |
| 17 | attempt56_star | 515,739 | 0.1990 | 0.6137 | -3.25% |
| 18 | attempt48_hull | 506,317 | 0.1946 | 0.6090 | -3.99% |
| 19 | attempt64_star | 506,051 | 0.1917 | 0.6089 | -4.00% |
| 20 | attempt131_star | 472,702 | 0.1675 | 0.5923 | -6.62% |

\*\* Top 5 carried through full 3-level convergence (2.0, 1.8 mm);
values shown are the converged 1.8 mm result. Ranks 6-20 are 3.0 mm
(screening-level) only.

**Pattern:** all 6 designs beating rect, and 12 of the top 12 overall,
are convex-hull-derived (triangle/quad) polygons. The 3 star-order
(more irregular, concave-capable) polygons in the batch rank 15th,
17th, and 20th — never in the top half. This is a systematic result,
not a coincidence of which 20 happened to be picked: smoother, more
convex pocket boundaries outperform irregular ones at this removal
level and load case, consistent with fewer/gentler stress
concentrations.

The winner (`attempt16_hull`) beats rect by a mesh-converged 0.46%.
Stress contours for the top 5 are in `figures/contours/config_b_
random_attempt{16,37,54,15,40}_hull_h1p8_vm_contour.png` — note that
peak von Mises does **not** track this ranking (attempt15_hull has the
highest peak among the 5 despite ranking 4th), since the score is
driven by net-section stress and stiffness, not the single worst
local point.

**Scope of this result:** best of 20 solved candidates out of 223
generated — not a proven global optimum over the full shape space.

## Limitations

- Prismatic vs. non-prismatic matters for this axial-only load case:
  Configuration A's stress/stiffness reduce exactly to a
  mass-versus-area trade (confirmed by the roller-base diagnostic), so
  the load case only discriminates between designs once panels are
  non-prismatic, as in Configuration B.
- Two mesh levels for Config A and Config B rect/circle/the random
  search's ranks 6-20; three for Config B cross/grid and the top 5
  random designs.
- One element spans each 2 mm wall at the coarsest mesh size; local
  stress gradients at junctions are not fully resolved there.
- Fully fused rails and panels: no bolted-joint contact or local joint
  stress, in either configuration.
- Net-section stress is a nodal-mean proxy with a corner-exclusion
  heuristic that does not work equally well for every pocket shape
  (see cross, above).
- Peak von Mises is never reported as a design value — it is
  boundary-condition- or singularity-dominated and non-converged at
  every mesh resolution used.
- Al 6061-T6 yield (276 MPa) and the 1.33 kg CDS mass budget are not
  yet verified against a primary source.

## Conclusion

Configuration A is a verified prismatic baseline. Among the 4 named
Configuration B shapes, rect and grid lead on stiffness and
stiffness-per-mass, circle trails modestly, and cross — despite
removing the same area — is a fundamentally weaker design due to its
throat geometry. All mass reductions are real (7.0/13.9/20.9% of
structure mass at 10/20/30% panel-area removal).

The primary result of the full study: a subsequent, unconstrained
random-polygon search at the best-performing removal level found
designs that beat the best hand-picked shape, and did so
systematically — every strong performer was a smooth, convex pocket
boundary, never an irregular one. Hand-picked, symmetric pocket shapes
are a reasonable starting point but are not guaranteed optimal even
under a simple removal-level and fairness constraint.
