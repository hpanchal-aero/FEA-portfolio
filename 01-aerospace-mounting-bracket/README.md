# Project 01 — Ribbed L-Bracket: Mass-Minimized Aerospace Mounting Bracket

**Status:** Verification plate stage complete. Progressive geometry —
holes, shoulder fillet, and right-angle frame (with and without
gusset) stages complete (all verified, with documented open items).
Full L-bracket not yet started.

## Research Question

How can the mass of a ribbed L-bracket mounting interface be minimized
while satisfying stress, displacement, and factor-of-safety constraints
under a representative static equipment load?

## Engineering Motivation

Mounting brackets are ubiquitous in aerospace structures and are a
common target for structural mass optimization, since bracket mass
is rarely mission-critical functionally but directly impacts overall
vehicle mass budget.

## Project Sequencing

The full ribbed L-bracket has no closed-form analytical solution, so
this project follows a staged verification approach:

1. **Verification plate** (complete) — plain flat cantilever plate,
   analytically verifiable via Euler-Bernoulli beam theory. Purpose:
   prove the OpenSCAD → Gmsh → CalculiX → Python pipeline is
   implemented correctly before trusting it on a geometry with no
   independent check.
2. Progressive geometry introduction — **holes (complete)** →
   **shoulder fillet (complete)** → **right-angle frame, with and
   without gusset (complete)**
3. Full L-bracket baseline model and mesh convergence, combining the
   fillet, holes, and gusset onto the bent frame geometry
4. Parametric mass-minimization study (the actual research question)

**Note on sequencing:** the original plan scoped "gusset" as a single
stage on a straight cantilever. Mid-project, the decision was made to
study the gusset on an angled right-angle frame instead (a closer
analogue to the final bracket), which necessarily introduced the bend
one stage earlier than originally planned. This is recorded here as
an explicit, discussed change to the plan, not a silent scope shift.

## Verification Plate — Geometry & Material

Plain rectangular plate, 60 × 40 × 4 mm:
- x — cantilever/load-path direction (length, L = 60 mm)
- y — width direction (B = 40 mm)
- z — thickness direction (H = 4 mm)

Material: Al 7075-T6, E = 71,700 MPa, ν = 0.33, ρ = 2810 kg/m³

## Loading

- Equipment mass: 2 kg; design acceleration: 12g quasi-static
  (**preliminary assumed value, not a sourced launch/vibration
  requirement**)
- Resultant static load: F = 235.44 N
- Applied as a statically-equivalent transverse force at the free-edge
  centroid, via a reference node coupled to the tip face through a
  `*DISTRIBUTING` constraint (area-weighted, not rigid — avoids
  artificially stiffening the tip cross-section)
- Root (x = 0): fully fixed (encastre)

## Analytical Reference (Euler-Bernoulli)

Primary target, per locked assumptions (L/h = 15, reasonable but
approximate for slender-beam theory):

- I = 2.1333 × 10⁻¹⁰ m⁴
- σₓₓ(x=0) = 132.4 MPa
- δ_tip = 1.108 mm

A wide-plate stiffness-corrected deflection estimate (0.988 mm,
via E → E/(1−ν²)) is reported as a secondary reference only — it does
not represent an exact target, since the true finite-width plate
response is not fully captured by a simple modulus substitution.

## Mesh Methodology

Two independent approaches were used and cross-checked:

1. **Unstructured tetrahedra** (C3D10) from an OpenSCAD STL export via
   Gmsh's STL→ClassifySurfaces→CreateGeometry pipeline. Used for the
   initial hand-verified baseline case (1.3 mm element size).
2. **Structured transfinite hexahedra** (C3D20), built directly from
   an exact Gmsh OpenCASCADE box (no STL faceting), used for the full
   5-level mesh convergence study.

**The structured hex approach is scoped to this simple rectangular
plate only.** It relies on the geometry being a trivial box; it is
*not* adopted as the general meshing strategy for this project. The
final ribbed/filleted L-bracket will require a separate mesh
methodology decision (very likely unstructured tetrahedra, since
transfinite hex meshing does not extend to arbitrary fillet/hole
geometry).

The unstructured-tet route required substantial debugging (STL
surface-splitting artifacts, CPS6 skin-element conflicts, node-type
vs. facial `*DISTRIBUTING` surfaces, and Jacobian-degenerate elements
near sharp corners) — all resolved and documented in the development
history. The structured-hex route avoided these issues by
construction, which is the reason it was adopted for the convergence
study specifically.

## Mesh Convergence Study

5 refinement levels (structured hex), element counts from 551 to
344,071:

| Level | Elements | Tip Uz (mm) | Root-band max VM (MPa) |
|---|---|---|---|
| 2.0mm | 551 | −1.0440 | 127.70 |
| 1.3mm | 2,700 | −1.0487 | 120.38 |
| 0.8mm | 14,504 | −1.0500 | 121.52 |
| 0.5mm | 65,807 | −1.0505 | 122.37 |
| 0.3mm | 344,071 | −1.0507 | 122.40 |

Convergence criterion (<2% change between final two levels): **met**
— 0.0245% (displacement), 0.0265% (stress) between 0.5mm and 0.3mm.
The 0.3mm result is treated as the converged reference.

**Non-monotonic stress behavior (2.0→1.3→0.8mm) investigated and
explained:** the windowed-max metric always samples the node nearest
the window's near-root edge; since discrete node x-positions differ
slightly between mesh densities, each level was measuring a marginally
different physical location, not indicating non-convergence of the
underlying field. Confirmed by direct node-location tracking across
all 5 levels — see development history.

## Verification: σₓₓ vs. Euler-Bernoulli Theory (Primary Comparison)

Beam theory directly predicts σₓₓ = M(x)c/I, so this — not a windowed
von Mises maximum — is the correct primary analytical comparison:

| x (mm) | Theory (MPa) | FEA (MPa) | % diff |
|---|---|---|---|
| 3 | 125.81 | 138.50 | +10.08% |
| 6 | 119.19 | 125.77 | +5.52% |
| 9 | 112.57 | 115.10 | +2.25% |
| 12 | 105.95 | 106.97 | +0.96% |
| 20 | 88.29 | 87.09 | −1.36% |
| 30 | 66.22 | 65.62 | −0.91% |
| 40 | 44.15 | 44.36 | +0.50% |

**Mid-span agreement (x = 12–40 mm) is within ±1.4%** — this is
excellent, near-exact confirmation of Euler-Bernoulli theory once away
from boundary-affected regions.

**Root-region deviation (x = 3–9 mm) is a Saint-Venant boundary-layer
effect**, not a modeling error: the rigid `*BOUNDARY` encastre
constraint at x=0 forces the full cross-section to zero displacement,
producing a local 3D stress concentration that decays with distance
from the boundary — consistent with classical Saint-Venant behavior,
fading within roughly 2–3 plate thicknesses (H=4mm; effect largely
gone by x≈12mm ≈ 3H).

This is corroborated directly by the transverse stress ratio
σᵧᵧ/σₓₓ, which decays smoothly from 0.323 (near ν=0.33) at x=3mm to
0.028 at x=40mm — confirming a biaxial stress state near the root
(cross-section restrained against anticlastic curvature by the rigid
BC) that relaxes toward the uniaxial state beam theory assumes, as
distance from the root increases.

**Tip displacement (−1.051 mm converged vs. −1.108 mm target, −5.2%)**
reflects the same 3D effects: real cross-sectional warping was
confirmed directly (0.4% std-dev in Uz across the tip face — beam
theory assumes perfectly plane sections) plus the accumulated
influence of the root-region stiffness effects on overall compliance.

## Verification vs. Validation

This stage constitutes **verification only** — confirming the FE model
solves the governing equations correctly and the solution is
mesh-independent. **No validation was performed or claimed** — no
physical test article or published dataset exists for this specific
geometry, and this is an explicit, permanent limitation of this
project stage.

## Cross-Discretization Corroboration

The unstructured-tet (1.3mm) and structured-hex (0.3mm converged)
results agree closely:

| Mesh type | Tip Uz | Root-band max VM |
|---|---|---|
| Tet (1.3mm) | −1.049 mm | 121.1 MPa |
| Hex (0.3mm, converged) | −1.051 mm | 122.4 MPa |

This demonstrates the FEA implementation itself (material properties,
boundary conditions, load application, coupling constraint, solver
configuration) is correct and consistent, independent of discretization
choice. It does **not** independently validate against physical
reality, and does not by itself explain the beam-theory departure —
both meshes solve the same 3D continuum model, so their mutual
agreement is expected regardless of that model's relationship to the
1D idealization.

## Figures (Verification Plate)

1. `mesh_visualization.png` — structured hex mesh (1.3mm shown for
   visual clarity; 0.3mm mesh used for all reported results)
2. `displacement_contour.png` — deformed shape (20× exaggerated),
   colored by displacement magnitude
3. `sigma_xx_verification.png` — **primary verification plot**: FEA
   σₓₓ vs. Euler-Bernoulli theory across the span
4. `sigma_xx_contour.png` — full-field σₓₓ spatial distribution
5. `biaxial_ratio_decay.png` — σᵧᵧ/σₓₓ ratio vs. x, showing Saint-Venant
   decay toward the free tip
6. `von_mises_contour.png`, `von_mises_root_zoom.png` — von Mises
   stress (combined biaxial effect); root-zoom shows the x=0
   singularity excluded from all quantitative comparisons
7. `convergence_displacement.png`, `convergence_stress.png` — mesh
   convergence trends

## Engineering Interpretation (Verification Plate)

The 3D solid FEA model reproduces Euler-Bernoulli beam theory to
within ~1.4% across the mid-span, and the root-region deviation is
fully explained by a known, physically expected 3D boundary effect
(Saint-Venant decay from the rigid encastre constraint), not by
modeling error. Two independent mesh discretizations agree closely.
The implementation is considered verified and the pipeline trustworthy
for the next stage.

## Limitations (Verification Plate)

- No physical validation performed (no test data available for this
  specific geometry)
- The 12g design acceleration is a preliminary assumed value, not a
  sourced requirement
- The structured hex meshing methodology is not general-purpose and
  will not extend to the filleted/holed L-bracket geometry
- L/h = 15 slenderness ratio is "reasonable but approximate" for
  strict Euler-Bernoulli applicability

---

# Progressive Geometry Stage 2 — Circular Hole

## Geometry

Same 60 × 40 × 4 mm plate as the verification stage, with a single
circular through-hole added:

- Hole diameter: d = 5 mm (matches the eventual bolt-hole spec)
- Hole location: x = 30 mm, y = 20 mm (mid-span, mid-width) — chosen
  to sit well clear of the root Saint-Venant region characterized in
  Stage 1, in a zone already confirmed to match beam theory closely
- d/B = 5/40 = 0.125, d/t = 5/4 = 1.25

**Why this stage matters:** unlike the fillet and gusset stages ahead,
a hole has an established analytical/numerical reference point (a
stress concentration factor), giving one more opportunity to check the
pipeline against something external before that check disappears at
the full-bracket stage.

Geometry built via Gmsh OpenCASCADE (`Box()` cut by `Cylinder()` via
`BooleanDifference`) — an exact circular boundary, no STL faceting.
Meshed with unstructured tetrahedra (C3D10), since structured hex
cannot wrap a circular cutout. Local refinement at the hole via a
Gmsh `Distance`+`Threshold` field pair, sized so `DistMax` is always
guaranteed greater than the hole radius (an earlier field-definition
bug that used an absolute-distance parameter incorrectly caused one
non-physical 2.57M-element mesh — caught and fixed before any solve).

## Mesh Convergence Study

Two refinement levels, both solved with the direct (spooles) solver
— the iterative solver was confirmed to stall on this mesh's local
size-gradient before either level below was attempted:

| Level | Hole element size | Elements | Peak σₓₓ (top surface) | Kt | % change |
|---|---|---|---|---|---|
| 1 | 0.25 mm | 62,306 | 134.966 MPa | 2.038 | — |
| 2 | 0.20 mm | 117,205 | 135.350 MPa | 2.044 | 0.28% |

Convergence criterion (<2%): **met**, decisively — 0.28% change in
peak stress between levels, with the peak location unchanged (same
z=4.00mm top-surface point, x≈29.8–29.95mm, y≈17.5/22.5mm — the
classical Kirsch-type location at 90° from the load axis). No third
level was needed.

**Converged result: Kt ≈ 2.044, peak σₓₓ ≈ 135.35 MPa**, against a
gross-section nominal bending stress (unperforated beam theory at
x=30mm) of ≈ 66.2 MPa. This is treated as a trustworthy,
mesh-independent FEA result, established independently of any
literature comparison.

## Literature Investigation

The naive first-pass estimate for this geometry — applying the
classical finite-width Kirsch/Howland uniaxial-tension chart
layer-by-layer through the thickness, using the local beam-bending
stress at each depth as that layer's "remote tension" — predicted
Kt ≈ 2.68 (peak ≈ 176 MPa), a −23.8% overprediction relative to the
converged FEA result. A literature investigation was carried out to
determine whether this gap is physically explainable, in three rounds:

**Round 1 — uniform-tension 3D-thickness literature** (Vaz et al. 2013;
Sternberg & Sadowsky; Folias & Wang). These are the standard citable
sources for "3D thickness effects on hole SCF," but all study a hole
under **uniform remote tension through the thickness** — a
fundamentally different loading case from ours, where the nominal
stress is a **linear bending gradient** (zero at mid-plane, maximum at
the surfaces). Applying Vaz et al.'s own surface-vs-midplane transition
criterion to our B/r=1.6 predicted the wrong peak location (mid-plane)
for the wrong reason — not because our result violates their finding,
but because their finding answers a different question than ours asks.
**Conclusion: not applicable, ruled out correctly rather than forced.**

**Round 2 — through-thickness bending literature** (Yang, Kim, Beom &
Cho, 2010, *Int. J. Mechanical Sciences* 52, 836–846; Peterson's
"Simple Transverse Bending" hole-in-plate chart). This is the correct
loading type — a linear stress gradient through the thickness, matching
our case exactly:

- Yang et al. (2010) study exactly this configuration (pure bending,
  m=−1) and report a transition thickness of t′≈3r, below which the
  peak SCF remains **on the free surface** and above which it migrates
  into the interior. Our t/r = 1.6 < 3 places us in the "peak-stays-
  on-surface" regime — **which matches our FEA observation** (peak at
  z=0/z=4, not z=2). This is a genuine, correctly-transferred
  qualitative confirmation. No tabulated Kt-vs-thickness data for pure
  bending could be extracted from the available text/figures, however
  — the underlying curve (their Fig. 16a) is described only in words.
- Peterson's "Simple Transverse Bending" chart, evaluated with zero
  extrapolation at our exact d/D=0.125 and d/t=1.25, gives Kt≈2.39 and
  a predicted peak of ≈181 MPa when matched to its own net-section
  moment-based nominal stress definition — a **larger** discrepancy
  from our converged 135.35 MPa than the original naive estimate.

**Round 3 — hole in a finite cantilever with a longitudinal moment
gradient.** Searched specifically for this combination (through-
thickness bending gradient + circular hole + finite beam length near a
free end). Nothing found addresses it. The one directly relevant hit
(an ASEE cantilever-beam-with-hole study) turned out to be the
classical **in-plane** bending case (hole diameter referenced against
beam height, not thickness) — a different configuration already
correctly ruled out earlier in this investigation. No source was found
studying a through-thickness bending gradient superimposed on a
longitudinal moment gradient near a beam's free end.

### Established vs. Not Established

**Established:**
- Kt ≈ 2.044 is mesh-converged (independently, via two refinement
  levels, 0.28% change).
- The stress peak occurs at the free surface, not the mid-plane.
- Yang et al. (2010) qualitatively supports surface-peak behavior for
  a plate this thin relative to the hole under pure through-thickness
  bending (t/r=1.6 < their transition thickness of ~3r).
- Every classical idealized correlation found for the correct loading
  type — evaluated at our exact geometry, no extrapolation —
  *overpredicts* our converged result (Kt=2.39–2.68 predicted vs.
  Kt=2.044 observed).

**Not established:**
- No universal quantitative correction exists in the literature found
  that converts any classical Kt to our observed 2.044 for this exact
  geometry.
- The hypothesis that the cantilever's local along-beam moment
  gradient (rather than pure 3D surface/thickness effects) accounts
  for some or all of the residual gap is a **plausible engineering
  interpretation, not a demonstrated cause**. No source addressing
  this specific combination was found in three rounds of targeted
  search.

This is documented as the final Stage 2 literature conclusion: the
FEA result is trusted on its own mesh-convergence merits; the
mechanism is qualitatively, but not quantitatively, supported by
published work; the residual gap to any classical prediction is real
and, as far as this search found, unresolved in the literature.

## Verification vs. Validation (Stage 2)

As with Stage 1, this stage is **verification only**. The FEA result
is confirmed mesh-independent and physically sensible (correct peak
location per the one directly-relevant reference found), but is not
validated against a physical test article, and is not fully
reconciled against any closed-form or tabulated literature prediction.

## Limitations (Stage 2)

- No closed-form or tabulated literature source was found that
  quantitatively reproduces the converged Kt ≈ 2.044 for this exact
  geometry (d/t=1.25, hole under a through-thickness bending gradient
  near a cantilever free end).
- The moment-gradient explanation for the residual gap is speculative
  engineering reasoning, not a demonstrated or literature-supported
  cause, and should not be cited as settled in any downstream summary
  (e.g. LinkedIn post, portfolio writeup) without that caveat.
- Only two mesh levels were run for convergence (justified by the
  decisively small 0.28% change, but a third level was not generated).
- The 2D iterative solver limitation (stalling on locally-refined
  meshes) means all hole-stage results used the direct solver only;
  no cross-solver corroboration exists for this stage, unlike Stage 1's
  cross-discretization check.

---

# Progressive Geometry Stage 3 — Shoulder Fillet

## Geometry

A distinct interpretation from the eventual full-bracket flange-
junction fillet, chosen deliberately: a **shoulder fillet** — a
step change in cross-section thickness on an otherwise straight
cantilever, blended by a fillet radius. This isolates the fillet
stress-concentration mechanism against one of the most thoroughly
documented SCF problems in mechanical engineering (Peterson's
stepped-flat-bar-in-bending chart), before the flange-junction fillet
is introduced later as part of the full angled-bracket geometry
(Stage 5), where it cannot be isolated from the gusset/bend geometry.

- 60 mm overall length, 40 mm constant out-of-plane width
- Root section (x=0–30mm): depth D = 6 mm
- Tip section (x=30–60mm): depth d = 4 mm
- Symmetric opposite shoulder fillets (both top and bottom faces,
  centroid preserved through the transition — required by the
  Peterson reference chart, which assumes this symmetry)
- Fillet radius: **r = 0.95 mm** (see construction history below for
  why this differs from the originally locked r = 1.0 mm target)
- Same material, load (F = 235.44 N transverse at tip), coordinate
  system, and encastre-root BC philosophy as Stage 1

**Non-dimensional parameters, as-built:** D/d = 1.5, r/d = 0.2375
(target 0.25), h/r = 1.0526 (target 1.0), where h = (D−d)/2 = 1.0mm —
all within the Peterson/ESDU chart's demonstrated small-r/d
calibration range, confirmed via a dedicated specification gate
before any geometry was built (checking the chart's actual parameter
definitions and range against secondary literature, not assumed).

### Verification reference

Peterson's "Opposite Shoulder Fillets in a Flat Bar" bending chart.
σnom = 6M/(t·d²), evaluated at the same beam-theory moment already
verified in Stage 1/2 (M(x=30) = 7063.2 N·mm, giving σnom = 66.22 MPa
— identical to the Stage 1/2 gross-section reference at this
location, confirming self-consistency of the mapping).

**Locked prediction (computed before any FEA, at h/r=1.0 target):**
Kt = 1.429, predicted peak σₓₓ ≈ 94.6 MPa.

### Construction history — why r = 0.95mm, not the locked r = 1.0mm

1. **r = 1.0mm exactly equals the riser height** h = (D−d)/2 = 1.0mm.
   The fillet's tangent point coincides exactly with the adjacent
   sharp-corner vertex — a genuine BRep degeneracy. Gmsh's OCC fillet
   solver fails outright ("Could not compute fillet"). This is a
   numerical-tolerance issue, not an engineering change — flagged and
   confirmed as such before adjusting.
2. **r = 0.99mm** succeeds numerically but leaves a 0.01mm residual
   sliver face (h−r) smaller than the finest intended mesh size
   (0.15mm), forcing pathological element aspect ratios (up to
   AR=62.77) in its immediate neighborhood.
3. `gmsh.model.occ.healShapes()` was attempted to remove the sliver.
   It did not repair the geometry in place — it duplicated the
   surface set (12→22 surfaces, orphaning the original topology
   alongside a new one) and shifted the bounding box by ~0.015mm.
   Reverted entirely rather than used.
4. **r = 0.95mm** (final): sliver grows to 0.05mm — large enough to
   mesh cleanly without healing. This is the geometry used for all
   subsequent work in this stage.

### Mesh methodology and a second defect, also fixed

Unstructured tetrahedra (C3D10), Gmsh OCC box-union + fillet, same
family of approach as Stage 2. Local refinement originally attempted
via a `Box` field covering the full cross-section over the transition
length — this produced a 1.44M-element mesh (far larger than anything
previously solved in this project) because it refined the *entire*
6×40×7mm slab around the transition, not just the fillet surfaces
themselves. Replaced with a `Distance`+`Threshold` field pair anchored
specifically on the two fillet (Cylinder-type) surfaces, identified
robustly by OCC surface type + z-sign — reducing the baseline mesh to
~11,000 elements.

This Distance+Threshold field's default **linear** interpolation was
then found to introduce a second, smaller defect: a hard slope
discontinuity at `DistMax` caused elements near that boundary to
stretch (AR up to 28.5). Switching to the field's built-in **Sigmoid**
interpolation (a smooth, continuous-slope transition) resolved this
decisively — mean aspect ratio dropped from 5.93 to 2.22, elements
with AR>20 dropped from 76 to a single isolated outlier, at every
subsequent mesh level.

## Mesh Convergence Study

| Level | FILLET_SIZE | Elements | Equations | Top σₓₓ (MPa) | Top Kt | Bottom σₓₓ (MPa) | Bottom Kt |
|---|---|---|---|---|---|---|---|
| Coarse | 0.30 mm | 10,819 | — (not solved to equilibrium check) | 86.54 | 1.307 | −91.49 | 1.382 |
| Medium (baseline) | 0.15 mm | 37,684 | 207,738 | 93.82 | 1.417 | −94.97 | 1.434 |
| **Fine (accepted)** | **0.075 mm** | **156,219** | **784,701** | **95.99** | **1.450** | **−94.98** | **1.434** |

**% change, medium → fine:**
- Bottom surface: +0.01% — converged, decisively.
- Top surface: +2.31% — marginally above the project's 2% criterion,
  and still trending in the same direction as the coarse→medium step
  (not yet clearly plateaued).

**Two further refinement attempts (FILLET_SIZE = 0.06mm and 0.05mm)
were mesh-quality-verified (both clean: zero inverted elements, max
AR < 15, consistent with the sigmoid-field fix) but rejected as
solver candidates on resource grounds** — extrapolated at ~1.22M and
~1.7M equations respectively, at or above the ~1.18M-equation point
that OOM-killed CalculiX's direct solver during Stage 2's hole study.
A quantitative mesh-growth diagnostic (comparing the fine and 0.06mm
meshes, bucketed by distance to the fillet surface) attributed ~70–73%
of the excess node/element growth to the immediate fillet-surface
refinement band itself (expected — halving element size on a
fixed-area curved surface roughly quadruples the elements needed to
resolve it) and ~27–29% to the fixed 1.5mm transition width, with
far-field mesh density essentially unaffected. This is recorded for
reference should further refinement be attempted in a future session,
but was not acted on in this one.

**Decision: the fine level (0.075mm, 156,219 elements) is accepted as
the baseline result for this stage**, with the top-surface 2.3%
residual documented explicitly as an open, unresolved limitation —
not silently treated as converged, and not treated as a failed
verification (the result itself is physically and numerically sound
by every other check performed).

## Baseline Solve — Full Sanity Check

All results below are from the accepted fine-level mesh (261,724
nodes / 156,219 elements / 784,701 equations, direct spooles solver,
44.1–44.5s across two identical-deck runs).

**Equilibrium check** (added via a minimal `*NODE FILE U,RF` deck
change, rerun on the identical mesh/BC/load — reproduced the same
784,701-equation solve, confirming no unintended model change):
- Applied load: Fz = −235.44 N
- Root ΣFz = +235.4397 N (opposite sign, equal magnitude, as required)
- **Relative equilibrium error: ≈0.0001%**
- Root ΣFx = +0.0002 N, ΣFy = −0.0002 N (both ≈0.0001% of applied
  load — no spurious lateral/axial reaction from the `*DISTRIBUTING`
  coupling)

**Displacement:**
- Peak |U| = 0.4120 mm (tip edge, x=60/y=20/z=2)
- Dominant Uz = −0.4112 mm; small Ux=0.0253mm (Poisson/bending
  coupling); Uy≈0
- Physically sane: smaller than Stage 1's ~1.05mm tip deflection
  under the same load/span — expected, since this geometry averages
  thicker (6mm root vs. Stage 1's uniform 4mm) over half its length

**Stress vs. Peterson reference (94.63 MPa predicted):**

| | σₓₓ | Observed Kt | vs. Peterson |
|---|---|---|---|
| Top fillet | +95.99 MPa | 1.450 | +1.44% |
| Bottom fillet | −94.98 MPa | 1.434 | +0.37% |

Both surfaces agree with the classical closed-form reference to
within ~1.5% — a substantially cleaner result than Stage 2's hole
study achieved against any available literature source.

**Sanity check status: 8/8 closed** — element type/count, node count,
material properties, root BC, load magnitude/direction/location,
reaction-force equilibrium, displacement/stress scale, and
rigid-body-motion (no singular-pivot solver error, and now also
confirmed explicitly via exact equilibrium) all verified.

## Verification vs. Validation (Stage 3)

As with Stages 1–2, this stage is **verification only**. Unlike
Stage 2, this stage's classical reference (Peterson's shoulder-fillet
chart) is a strong, closely-matching quantitative check (<1.5% on
both surfaces) — but this remains a check against another idealized
model, not against physical test data, and no physical validation is
claimed.

## Limitations (Stage 3)

- The fillet radius used (r=0.95mm) differs from both the originally
  locked verification target (r=1.0mm, infeasible due to a genuine
  BRep degeneracy) and the eventual full-bracket fillet spec
  (r=3.0mm) — this stage verifies the pipeline and the general
  shoulder-fillet SCF mechanism at a chart-validated radius, not the
  bracket's final geometry.
- The top-surface stress result has **not** been demonstrated
  mesh-converged: medium→fine changed by +2.31%, marginally above
  the project's 2% criterion, and still trending upward. The
  bottom-surface result (+0.01% change) is converged. This asymmetry
  between two nominally-symmetric surfaces is itself unexplained —
  documented as an open item, not resolved.
- Two further refinement levels were mesh-quality-verified but never
  solved, due to resource constraints relative to this project's
  documented direct-solver memory ceiling. A path to a safer further
  refinement (reducing `TRANSITION_WIDTH` proportionally, per the
  mesh-growth diagnostic) was identified but not attempted.
- The Peterson chart's underlying calibration is presumed to assume a
  spatially uniform remote bending moment; this stage's cantilever
  geometry has a real (if small, ~3.3% over the fillet's immediate
  neighborhood) moment gradient along x that the classical chart does
  not account for. Unlike Stage 2, this was not found to produce a
  large discrepancy — the <1.5% agreement suggests this effect is
  genuinely small here, but it was not independently isolated or
  quantified.

---

# Progressive Geometry Stage 4a — Right-Angle Frame (Sharp Corner, No Gusset)

## Geometry

A deliberate expansion of scope: rather than study the gusset on a
straight cantilever (the original plan), the frame/bend geometry was
introduced here, one stage earlier than originally planned, so the
gusset (Stage 4b) could be studied on a closer analogue of the final
bent bracket. This sub-stage (4a) establishes and verifies the bent
frame *without* a gusset, as the baseline the gusset's effect is
measured against.

- **Flange A** (fixed): x ∈ [0, 60mm], z ∈ [−2, 2] (centered),
  40mm constant width (y). Root (x=0) is the fixed encastre boundary.
- **Corner**: sharp 90° reentrant junction at (x=58, z=2) — no
  fillet, no gusset in this sub-stage, by design.
- **Flange B** (loaded): x ∈ [58, 62], z ∈ [0, 60], same 40mm width
  and 4mm thickness. Free tip at z=60.
- Same material (Al 7075-T6) and cross-section as every prior stage.

**Load:** F = 235.44 N applied at the Flange B tip, in the **−x
direction** — chosen specifically so the whole problem stays planar
(pure bending in Flange B, combined bending+axial in Flange A), a
direct two-member extension of Stage 1's single-flange bending case
rather than a new loading mode.

**A known construction detail, flagged rather than hidden:** the
outer (convex) corner comes out as a small unmitered step, since the
geometry is built from two axis-aligned boxes rather than a mitered
join. This does not affect verification, since every quantitative
comparison point sits away from that region.

### Verification strategy — two independent closed-form checks

1. **Internal force diagram (exact statics)**: axial force and
   bending moment derived directly from free-body equilibrium, not a
   beam-theory approximation. Flange B: M(s) = F·s (s = distance from
   tip), pure bending. Flange A: constant M = F·L = 14,126.4 N·mm and
   constant axial N = −F = −235.44 N (compressive) along its entire
   length.
2. **Tip deflection (Castigliano / unit-load method)**: closed-form
   strain-energy calculation combining Flange B's bending, Flange A's
   bending, and Flange A's (much smaller) axial contribution.

**A genuine sign error, caught and corrected via the FEA comparison
itself — documented transparently, not hidden:** the original
hand-derivation had the Flange A axial force sign backwards (called
it tensile; it is compressive, confirmed independently via free-body
equilibrium and via "bend-closing" physical intuition — pulling the
tip of the vertical flange toward the fixed flange closes the frame's
interior angle, putting the concave/inside surface in compression and
the convex/outside surface in tension, the same behavior as bending
an elbow) and the tension/compression face labels swapped accordingly.
Magnitudes were unaffected; only the labeling was wrong. Corrected
values, locked before further FEA comparison:

| | σ_bending | σ_axial | σ_total |
|---|---|---|---|
| Flange A, concave (z=+2, inside) | −132.4350 | −1.4715 | **−133.9065 MPa** |
| Flange A, convex (z=−2, outside) | +132.4350 | −1.4715 | **+130.9635 MPa** |

Flange B target magnitudes (concave = compression, convex = tension):
0 (tip) → 33.11 → 66.22 → 99.33 → 132.44 MPa at s=0/15/30/45/60mm.

**Locked tip deflection target:** 4.4342 mm (in the −x direction).

**A known, expected, unavoidable limitation, stated before any FEA
was run:** the sharp reentrant corner produces a classical FEA stress
singularity — not a modeling error, a known mathematical feature of a
zero-radius interior corner. Quantitative comparison points are all
chosen away from the corner (mid-span stations); corner behavior is
reported separately and excluded from the convergence criterion.

## Mesh Convergence Study

Three levels (unstructured tetrahedra, Distance+Threshold sigmoid
field anchored on the reentrant corner edge):

| Level | Corner size | Elements | Tip Ux (mm) | % change |
|---|---|---|---|---|
| Coarse | 1.0mm | 5,215 | −4.1986 | — |
| Medium (baseline) | 0.5mm | 8,124 | −4.2128 | — |
| Fine | 0.25mm | 20,910 | −4.2180 | **0.12%** |

Every mid-span stress comparison point (Flange A both faces at x=30;
Flange B both faces at s=15/30/45) changed by **at most 1.61%**
between medium and fine — **decisively converged**, well under the
2% criterion, and better-behaved than either Stage 2 or Stage 3.

**Corner (s=60) behavior, confirmed non-convergent by direct
evidence, excluded from the convergence criterion by design:** peak
values at the corner do not approach the idealized target with
refinement (coarse: −72.4%/−93.6% vs. theory; fine: −58.1%/−109.2%,
convex face even changing sign) — the expected fingerprint of a
genuine stress singularity, confirmed rather than assumed.

## Baseline Results (fine level, accepted)

- **Equilibrium**: essentially exact at all three levels
  (0.0004–0.0017% error)
- **Tip deflection**: −4.2180mm FEA vs. −4.4342mm target, **+4.88%**
  — consistent with the same order of 3D/boundary-effect gap seen in
  Stage 1 (−5.2%), not a red flag
- **Mid-span stress, both flanges, both faces**: converged (<2%) and
  matches the corrected analytical reference to within **0.02–3.30%**
  — strong agreement across every valid comparison point

## Established vs. Not Established (Stage 4a)

**Established:**
- The two-member frame's internal force diagram and tip deflection
  are verified via independent closed-form statics/Castigliano
  checks, matching FEA to within 0.02–3.30% (stress) and 4.88%
  (deflection).
- Mesh convergence at every mid-span comparison point is decisive
  (≤1.61% change), independently of the corner region.
- The reentrant-corner singularity is confirmed non-convergent by
  direct three-level evidence, not merely predicted.

**Not established:**
- No finite, mesh-independent stress value exists at the sharp
  reentrant corner for this geometry; any single reported corner
  value is not physically meaningful.

---

# Progressive Geometry Stage 4b — Right-Angle Frame WITH Gusset

## Geometry and Assumptions

Adds a flat triangular-prism gusset to the exact Stage 4a frame
geometry, bridging the concave (inside-of-bend) faces:

- **Leg 1**: along Flange A's concave face (z=2), from the corner
  (x=58) to x=38 (20mm leg) — flush/coplanar with Flange A's top
  surface.
- **Leg 2**: along Flange B's concave face (x=58), from the corner
  (z=2) to z=22 (20mm leg) — flush/coplanar with Flange B's inner
  surface.
- **Hypotenuse**: connects (38,y,2) to (58,y,22) directly — the only
  genuinely new surface feature this geometry introduces, since both
  legs are coplanar with existing flange faces.
- Gusset thickness: 4mm, matching the flanges — the simplest choice,
  avoiding a second thickness-mismatch problem stacked on the one
  under study.
- Full 40mm width, sharp edges at both the leg/flange junctions and
  the toe — **no fillet, no taper at the gusset termination**, by
  deliberate design choice for this isolated verification sub-stage.
  This is a simplification, explicitly not revisited as part of this
  stage (a tapered/filleted redesign is a natural follow-up, not
  performed here).
- Same material, load (F=235.44N at Flange B tip, −x direction),
  coordinate system, and root BC as Stage 4a.

**Verification strategy, scoped from the outset:** no independent
analytical/chart reference exists for a gusseted corner's local
stress field. This stage's deliverable is explicitly comparative
(FEA-vs-FEA against the Stage 4a baseline), not closed-form, with one
exception: total tip deflection is expected to *decrease* (a
directly predictable, if not quantitatively closed-form, sanity
check).

### A prediction, checked before being relied upon

Because both gusset legs are flush/coplanar with the existing flange
faces, they should disappear as distinct surface features after the
boolean union, fully enclosing Stage 4a's reentrant corner. **This
was confirmed by direct geometric inspection** (no distinct edge
remains at the old corner location, x=58/z=2) — not assumed.

## Mesh methodology

Same family of approach as 4a: unstructured tetrahedra, Gmsh OCC
union (three solids), sigmoid Distance+Threshold refinement — now
targeted at the hypotenuse **face** (a `Plane`-type surface; an
initial attempt to target a corresponding "edge" failed, since a
straight edge cannot simultaneously span the full width and be
diagonal in x-z — only the flat sloped face has both properties,
the same correction Stage 3 needed when it moved from edge-based to
surface-based field targeting for the curved fillet).

## Baseline Results vs. Stage 4a

**Equilibrium**: essentially exact (−0.0021% error).

**Tip deflection**: −2.4628mm (4b) vs. −4.2180mm (4a fine level) —
a **41.6% reduction**. This is a global response quantity, verified
via near-exact equilibrium on the solved model, and is the primary,
high-confidence quantitative result of this stage: the gusset has a
substantial global stiffening effect.

**Valid mid-span comparisons** (stations clear of the gusset
footprint — Flange B s=15, Flange A x=30): stress still matches the
original 4a beam-theory targets to within **0.27–3.76%**, confirming
the measured stress disturbance associated with the gusset
termination is spatially localized and does not disturb the far-field
bending behavior. (Stations at s=45/60 on Flange B fall inside or
adjacent to the gusset's own footprint and are not valid comparisons
against the unreinforced 4a target — they were excluded from this
comparison rather than misreported.)

## Targeted Toe Convergence Study

The whole-model peak stress in the baseline 4b solve (−274.4 MPa) was
located not at the old corner (confirmed eliminated) but at the
**gusset's far toe** (x=38, z=2) — the line where the un-tapered,
un-filleted gusset termination meets the plain flange, a classic
sharp-transition stress-riser geometry. A dedicated three-level
convergence study, refining specifically at this toe edge (geometry,
material, BCs, and loading held fixed throughout), was performed to
determine whether this value is a bounded, mesh-convergent stress
concentration or a non-convergent singularity.

| Level | Toe size | Elements | Equations | Peak σₓₓ at toe | % change |
|---|---|---|---|---|---|
| Baseline | ~0.5mm (implicit) | 44,759 | 223,995 | −274.39 MPa | — |
| toe_fine | 0.15mm | 83,082 | 393,591 | −422.01 MPa | +53.8% |
| toe_intermediate | 0.10mm | 139,252 | 633,180 | −483.22 MPa | +14.5% |

**The peak toe stress increases monotonically across the tested mesh
levels without evidence of convergence** — failing the 2% criterion
by a wide margin, with no plateau across three genuinely different,
independently solved mesh densities. (A further-refined level,
0.05mm, was mesh-quality-verified — zero inverted elements, clean
aspect ratios — but rejected as a solver candidate: at 582,770 nodes
it extrapolates to ~1.75M equations, above the ~1.18M-equation point
that OOM-killed CalculiX's direct solver in Stage 2. An intermediate
0.06mm attempt also overshot its intended target due to the same
nonlinear mesh-growth behavior documented in Stage 3, and was
likewise not solved.)

**Critically, this is spatially confined, confirmed by a separate
check moving away from the toe along Flange A:**

| Distance from toe | Baseline | toe_fine | toe_intermediate |
|---|---|---|---|
| 0mm (at toe) | −164.7 | −153.3 | −152.5 |
| 2mm | −143.3 | −140.9 | −143.3 |
| 4mm | −141.1 | −140.9 | −141.0 |
| 6mm | −138.5 | −138.4 | −138.4 |
| 8mm | −136.7 | −136.6 | −136.5 |
| 13mm | −133.3 | −133.2 | −133.3 |
| 18mm | −132.2 | −132.2 | −132.2 |

Every station 2mm or more from the toe is converged to within a
fraction of a percent across all three levels, decaying smoothly
toward the plain-flange far-field value (~−132 MPa, matching 4a's
Flange A convex-side result) — a Saint-Venant-type decay pattern,
structurally identical to the boundary effect already characterized
at Stage 1's root.

## Established vs. Not Established (Stage 4b)

**Established:**
- The gusset reduces tip deflection by 41.6%, verified via near-exact
  equilibrium on the solved model.
- The original Stage 4a reentrant-corner singularity is eliminated as
  a surface feature (confirmed by direct geometric inspection).
- Mid-span stress at valid (non-gusset-footprint) comparison stations
  is unchanged from 4a, confirming the measured stress disturbance
  associated with the gusset termination is spatially localized.
- The gusset's flat, un-tapered, un-filleted termination (as built in
  this isolated verification geometry) introduces a **new stress
  concentration at its toe that increases monotonically across the
  tested mesh levels without evidence of convergence** — demonstrated
  directly via a three-level targeted convergence study, not
  inferred.
- Stress is fully converged (sub-percent) everywhere more than ~2mm
  from the toe, confirming this behavior is spatially confined and
  does not compromise the rest of the model.

**Not established:**
- No finite, mesh-independent peak stress value exists at the gusset
  toe for this geometry; any single reported number there is not
  physically meaningful and must not be used for a factor-of-safety
  calculation.
- Whether a tapered or filleted gusset termination would resolve this
  secondary stress concentration was **not tested in this stage** —
  flagged as a natural follow-up question raised by this result, not
  investigated here. No taper/fillet redesign was introduced.
- This isolated two-flange verification geometry is not the final
  bracket geometry; the gusset's behavior in the actual angled,
  filleted, holed full bracket (Stage 5) is not established by this
  result alone, though the underlying mechanism (untapered
  terminations concentrate stress) should transfer.

## Engineering Takeaway and Limitation

This stage demonstrates, numerically rather than by assumption, why
real gusset designs taper or fillet their terminations: the gusset,
**as built in this simplified isolated-verification form**, trades
one singularity (4a's sharp reentrant corner) for a different one
(the un-tapered toe), rather than eliminating stress concentration
altogether, even while delivering a substantial and genuine (41.6%)
global stiffness benefit. The verification process — specifically,
the targeted convergence study — caught this rather than allowing a
single baseline-mesh reading to be misreported as a converged,
actionable stress value.

## Figures (Stage 4a/4b)

1. `stage4_frame_mesh_visualization.png`, `stage4_frame_corner_slice.png`
   — 4a geometry/mesh overview and reentrant-corner cross-section
2. `stage4_frame_deformed.png` — 4a deformed shape, used to visually
   resolve the axial/bending sign question during the sign-error
   investigation
3. `stage4b_gusset_mesh_overview.png`, `stage4b_gusset_corner_slice.png`
   — 4b geometry/mesh overview and gusset cross-section
4. `stage4b_gusset_full_contour.png` — full-model von Mises contour
   (5× deformed); von Mises used rather than a single raw stress
   component so both flanges display correctly, since each bends
   about a different axis
5. `stage4b_toe_nonconvergence.png` — peak toe stress vs. mesh
   refinement (primary non-convergence evidence)
6. `stage4b_toe_decay_profile.png` — stress vs. distance from the
   toe, overlaid across all three mesh levels (primary evidence of
   spatially-confined behavior)

---

# Progressive Geometry Stage 5 — Full L-Bracket Integration

## Geometry and Assumptions

Combines every prior geometric feature onto the bent Stage 4a/4b
frame, forming the final bracket geometry:

- **R3mm interior flange-junction fillet** at the true 90° bend
  (Option A construction: fillet built into the base frame first,
  gusset legs attach at the fillet's tangent points — 17mm effective
  leg length, rather than the original 20mm sharp-corner legs used
  in Stage 4b's isolated gusset study).
- **Quadrilateral gusset cross-section**, not triangular — the R3
  fillet consumes the original sharp-corner vertex the Stage 4b
  gusset attached to. Vertices: Toe (38,2), tangent point on Flange A
  (55,2), tangent point on Flange B (58,5), far corner (58,22). Full
  40mm frame width (matching every prior stage).
- **R1mm toe fillet** at the gusset toe (38,2) — a deliberate response
  to Stage 4b's finding that an un-tapered gusset termination produces
  a non-convergent stress concentration. This stage tests whether a
  fillet, rather than the untested taper/fillet redesign Stage 4b
  explicitly deferred, resolves that behavior (see Toe Fillet Stress
  Convergence, below).
- **Ø5mm through-hole** in Flange A at (x=20, y=20), matching Stage
  2's hole specification and location logic (clear of root
  Saint-Venant effects, mid-width).
- Same material (Al 7075-T6), load (F=235.44N at Flange B tip,
  −x direction), and root encastre BC as every prior stage.
- **A known, documented, intentional design consequence**: a small
  "corner-bite void" (~102.7mm³) exists between the R3 fillet's arc
  and the gusset's flat tangent-to-tangent chord face. This is a
  direct geometric consequence of keeping the gusset flat rather than
  lofting it to match the cylindrical fillet — not a defect, and not
  revisited in this stage.
- Final locked geometry: `full_bracket_geometry_final.brep`, mass =
  26,859.58 mm³, single valid volume, verified clean via BRep
  diagnostics (no slivers, no duplicate/non-manifold surfaces).

**Unlike Stages 1–4b, no single closed-form or chart reference exists
for this combined geometry.** This stage is explicitly framed as an
*integration and consistency study*: checking equilibrium, checking
global deflection against the Stage 4a/4b bracketing values (not
assumed greater or lesser), and checking that feature-level behavior
(hole SCF, fillet SCF) remains consistent with each feature's
isolated characterization in Stages 2–4b, rather than claiming a new
independent verification.

## Toolchain Deviation — Netgen Replaces Gmsh for This Project, Locked Going Forward

**This is a material methodology change from every prior stage,
recorded explicitly rather than silently adopted.**

Every previous stage used Gmsh (Python API) for meshing and CalculiX
2.21 (apt-installed) as the solver. For Stage 5, meshing was
performed with **Netgen** (`python3-netgen`, v6.2.2401) instead of
Gmsh, and solving with **CalculiX 2.23** (installed via conda-forge
into a dedicated `ccx223` environment) instead of 2.21. Following
discussion, **this is now the locked methodology for the remainder of
the project** (Stage 6 onward), not a Stage-5-only exception.

### Why: an extensively investigated, unresolved Gmsh defect

Meshing Stage 5's geometry with Gmsh and solving with C3D10 (10-node
quadratic tetrahedra) produced a persistent, non-physical failure:*ERROR in e_c3d: nonpositive jacobian
determinant in element [tag] 
affecting a small fraction of elements (~0.3–1.4% depending on mesh
density), consistently concentrated at the tail of Gmsh's internal
element-numbering range, near — but not exclusively at — the tightly
curved fillet and hole-bore regions.

**Sixteen distinct mechanisms were investigated and ruled out, each
with direct evidence, not assumption:**

1. Curvature at known features — isolated single-element and
   83-element neighbor-patch extracts passed; the identical element
   failed only in the full-model context (a result that should be
   mathematically impossible for a per-element geometric check).
2. C3D10 midside-node ordering/connectivity — confirmed already
   correct; a "fix" was applied, had zero effect, and was reverted.
3. Toe-fillet relief mesh-sizing fields — ruled out; the same defect
   appeared in meshes generated before these fields existed.
4. Mesh periodicity from OCC extrude operations — checked directly
   via `gmsh.model.getPeriodic()`; none exists anywhere in the model.
5. Corner-node degeneracy / inverted neighbor elements — all 462
   failing elements' linear (corner-only) volumes checked directly:
   all positive, comparable in scale to neighboring good elements.
6. Mixed element types / section controls / reduced integration —
   confirmed only one `*ELEMENT` block exists, correctly declared,
   no conflicting keywords present.
7. `*COUPLING`/`*KINEMATIC` reference-node dependency — rebuilt with
   coupling removed entirely (direct nodal load instead): failure
   persisted, with a different but still tail-concentrated element set.
8. Duplicate element tags — directly counted: zero duplicates across
   157,467 elements.
9. Direct vs. iterative solver — both failed, same signature.
10. Coordinate-rounding node collisions (6-decimal `.inp` export) —
    checked directly: zero collisions in the densest mesh.
11. Coordinate precision itself, independent of collisions — rewrote
    the `.inp` export at full floating-point precision: failure
    persisted identically.
12. CalculiX version — reproduced identically on both 2.21 and 2.23
    against the exact same `.inp` file.
13. Meshing algorithm — switched Gmsh's 3D algorithm from Delaunay to
    HXT; produced a mesh with zero bad elements by Gmsh's own quality
    metric, but CalculiX still failed, at the new mesh's own tail.
14. Model scale — a deliberately coarsened mesh (61,718 vs. ~157,000
    elements) failed at the same relative signature, ruling out a
    simple large-model memory/precision threshold. (Separately
    confirmed via the CalculiX user forum that a different, unrelated
    large-model failure mode exists near 1M+ DOF due to a documented
    32-bit integer limitation in the Spooles solver source — not
    applicable here, since the smallest failing case was 277,971 DOF.)
15. Gmsh's internal mesh-optimization pass — disabled entirely:
    quality got dramatically worse and more broadly distributed,
    confirming the optimizer performs necessary work rather than
    causing this specific, narrow defect.
16. A documented forum precedent describing an identical
    isolated-vs-assembled Jacobian-failure symptom, resolved by the
    original poster via remeshing with different internal numbering —
    directly informed and motivated mechanism #13, which was tested
    and also failed to resolve the issue.

**Established conclusion:** this is a genuine limitation or edge case
in Gmsh's Delaunay-family tetrahedral generation (both default and
HXT algorithms) for this geometry's tightly-refined curved-feature
regions. It is not a defect in the project's geometry (independently
verified clean via BRep diagnostics), not a `.inp`-writing bug
(multiple precision/connectivity fixes had no effect), not a CalculiX
version issue (reproduced on 2.21 and 2.23), and not a simple scale
effect. The repeatedly-observed context-dependence (identical node
coordinates producing different Jacobian verdicts depending on
whether the element is isolated or embedded in the full mesh) strongly
suggests Gmsh silently produces subtly different node placement in
the full-model context than in isolated extracts — the exact
mechanism was not further identified, given the exhaustiveness of the
above investigation and the availability of a working alternative.

### The Netgen-based workaround

Rather than continue debugging Gmsh, meshing was rebuilt from scratch
using Netgen, with a deliberately independent method chosen to
sidestep every failed Gmsh mechanism above, not merely to try a
different tool:

1. Load the verified BRep directly via `netgen.occ.OCCGeometry`.
2. Generate an order-1 (4-node, corner-only) linear tetrahedral mesh
   using Netgen's own native algorithm — entirely independent of Gmsh.
3. **Verify every element has strictly positive linear volume** via
   direct signed-volume computation in Python, before proceeding —
   a hard mathematical check, not a solver-side hope.
4. **Promote to C3D10 (quadratic) manually**, placing every midside
   node at the exact geometric midpoint of its parent edge, rather
   than using either tool's internal curving/order-2 logic. This
   exploits a proven identity: a straight-sided quadratic tet's
   Jacobian is identical everywhere to its linear parent's Jacobian.
   Since step 3 already guarantees a positive linear-parent volume,
   this **guarantees a valid quadratic element by construction** —
   sidestepping the unresolved Gmsh/curving question entirely rather
   than resolving it.
5. Local mesh-size control via Netgen's `SetFaceMeshsize()` on the
   three feature faces (hole bore, R3 fillet, R1 toe fillet),
   identified robustly by area fingerprint — the same identification
   method used with Gmsh throughout this project.

This mesh generation and quadratic-promotion pipeline was successful
at all three levels solved (see below), with zero Jacobian failures
at any level — confirming the workaround, not merely masking the
original defect.

## Mesh Convergence Study

Three refinement levels, controlled via face-level mesh sizing at the
hole bore, R3 fillet, and R1 toe fillet (global background size held
constant at 6.0mm across all levels):

| Level | Nodes | Elements | Equations | Equil. error | Tip \|Ux\| (mm) | Peak VM near toe (MPa) |
|---|---|---|---|---|---|---|
| 1 | 67,720 | 43,686 | 202,896 | 0.0002% | 2.6367 | 226.21 |
| 2 | 132,013 | 86,232 | 395,775 | 0.0008% | 2.6397 | 230.01 |
| 3 | 270,250 | 179,242 | 810,486 | 0.0008% | 2.6401 | 230.55 |

**% change, L2→L3:** displacement 0.02%, peak near-toe stress 0.24% —
both decisively converged, with the stress delta shrinking sharply
from the L1→L2 step (1.68%), the signature of genuine convergence to
a finite value rather than the ever-growing, non-convergent pattern
documented at Stage 4b's sharp (un-filleted) toe.

**A fourth level was considered and deliberately not run.** Continuing
the same refinement trend would have required ~1.6–1.7M equations,
exceeding this project's documented ~1.18M-equation direct-solver
ceiling; an iterative solver could reach that range but would
introduce a solver-methodology variable into what is otherwise a
single-solver convergence study. Given L1–L3 already show a clean,
sharply shrinking convergence trend, L3 was accepted as sufficient —
a defensible but explicit engineering judgment call, not a default.

## Baseline Results (Level 3, accepted)

**Global equilibrium**: 0.0008% error — excellent at all three levels,
confirming the `*COUPLING`/`*KINEMATIC` setup (with the corrected
`1, 3` DOF specification for solid elements) is correct across mesh
densities.

**Tip deflection**: −2.6401mm, bracketed as expected between Stage 4a
(−4.2180mm, no gusset) and Stage 4b (−2.4628mm, gusset with 20mm
sharp-corner legs) — sitting closer to 4b, consistent with Stage 5
retaining the same dominant gusset-stiffening mechanism, offset
slightly by the reduced 17mm effective leg length (fillet consumes
3mm of each leg) and the added hole. Reported as an observed,
internally consistent FEA comparison against the established
bracketing values, not a predetermined target.

## Stress Results — Two Distinct Governing Locations

Direct inspection of the full nodal stress field (not just the
toe-local search used for the convergence table above) identified
**the global peak stress at 286.22 MPa, located at the Flange A hole
edge — not at the toe fillet**:

| Feature | Peak von Mises (L3) | Character |
|---|---|---|
| Flange A hole edge | **286.22 MPa** (global peak) | Classic bore concentration, located at the two points ~90°/270° around the bore relative to the load axis (x≈20±0.4, y≈20±2.5, matching the bore radius) — the textbook Kirsch-type pattern, also visually confirmed as a clean two-lobe contour |
| Gusset toe fillet | 230.55 MPa (local peak) | Forms a continuous ridge along the full 40mm width, not a point concentration — confirmed by direct nodal inspection (a straight line of near-peak nodes spanning the width) and visually confirmed via a top-down contour view |

**The hole, not the toe, governs Stage 5's peak stress** — a result
this stage's structure did not originally set out to test for, but
was surfaced directly by the L3 full-field inspection rather than
assumed from the narrower toe-focused search that motivated the
convergence study.

### Toe fillet: does the R1 fillet resolve Stage 4b's non-convergence?

Stage 4b found that an un-tapered, un-filleted gusset toe produces a
stress concentration that grows without bound under mesh refinement —
a genuine singularity. Stage 5 adds an R1mm fillet at the same
location specifically to test whether this resolves that behavior.

**Result: yes, on the evidence of this study.** The toe stress
converges cleanly across L1→L3 (1.68% then 0.24% change), in sharp
contrast to Stage 4b's monotonic, non-convergent growth. This is
physically expected — a finite fillet radius regularizes what would
otherwise be a sharp re-entrant corner — but Stage 5 is the first
point in this project where that expectation is confirmed
numerically for this bracket's specific geometry, rather than merely
argued as a general principle.

## Established vs. Not Established (Stage 5)

**Established:**
- Global equilibrium is verified to ≤0.0008% error at all three mesh
  levels, confirming the coupling/BC/load setup is correct.
- Tip deflection is converged (0.02% L2→L3) and falls consistently
  between the Stage 4a and 4b bracketing values.
- The R1mm toe fillet resolves Stage 4b's non-convergent toe
  singularity — the toe stress is now mesh-convergent (0.24% L2→L3),
  demonstrated directly via the same three-level methodology used to
  characterize Stage 4b's non-convergence, not merely predicted.
- The Flange A hole edge, not the gusset toe, is the bracket's
  governing (global peak) stress location — identified via direct
  full-field nodal inspection, not assumed from a narrower search.
- The Gmsh C3D10 Jacobian failure is a genuine, reproducible defect
  independent of geometry validity, `.inp` precision/connectivity, and
  CalculiX version, established via sixteen independently tested and
  ruled-out mechanisms.
- The Netgen-based straight-sided quadratic promotion produces
  mathematically guaranteed valid C3D10 elements by construction, and
  eliminated the Jacobian failure entirely across all three mesh
  levels solved.

**Not established:**
- No independent analytical, chart, or published reference exists for
  this combined geometry's stress field — this stage is an internal
  consistency and convergence study, not a verification against an
  external ground truth, unlike Stages 2 and 3.
- The hole-edge peak (286.22 MPa) has not been independently checked
  against a Peterson-style Kt reference for this specific combined
  loading state (bending + local bracket flexibility, rather than the
  simpler cantilever bending case characterized in Stage 2) — a
  natural follow-up, not performed in this stage.
- A fourth mesh level was not run; the L1–L3 convergence trend is
  clean and decisively shrinking, but this is a smaller evidentiary
  base than, for example, Stage 4b's toe non-convergence study (which
  itself used three levels to establish non-convergence, an inherently
  easier claim than establishing convergence).
- The exact mechanism behind Gmsh's context-dependent Jacobian failure
  (why isolated and embedded instances of the same element differ)
  was not identified, only its existence and independence from every
  tested alternative explanation.

## Limitations (Stage 5)

- **Toolchain methodology changed mid-project** (Gmsh/CalculiX 2.21 →
  Netgen/CalculiX 2.23), now locked for Stage 6 onward. Any future
  return to Gmsh for a similarly tightly-refined curved-feature
  geometry should anticipate the same class of failure documented
  here, unless the underlying Gmsh defect is independently resolved
  upstream.
- The corner-bite void (~102.7mm³) between the R3 fillet and the
  gusset's flat chord face is a permanent, intentional feature of the
  Option A geometry decision — not evaluated for its own local stress
  effect in this stage.
- The "near-toe" stress search used a fixed sampling radius (3mm in
  the convergence table, 8mm in the visualization crop); the toe
  concentration is a line feature along the full 40mm width, so this
  metric reports the maximum found within that window, not a fixed
  material point tracked across levels.
- The Stage 5 global stress peak (hole edge, 286.22 MPa) has not been
  checked against Stage 2's isolated hole-in-cantilever result (135.35
  MPa observed / no reconciled literature Kt); the two loading states
  differ substantially (Stage 2: simple cantilever bending; Stage 5:
  bent-bracket combined load path), so a direct comparison was not
  attempted and should not be assumed valid.

## Figures (Stage 5)

1. `stage5_full_bracket_vm_iso.png` — full bracket, isometric, von
   Mises, global 0–286.22 MPa scale
2. `stage5_full_bracket_vm_top.png` — plan view showing Flange A, hole,
   and the transition to Flange B in one frame
3. `stage5_full_bracket_vm_front.png`, `stage5_full_bracket_vm_back.png`
   — orthogonal elevation views
4. `stage5_hole_vm_top.png`, `stage5_hole_vm_angle.png` — hole
   close-up, showing the two-lobe Kirsch-type concentration pattern
   (governing global peak, 286.22 MPa)
5. `stage5_toe_vm_top.png`, `stage5_toe_vm_angle.png` — toe fillet
   close-up; the top-down view confirms the concentration is a
   continuous ridge along the full width, not a localized artifact

## Reproducibility

- Mesh generation: `scripts/mesh_netgen_level1_v2.py [1|2|3]`, run
  with system Python (`/usr/bin/python3` — must not be run inside the
  `ccx223` conda environment, which does not have `netgen` installed)
- Solve: `conda activate ccx223`, then
  `ccx -i simulation/stage5_netgen_L[N]_solve`
- Result extraction: `scripts/extract_convergence_results.py`
  (`/usr/bin/python3`)
- Stress-location diagnostics: `scripts/diagnose_peak_stress.py`
- Visualization: `scripts/visualize_stage5_stress.py`
  (`/usr/bin/python3`; requires PyVista, offscreen rendering)

## Status Log

- [x] Specification defined and approved
- [x] Baseline geometry (OpenSCAD reference + Gmsh OCC verification mesh)
- [x] Baseline mesh (tet, hand-verified; hex, automated + validated)
- [x] Baseline CalculiX model
- [x] Baseline solve + sanity check
- [x] Mesh convergence study (5 levels, hex, converged)
- [x] Analytical verification (σₓₓ primary comparison, root-cause
      investigation of boundary effects complete)
- [x] Progressive geometry: holes
  - [x] Geometry and mesh methodology (Gmsh OCC boolean, unstructured
        tet, local refinement)
  - [x] Mesh convergence study (2 levels, 0.28% change, converged)
  - [x] Literature investigation (3 rounds; qualitative mechanism
        confirmed, quantitative gap unresolved and documented as such)
  - [x] Stage conclusion documented (established vs. not established)
- [x] Progressive geometry: shoulder fillet
  - [x] Specification gate (Peterson chart applicability and range
        verified before geometry construction)
  - [x] Geometry and mesh methodology (Gmsh OCC fillet, sigmoid
        Distance+Threshold refinement, two defects diagnosed and fixed)
  - [x] Mesh convergence study (3 solved levels; bottom surface
        converged, top surface accepted with documented 2.3% residual)
  - [x] Baseline sanity check (8/8 closed, including equilibrium
        verification to ≈0.0001% error)
  - [x] Stage conclusion documented (established vs. not established,
        open convergence item preserved)
- [x] Progressive geometry: right-angle frame (4a, sharp corner)
  - [x] Specification gate (scope change to introduce the frame here
        rather than at the originally-planned full-bracket stage,
        discussed and approved)
  - [x] Closed-form statics + Castigliano reference locked before FEA
  - [x] Sign-error caught via FEA comparison, root-caused via two
        independent physical arguments, corrected and documented
  - [x] Mesh convergence study (3 levels, all mid-span points <2%,
        decisively converged)
  - [x] Corner singularity confirmed non-convergent by direct evidence
  - [x] Stage conclusion documented (established vs. not established)
- [x] Progressive geometry: right-angle frame with gusset (4b)
  - [x] Geometry built on the verified 4a frame; old singularity
        confirmed eliminated by direct inspection
  - [x] Global stiffening effect verified (41.6% deflection reduction,
        near-exact equilibrium)
  - [x] Valid mid-span comparisons confirmed unchanged from 4a
  - [x] New gusset-toe stress concentration identified and
        characterized via targeted 3-level convergence study
        (non-convergent, spatially confined)
  - [x] Stage conclusion documented (established vs. not established,
        taper/fillet redesign explicitly deferred, not performed)
- [x] Full L-bracket integration (Stage 5)
  - [x] Final geometry locked (R3 fillet + quadrilateral gusset + R1
        toe fillet + Ø5mm hole), verified clean via BRep diagnostics
  - [x] Gmsh C3D10 Jacobian failure investigated and characterized
        (16 mechanisms tested and ruled out); Netgen-based
        straight-sided quadratic promotion adopted as the resolution
  - [x] Toolchain methodology change (Netgen + CalculiX 2.23) locked
        for the remainder of the project, documented explicitly
  - [x] 3-level mesh convergence study (displacement and toe stress
        both converged; equilibrium ≤0.0008% at all levels)
  - [x] Global stress field inspected directly; hole edge identified
        as the true governing peak (286.22 MPa), toe fillet confirmed
        as a separate, now-convergent local concentration (230.55 MPa)
  - [x] R1 toe fillet confirmed to resolve Stage 4b's non-convergent
        toe singularity
  - [x] PyVista von Mises visualization (multiple angles, full bracket
        + hole + toe close-ups)
  - [x] Stage conclusion documented (established vs. not established)
  - [ ] Git commit and push
- [ ] Parametric mass-minimization study
- [ ] Final documentation and figures
