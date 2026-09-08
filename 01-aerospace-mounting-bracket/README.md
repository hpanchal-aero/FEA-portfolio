# Project 01 — Ribbed L-Bracket: Mass-Minimized Aerospace Mounting Bracket

**Status:** Verification plate stage complete. Progressive geometry —
holes and shoulder fillet stages complete (both verified, with
documented open items). Gusset, full L-bracket not yet started.

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
   **shoulder fillet (complete)** → gusset
3. Full L-bracket baseline model and mesh convergence
4. Parametric mass-minimization study (the actual research question)

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
- [ ] Progressive geometry: gusset
- [ ] Full L-bracket baseline
- [ ] Parametric mass-minimization study
- [ ] Final documentation and figures
