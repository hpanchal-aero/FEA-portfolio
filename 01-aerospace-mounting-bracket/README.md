# Project 01 — Ribbed L-Bracket: Mass-Minimized Aerospace Mounting Bracket

**Status:** Verification plate stage complete. Progressive geometry —
holes stage complete (verified, literature-investigated, unresolved
magnitude gap documented). Fillet, gusset, full L-bracket not yet
started.

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
2. Progressive geometry introduction — **holes (complete)** → fillet
   → gusset
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
- [ ] Progressive geometry: fillet
- [ ] Progressive geometry: gusset
- [ ] Full L-bracket baseline
- [ ] Parametric mass-minimization study
- [ ] Final documentation and figures
