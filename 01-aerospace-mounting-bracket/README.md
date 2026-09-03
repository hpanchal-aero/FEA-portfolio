# Project 01 — Ribbed L-Bracket: Mass-Minimized Aerospace Mounting Bracket

**Status:** Verification plate stage complete. L-bracket geometry (holes,
fillet, gusset) not yet started.

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

1. **Verification plate** (this stage) — plain flat cantilever plate,
   analytically verifiable via Euler-Bernoulli beam theory. Purpose:
   prove the OpenSCAD → Gmsh → CalculiX → Python pipeline is
   implemented correctly before trusting it on a geometry with no
   independent check.
2. Progressive geometry introduction — holes → fillet → gusset
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

## Figures

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

## Engineering Interpretation

The 3D solid FEA model reproduces Euler-Bernoulli beam theory to
within ~1.4% across the mid-span, and the root-region deviation is
fully explained by a known, physically expected 3D boundary effect
(Saint-Venant decay from the rigid encastre constraint), not by
modeling error. Two independent mesh discretizations agree closely.
The implementation is considered verified and the pipeline trustworthy
for the next stage.

## Limitations

- No physical validation performed (no test data available for this
  specific geometry)
- The 12g design acceleration is a preliminary assumed value, not a
  sourced requirement
- The structured hex meshing methodology is not general-purpose and
  will not extend to the filleted/holed L-bracket geometry
- L/h = 15 slenderness ratio is "reasonable but approximate" for
  strict Euler-Bernoulli applicability

## Status Log

- [x] Specification defined and approved
- [x] Baseline geometry (OpenSCAD reference + Gmsh OCC verification mesh)
- [x] Baseline mesh (tet, hand-verified; hex, automated + validated)
- [x] Baseline CalculiX model
- [x] Baseline solve + sanity check
- [x] Mesh convergence study (5 levels, hex, converged)
- [x] Analytical verification (σₓₓ primary comparison, root-cause
      investigation of boundary effects complete)
- [ ] Progressive geometry: holes
- [ ] Progressive geometry: fillet
- [ ] Progressive geometry: gusset
- [ ] Full L-bracket baseline
- [ ] Parametric mass-minimization study
- [ ] Final documentation and figures
