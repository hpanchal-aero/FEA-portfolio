# Project 03 — UAV Wing Structural Analysis — Development Log

## Specification (locked)

**Research question**: How does rib configuration (spacing, then
non-uniform placement at a fixed count) affect a tapered UAV wing's
tip deflection, skin buckling resistance, and structural margin under
a combined self-weight + elliptical aerodynamic load, built up from a
verified spar-to-skin structural model with independently selected
skin thickness?

**Geometry**: NACA 0012 airfoil, tapered planform, root chord 300 mm,
tip chord 150 mm, 1 m half-span.

**Material**: Al 2024-T3, chosen as a standard isotropic aerospace-
skin alloy — deliberately distinct from Al 6061-T6 used in Projects
01/02, without taking on composite-laminate modeling complexity.
Property values (E, nu, rho, yield) to be sourced and cited before
being locked into the first `.inp` file.

**Loading**:
- Combined self-weight (gravity) + elliptical spanwise aerodynamic
  load, applied together as ONE superposed static case (not screened
  as two separate cases) — linear elastic, small-deflection theory,
  so superposition is valid. A candidate configuration fails the
  structural check if this combined stress exceeds yield/1.5, i.e. the
  self-weight case is not judged in isolation; the question is whether
  there is enough margin left for the aero load once self-weight is
  already accounted for.
- A linear eigenvalue buckling analysis (`*BUCKLE`) on the skin under
  the same combined load — this is a NEW analysis type for this
  portfolio (Projects 01/02 were linear static only) and requires its
  own closed-form-verified sub-stage (a simple flat-plate buckling
  case checked against a closed-form critical buckling load) before
  being trusted as a screening metric across candidate designs.

**Constraints**: von Mises stress <= yield/1.5 (FoS convention
consistent with Project 02). No separate displacement limit set at
this stage; tip deflection is tracked as a comparison metric, not
(yet) a hard constraint.

**Explicitly out of scope**: composite/laminate material modeling;
CFD-derived pressure loading (the elliptical distribution is an
analytical idealization, not from Project 04's future CFD work, though
this project is intended as Project 04's eventual structural
foundation per the original portfolio sequencing); aileron/control-
surface modeling; torsional/aileron-reversal load cases.

**Build sequence** (locked, six stages):
1. Spar alone, verified against Euler-Bernoulli beam theory (own I).
2. + ribs at a placeholder spacing (structure only, not the rib
   study — establishes the assembly before skin is added).
3. + skin (full NACA 0012 airfoil, skin structurally participating,
   NOT simplified to a non-load-bearing shell, per explicit decision).
   Skin thickness selection study runs here: multiple candidate
   thicknesses (count and values TBD when this stage starts, current
   proposal ~3 candidates spanning roughly 0.5-1.5mm, to be confirmed),
   scored by a combined metric (tip deflection, buckling load factor,
   mass) under the real combined self-weight+aero load. Winning
   thickness locked for every subsequent stage.
4. Load application refinement: distributing the elliptical spanwise
   load properly over the now-complete skin surface (Stages 1-2 will
   necessarily apply the load in a simplified/direct-to-spar way,
   since no skin exists yet to distribute it over -- this stage
   revisits that, analogous to how Project 02's load mechanism changed
   between Configuration A and B).
5. Rib spacing study (Phase A): three uniform spacings (100mm, 50mm,
   20mm), spar position/dimensions and the now-locked skin thickness
   held fixed. Each spacing determines an actual rib count once
   applied to the 1m half-span (exact counts TBD once end-condition
   details, i.e. whether ribs sit exactly at root/tip or are offset,
   are settled at that stage). Scored by the same combined metric.
   Winning spacing's rib count N is locked.
6. Rib placement study (Phase B): root and tip ribs held fixed as
   anchor points; the N-2 interior ribs' positions are varied across
   many generated candidate placements (generation method TBD at that
   stage -- likely a random/parametric search, mirroring Project 02's
   random pocket-shape approach). All candidates scored by the
   combined metric under the combined load; the top 25 are then
   checked individually for the self-weight+aero combined-stress
   margin criterion above. The best surviving design is reported as
   the project's primary result.

**Decision log** (major decisions made during specification, with
reasoning, for future reference):
- Distributed elliptical load chosen over a point/tip load: a point
  load only tests spar bending stiffness and would not meaningfully
  discriminate between rib/skin configurations, unlike a distributed
  load -- directly informed by Project 02's finding that load-case
  choice determines whether a comparative study has real
  discriminating content (Config A's axial-only case reduced to a
  pure mass-vs-area trade for prismatic geometry).
- Real NACA 0012 airfoil chosen over a simplified rectangular wing-box,
  with skin structurally participating (not idealized as non-load-
  bearing) -- explicitly chosen despite the added geometric and
  verification complexity, because the user judged that simplifying
  would invalidate the results/values for this specific study's
  purpose.
- Staged build-up (spar -> +ribs -> +skin) chosen specifically so each
  structural element's individual contribution to stiffness/strength
  can be isolated and quantified, not just to manage implementation
  complexity -- this quantification is itself a stated goal of the
  research question, not merely a verification convenience.
- Spar position/dimensions held fixed for the ENTIRE project (not
  varied or optimized) -- a deliberate scope-limiting decision, after
  discussing the alternative (a value-selection step at every stage)
  and judging that optimizing spar geometry before ribs/skin exist
  could bias results toward a baseline that isn't representative of
  the final structure. A follow-up spar-position/dimension study on
  the completed structure is noted as a possible future project rather
  than in-scope here.
- Skin thickness selection is a study (multiple candidates scored),
  not a single assumed value -- distinguishing it from spar
  position/dimensions, which remain a fixed, documented assumption.
- Combined-metric scoring approach (tip deflection, buckling load
  factor, mass) carried over from Project 02's combined-score
  methodology, adapted with buckling load factor replacing net-section
  stress as the "structural margin" term, since skin buckling
  resistance is the more physically relevant failure mode for a thin
  wing skin under bending-induced compression, per the user's explicit
  request that ribs' role in preventing skin buckling be captured.
- Self-weight is NOT screened in isolation with an arbitrary threshold
  -- initial framing considered a standalone self-weight deflection/
  stress limit, but was revised to a physically motivated criterion:
  self-weight and aero load are superposed, and a design fails only if
  the combined stress leaves no margin before the aero load is even
  applied on top of self-weight sag. This was an explicit correction
  during specification, recorded so the reasoning is not lost.

**Toolchain**: Gmsh OCC kernel (geometry), Netgen (meshing,
straight-sided C3D10 via manual midside-node promotion), CalculiX 2.23
(solver, `ccx223` conda environment), Python (NumPy, PyVista,
Matplotlib) -- unchanged from Projects 01/02.

**Directory structure**:03-uav-wing-structure/
├── README.md
├── development-log.md
├── scripts/
├── mesh/
├── simulation/
├── results/
└── figures/
└── contours/

## Stage 1: spar alone

### Geometry (scripts/build_stage1_spar.py)

Tapered I-beam, strong axis vertical (matching real wing-spar
orientation: flanges resist bending from lift/gravity, both acting
perpendicular to the chord plane). Root (z=0) height 28.8 mm = 80% of
NACA0012 local thickness at 30% chord for a 300mm root chord (thickness
formula: y_t/c = 5*0.12*(0.2969*sqrt(x/c) - 0.1260*(x/c) -
0.3516*(x/c)^2 + 0.2843*(x/c)^3 - 0.1015*(x/c)^4), giving 0.11997*chord
at x/c=0.30). Tip (z=1000mm) height 14.4mm, same formula at 150mm tip
chord. Web 2mm, caps 8mm x 3mm, both constant along span (only overall
height H(z) tapers linearly). Built via Gmsh OCC addThruSections
(loft) between root and tip I-beam profile wires.

Validated: root/tip cross-sectional areas matched an independent
analytic I-beam-area formula exactly (93.6000 / 64.8000 mm^2, shoelace
vs. analytic agreement to machine precision); total volume matched a
closed-form linear-taper expression (SPAN*(A_root+A_tip)/2 = 79,200
mm^3, exact -- this closed form is valid specifically because only H
varies along the taper, with web/cap dimensions held constant, making
cross-sectional area linear in span position for this profile; it is
not a generic loft-volume approximation). Bounding box exact, 14
surfaces, 0 slivers, cleanliness confirmed (14 model surfaces = 14
solid-boundary surfaces, no stray faces).

Mesh feasibility: 3.0/2.0/1.5mm all meshed cleanly (24,879 / 60,676 /
117,143 tets), no thin-wall failures despite the 2mm web -- unlike
Project 02's flat panels, this is a closed I-beam profile, not an open
thin plate, so it behaved better at fine mesh sizes.

### Load (scripts/mesh_stage1_spar.py)

PLACEHOLDER total lift: 24.5 N, from an assumed 5 kg total aircraft
mass at load factor n=1 (straight-and-level flight), i.e. total lift =
weight = 5*9.81 = 49.05N over the full wing, 24.5N per half-span
model. Explicitly flagged as not a sourced design requirement; to be
revisited in Stage 4.

Elliptical spanwise distribution: q(z) = q0*sqrt(1-(z/b)^2), q0 solved
so integral over [0,b] equals the total lift. Applied as nodal forces
on the spar's top-cap OUTER FACE (the full 8mm-wide x 1000mm-long
face, not a thin edge line as first intended -- an early implementation
mistake was caught and, on reflection, kept: distributing over a face
is a more physically reasonable Stage-1 approximation for where an
eventual skin's aero pressure would act, given no skin exists yet).
Nodal weights are proportional to local q(z) at each node's span
position, normalized so the total sums to exactly the target lift.

Self-weight: `*DLOAD, GRAV` at 9810 mm/s^2 (i.e. 9.81 m/s^2 in the
mm-N-s-tonne-MPa consistent unit system), material density from
`*DENSITY`.

Boundary condition: root face (z=0) nodes fully fixed (cantilever).

### MAJOR BUG: density unit error (1000x), found and fixed

**Symptom**: first solve gave equilibrium error -199.97%, tip
deflection 454mm (45% of the 1000mm span -- physically absurd for
small-deflection linear elastic theory), and peak von Mises 1771 MPa
(over 6x the material's yield). The root-reaction sum showed
self-weight of 2159.9 N, i.e. an implied spar mass of ~220 kg.

**Diagnosis**: CalculiX's consistent unit system (mm, N, s, MPa)
requires density in tonne/mm^3, not kg/mm^3. The correct conversion is
density(kg/m^3) * 1e-12 = density(tonne/mm^3). The `RHO` constant used
throughout this project's scripts (`2780e-9`) was off by exactly
1000x -- it should have been `2780e-12`. Confirmed independently three
ways: (1) a fully SI-unit calculation (no portfolio unit system
involved) gave spar mass = 220.176 g, weight = 2.16 N; (2) CalculiX's
own `*DLOAD GRAV` documentation confirms it expects a physical
acceleration magnitude and computes mass loading directly from
density x volume x acceleration, with no separate scaling; (3) the
measured root reaction (2134.8 N) and the formula-predicted self-weight
(RHO*A*G integrated over span, 2159.9N) agreed with each other almost
exactly, both being ~1000x the correct physical value (2.16N),
confirming the error was in the shared input (RHO), not the formula
structure or CalculiX's handling of it.

**This same 1000x error exists in Project 01 and Project 02's scripts**
(e.g. Project 02's `RHO = 2700e-9` for Al 6061-T6), but never affected
any actual FEA result there, because neither project ever invoked
CalculiX's mass-dependent physics (`*DLOAD GRAV`, modal analysis,
dynamics) -- density in those projects only ever fed a Python-side
print statement (`volume * RHO * 1000`, labeled "mass in grams"). That
print statement's extra `*1000` happened to exactly cancel the
density's 1000x error, so printed masses (e.g. Project 02's Config A:
257.31 g) were numerically correct by coincidental cancellation, not
because RHO was actually right. No re-run or correction of Project
01/02 results is needed, since density was never used for anything
except that display-only calculation there -- this note is recorded
here as a caution for any future work in this portfolio that DOES use
gravity loading, modal analysis, or dynamics with an `RHO` constant
copied from those projects' pattern.

**Fix**: `RHO` corrected from `2780e-9` to `2780e-12` in
`build_stage1_spar.py`, `mesh_stage1_spar.py`, and
`extract_stage1_spar.py`. The build script's mass-print formula was
also corrected from `volume*RHO*1000` (labeled grams) to
`volume*RHO*1e6` (tonnes to grams is x1e6, not x1000), now giving the
correct 220.18 g for the right reason rather than by cancellation. A
cosmetic display-only follow-up bug (the same print line's density
echo, `RHO*1e9` re-labeled as kg/m^3, no longer correct after the RHO
fix) was also caught and corrected to `RHO*1e12`.

Also found during this same debugging pass: the extraction script's
equilibrium check had an unrelated sign error (compared the root
reaction directly against the net applied force, rather than against
its equal-and-opposite value) -- this made the (now-correct) physics
look like a "-199.86%" equilibrium error even after the density fix,
which was actually just a factor-of-(-2) discrepancy from the sign
mistake, unrelated to density. Fixed by explicitly computing
`expected_reaction_Y = -(applied_lift - self_weight)`.

### Verification results (post-fix)

Reference: numerical tapered Euler-Bernoulli beam theory
(scripts/extract_stage1_spar.py, `tapered_beam_tip_deflection()`),
integrating the true tapered I(z) under the actual elliptical lift +
self-weight distribution via double numerical integration (moment,
then curvature, then slope, then deflection, on a 20,000-point grid) --
no closed form exists for this taper + combined load, so this is a
numerical, not closed-form, reference, distinguished explicitly from
the sourced closed-form checks used in Projects 01/02.

A supplementary Tier-1 sanity bound (prismatic root-section
approximation, uniformly-distributed-load simplification, no
self-weight) was computed as a quick order-of-magnitude gut check:
4.19mm, correctly LARGER (stiffer, i.e. less deflection... wait, this
bound predicts less deflection, meaning it should be smaller than the
true tapered result if it is "conservative/stiffer" -- confirmed: 4.19
Tier-1 vs 3.84mm actual is NOT smaller, it is larger. This means Tier 1
actually predicts MORE deflection despite using the stiffer root
section throughout, which at first reading seems contradictory to
"stiffer section = less deflection" -- the resolution is that Tier 1's
UDL approximation OVER-estimates the effective loading relative to the
true elliptical distribution's actual moment arm effect combined with
omitting the (deflection-reducing... no, self-weight ADDS deflection)
self weight. Given Tier 1 omits self-weight (which adds deflection) yet
still predicts MORE deflection than the full tapered+self-weight FEA
result, the two simplifications (stiffer prismatic root section,
UDL-vs-elliptical load shape) do not cleanly bound the answer in one
direction -- Tier 1 is confirmed here to be a rough order-of-magnitude
check only, not a validated bound, and its relationship to the true
answer should not be over-interpreted beyond "same order of magnitude."

| | 2.0mm | 1.5mm |
|---|---|---|
| Tets / C3D10 nodes / equations | 60,676 / 118,841 / 356,523 | 117,143 / 228,813 / 686,439 |
| Element / equation % of ceiling | 26.4% / 30.2% | 50.9% / 58.2% |
| Peak RSS / wall time | 1.21 GB / 13.2s | 2.59 GB / 32.3s |
| Root reaction vs expected | -0.1363% | -0.0860% |
| Tip Uy vs tapered theory | -0.1625% (3.841mm vs 3.847mm) | -0.2084% (3.839mm vs 3.847mm) |
| Peak von Mises | 16.60 MPa (root flange-web corner) | 18.24 MPa (same location) |

Tip deflection and equilibrium are converged between mesh levels (both
shifts <0.05 percentage points). Peak von Mises shows a real +9.9%
mesh sensitivity, consistent with the discontinuity-driven,
slower-converging peak-stress pattern already documented in Project
02 (e.g. Configuration A's base-corner peaks) -- not treated as a
design value at this stage, and not a cause for concern given both
values sit far below yield.

## Status Log

- Specification defined and approved
- Stage 1 (spar alone): geometry validated, meshed (2 levels), solved,
  verified against numerical tapered beam theory (-0.16 to -0.21%
  agreement), one major bug (density units, 1000x) and one minor bug
  (equilibrium sign error) found and fixed
- Stage 2 (+ ribs, placeholder spacing): not started
- Stage 3 (+ skin, skin thickness study): not started
- Stage 4 (load application refinement): not started
- Stage 5 / Phase A (rib spacing): not started
- Phase B (rib placement search): not started
- Documentation: this entry
- Git commit and push: pending (this stage)
