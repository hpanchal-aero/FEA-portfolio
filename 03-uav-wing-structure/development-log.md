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

A supplementary Tier-1 sanity bound was also computed: prismatic root
cross-section, a uniformly distributed load in place of the elliptical
distribution, and no self-weight. It gives 4.19 mm versus 3.84 mm from
the FEA. It is NOT a bound. Decomposing it during the Stage 2
documentation pass (same beam integration, one simplification at a
time, starting from the true model at 3.847 mm) gives:

- omitting self-weight: 3.847 -> 4.306 mm (+11.9%), because self-weight
  acts opposite to the lift and reduces the net moment;
- a uniform load of the same total in place of the elliptical
  distribution: 4.306 -> 5.723 mm (+32.9%, tapered I, no self-weight),
  because it carries more moment about the root (12,250 vs 10,398 N*mm);
- holding the root cross-section constant along the span:
  4.306 -> 3.242 mm (-24.7%, elliptical, no self-weight), because the
  root section is stiffer than the tapered section.

Two simplifications raise the prediction and one lowers it, so the
+8.9% net is partly cancellation of errors. Tier 1 is a
same-order-of-magnitude check only. An earlier version of this section
stated that self-weight adds deflection, and the note printed by
scripts/extract_stage1_spar.py said Tier 1 should underestimate the FEA
result; both statements were wrong and are corrected here (the script's
printed note is corrected in the same commit).

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

## Addenda to Stage 1 (added during Stage 2)

### Addendum A: the Stage 1 root-reaction residual is a reporting convention

Stage 1 recorded equilibrium errors of -0.1363% (2.0 mm) and -0.0860%
(1.5 mm). Stage 2 showed a mesh-independent residual of about 0.24%, and
isolating it (scripts/equilibrium_split_test.py, then
scripts/check_equilibrium_constrained_loads.py) showed that the
reaction CalculiX reports at constrained nodes (FORC) does not include
external load applied directly on those nodes. Adding that load back
(gravity share via the consistent nodal load of a straight-sided
C3D10, -V/20 per corner node and +V/5 per midside node, plus any CLOAD
entries on the fixed nodes) collapses the residual to 1e-5 N or less
on all four runs:

| Run | Load on fixed nodes | Naive error | Corrected error |
|---|---|---|---|
| Stage 1, 2.0 mm | lift CLOAD +0.031165 N, gravity -0.000681 N | +0.1363% | -0.0001% |
| Stage 1, 1.5 mm | lift CLOAD +0.019773 N, gravity -0.000552 N | +0.0860% | -0.0000% |
| Stage 2, 2.0 mm | gravity -0.050287 N | -0.2363% | -0.0000% |
| Stage 2, 2.5 mm | gravity -0.050370 N | -0.2366% | +0.0000% |

Sign convention: the original Stage 1 extractor divided by the signed
(negative) expected reaction, so its recorded -0.1363% / -0.0860%
correspond to +0.1363% / +0.0860% under the |expected| convention used in
the table above and in the Stage 2 extractor. In magnitude terms the
measured Stage 1 reaction was smaller than expected (lift load applied
directly on the fixed nodes is missing from the reported reaction) and
the Stage 2 reaction was larger (the gravity share on the fixed nodes is
missing). The CalculiX reporting behaviour was
established empirically on these four runs, not checked against the
CalculiX documentation. Loads on the fixed face have essentially zero
moment arm about the root, so tip deflection is unaffected. No Stage 1
result changes.

### Addendum B: Stage 1 stress verification

Stage 1 verified deflection and equilibrium but never stress. Added
with scripts/verify_stage2_stress.py (spar-only mode): axial stress
sigma_zz on the outer flange faces at z = 100, 300 and 500 mm (nodes
with |x| < 3.5 mm, within 0.15 mm of the flange surface) against
Euler-Bernoulli bending stress -M(z)*y/I(z), with M(z) from the same
beam integration as the deflection reference. Criterion set in
advance: mean deviation within 5%.

| Mesh | Largest mean deviation over 6 station/side cases |
|---|---|
| 2.0 mm | 0.14% |
| 1.5 mm | 0.22% |

Signs are correct (top flange in compression, bottom in tension). The
root flange-web corner peak (16.60 / 18.24 MPa) remains unconverged and
is not covered by this check.

## Stage 2: spar + ribs (200 mm placeholder spacing, 6 ribs)

Purpose: add ribs to the verified spar and quantify their contribution
before skin is introduced. This is NOT the rib-spacing study (Stage 5).

### Decisions (locked)

- Ribs are true NACA 0012 airfoil-shaped thin plates, scaled to the
  local chord (300 mm root to 150 mm tip), not simplified rectangles.
- Stations z = 0, 200, 400, 600, 800, 1000 mm (root + tip + 4 interior).
- Rib thickness 1.5 mm, centred on the station, so the root rib spans
  z = -0.75 to +0.75 mm and the tip rib 999.25 to 1000.75 mm.
- Ribs fused to the spar by boolean union.
- Same placeholder load as Stage 1 (self-weight + 24.5 N elliptical lift
  on the spar top-cap outer face).
- Fixed root = nodes at the assembly's true minimum z (-0.75 mm), computed
  from the mesh, because z = 0 is now inside solid material (confirmed).
- Mesh sizes 2.0 mm (primary) and 2.5 mm (convergence check).

### Geometry (scripts/build_stage2_ribs.py)

Each rib outline is 40 cosine-spaced points per surface, extruded 1.5 mm.
Per-rib area x thickness matches the extruded volume exactly at all six
stations. Total rib volume 39,345.999 mm^3; fused assembly volume
117,951.9990 mm^3; mass 327.91 g (Stage 1 spar alone: 220.18 g, +107.73 g,
+49%). Final model: 1 volume, 540 surfaces, all 540 on the solid boundary,
0 slivers. A PyVista render of the fused assembly shows six correctly
tapered ribs correctly placed along the spar.

Overlap cross-check (done with a read-only diagnostic, not inside the
build script): boolean intersection of the spar with each rib gives
70.1919, 131.7600, 123.1200, 114.4800, 105.8400, 48.6081 mm^3, total
594.0000 mm^3, and 79,200.0000 + 39,345.9990 - 594.0000 = 117,951.9990,
which reproduces the fused volume exactly.

Three geometry bugs were found and fixed:

1. **Duplicate leading-edge point.** airfoil_outline() returned upper +
   lower, where lower[-1] and upper[0] were the same point differing only
   by 0.0 vs -0.0, so the closing edge had zero length and Gmsh raised
   "Could not create line". Fixed with upper + lower[:-1].
2. **Fuse producing 2 volumes.** A single 6-way fuse returned 2 volumes.
   A read-only diagnostic showed every rib overlaps the spar (no
   positioning gap), and fusing one rib at a time gave 1 volume at every
   step. Even a sequential loop run after building all six ribs still left
   an orphaned volume at model level (every fuse reported 1 while the model
   held 2). Restructuring to build-one-fuse-one, fusing each rib
   immediately after creating it, fixed it. The underlying cause was not
   isolated; the working hypothesis is OCC boolean fragility when ribs sit
   flush with the spar's exact end faces (the run logged "Cannot bind
   existing OpenCASCADE volume" warnings).
3. **Bounding-box check.** The z-extent check failed by a 5e-6 mm residual
   (expected [-0.75, 1000.75]) accumulated over five sequential fuses. It
   now uses a separate BBOX_FUSE_TOL of 1e-3 mm; Stage 1's single-extrusion
   checks keep 1e-6 mm.

### Mesh feasibility (mesh-only, Netgen)

| maxh | Tets | Linear nodes | Time | Tets vs 230,000 ceiling |
|---|---|---|---|---|
| 3.0 mm | 54,086 | 19,870 | 35.9 s | 23.5% |
| 2.0 mm | 108,415 | 36,661 | 69.1 s | 47.1% |
| 1.5 mm | 196,862 | 65,645 | 141.5 s | 85.6% |

No intersection warnings at any size. 1.5 mm was rejected: the estimated
equation count was about 97.6% of the empirical 1.18M-equation ceiling
(the level at which a Project 01 job crashed the WSL host), so the
convergence pair is 2.0 / 2.5 mm.

### Load and boundary conditions (scripts/mesh_stage2_ribs.py)

Meshed volume matches the geometry volume to -0.0000% at both sizes. The
fixed set is the 9,410 (2.0 mm) / 7,147 (2.5 mm) nodes at z = -0.75 mm.
Load nodes are the spar top-cap outer-face nodes: |y - H(z)/2| < 0.15 mm
and |x| < 4.5 mm, 0 <= z <= 1000 mm. The elliptical lift is assigned slab
by slab (5 mm slabs in z): each slab receives the exact analytic integral
of q(z) over that slab, shared equally by the load nodes inside it.

### Debugging history of the load and the tip metric

The first Stage 2 solves looked plausible but did not converge, and the
cause took four steps to find. Each step is recorded because the
intermediate numbers were misleading.

1. **Rib-surface nodes in the load selection.** The y-only top-face test
   also caught rib-surface nodes (706 of 10,113 at 2.0 mm, 499 of 6,136 at
   2.5 mm, with |x| outside the spar cap), because the ribs are taller than
   the spar. Restricting |x| fixed the contamination.
2. **Per-node load weighting versus mesh refinement.** Weighting each node
   by q(z) at its own z (not by tributary length) over-weights regions the
   mesh refines locally around the 1.5 mm ribs. The applied load moment
   was -4.29% (2.0 mm) and -0.35% (2.5 mm) against the intended
   10,398.1 N*mm, and 10 load nodes sat on the fixed face, where their
   force goes straight into the support. Stage 1's control meshes were
   within -0.09% / -0.10%. Slab assignment brought both Stage 2 meshes to
   10,398.0 and 10,398.2 N*mm.
3. **Tip metric.** The tip deflection was first taken as the mean over all
   nodes within 2 mm of the tip. Now that the tip rib dominates that band
   (7,969 nodes, spread across the chord), the mean depends on how nodes
   fall across the chord. Uy is linear in x across the tip (twist slope
   about -0.00165 mm/mm, roughly 0.096 degrees, residual 0.003 mm),
   caused by the off-axis rib weight acting on the low-torsion open
   I-section. The primary metric is now the spar-axis (|x| < 4.5 mm)
   deflection. The twist magnitude has not been checked against torsion
   theory.
4. **Equilibrium residual.** See Addendum A: a reporting convention, not a
   load error.

Retractions. Three statements made during this stage were wrong:

- The extractor's docstring said ribs, being thin plates normal to the
  span, should not meaningfully change tip deflection. That was an
  expectation stated without computing it.
- When the first (contaminated) results showed tip deflection 10-12%
  below Stage 1, I attributed it to ribs stiffening the root because they
  are taller than the spar. That explanation was never verified, and the
  effect changed with mesh size.
- The Stage 1 vs Stage 2 comparison labelled the change "ribs'
  contribution to stiffness". It conflated added self-weight with
  stiffness.

### Verification results (2.0 mm | 2.5 mm)

| | 2.0 mm | 2.5 mm |
|---|---|---|
| Tets / C3D10 nodes / equations | 108,415 / 216,583 / 649,749 | 67,696 / 141,323 / 423,969 |
| Equation % of ceiling | 55.1% | 35.9% |
| Wall time / peak RSS | 107.6 s / 2.02 GB | 62.6 s / 1.19 GB |
| Equilibrium, fixed-node load included | -0.0000% | +0.0000% |
| Equilibrium, naive | -0.2363% | -0.2366% |
| Spar-axis tip Uy | 3.6791 mm | 3.6795 mm |
| vs ribs-as-weights beam reference (3.6836 mm) | -0.123% | -0.114% |
| Twist slope | -0.001677 mm/mm | -0.001642 mm/mm |
| Flange stress vs beam theory (z = 100/300/500) | largest 0.06% | largest 0.06% |
| Peak von Mises (root corner, unconverged) | 16.96 MPa | 16.16 MPa |

The reference is the Stage 1 tapered-beam integration with the ribs added
as point weights carrying zero bending stiffness (their volumes net of
the spar overlap). The spar-axis deflection agrees between meshes to
0.01% and with the reference to about 0.1%, inside Stage 1's -0.16% /
-0.21% verification band.

### Findings

- The ribs add no measurable spanwise bending stiffness in this load
  case. The whole change from Stage 1 (3.841 -> 3.679 mm, -4.21%) is the
  ribs' added self-weight, which alone predicts -4.25%.
- The ribs cost +107.7 g (+49% mass) and add twist (about 0.1 degrees at
  the tip) without a bending benefit. Their value has to come from the
  skin (buckling resistance, shape retention), which is Stage 3 onward.
- Von Mises in the ribs is near zero everywhere; the spar carries the
  bending stress, decaying from root to tip.

### Figures

figures/stage2_ribs_preview.png (assembly),
figures/contours/stage2_h2p0_uy_iso.png and _uy_side.png (Uy, deformation
x50, the colorbar maximum 3.76 mm is the leading edge of the tip rib
where the twist adds displacement; the spar-axis value is 3.679 mm),
figures/contours/stage2_h2p0_vonmises.png (colour scale capped at the
99th percentile, 11.6 MPa) and _vonmises_root.png (root close-up).
scripts/view_stage2_interactive.py opens a rotatable view with a
deformation slider.

### Limitations and open items

- The load remains the 24.5 N placeholder, applied to the spar top cap.
- Convergence is shown between two levels (2.0 / 2.5 mm, ratio 1.25); a
  third level was excluded on resource grounds.
- The root corner peak stress (16.96 vs 16.16 MPa, 5% apart, at a
  different z on each mesh) is unconverged and is not a design value.
- The twist magnitude is not verified against torsion theory.
- The CalculiX reaction-force convention was established empirically.
- Rib stiffness contribution was measured for this load case only.

## Status Log

- Specification defined and approved
- Stage 1 (spar alone): complete and verified; committed locally
  (1f64068). Bugs: density units (1000x), equilibrium sign. Addenda: the
  root-reaction residual is a reporting convention; stress verified.
- Stage 2 (+ ribs, placeholder spacing): complete and verified
  (documentation in this entry); ribs contribute weight, not bending
  stiffness
- Stage 3 (+ skin, skin thickness study): not started
- Stage 4 (load application refinement): not started
- Stage 5 / Phase A (rib spacing): not started
- Phase B (rib placement search): not started
- Git: Stage 2 commit pending; not pushed
