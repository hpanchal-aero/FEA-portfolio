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

## Stage 3: spar v2 + ribs + skin (skin thickness locked at 1.0 mm)

Status: numerically complete; skin thickness locked at 1.0 mm by the
user. Stage 1/2 describe the old spar and the 24.5 N placeholder load;
this section supersedes the load case and the spar for Stage 3 onward
without altering the earlier records.

### Decisions (locked)

- **Spar redesign (Option C).** The Stage 1/2 spar depth was 80% of the
  local airfoil thickness, leaving the caps 3.6 mm (root) / 1.8 mm (tip)
  inside the airfoil outline, so no skin could reach the spar. Options:
  (A) skin attached at the ribs only, (B) same-material filler strips on
  the caps, (C) spar depth = airfoil depth. The user first chose B on my
  description ("a little mass"); I then retracted that description: the
  filler would add about 120 g and roughly double the spar bending
  inertia (hand estimate x2.1-2.5), i.e. not a thin shim. With the
  corrected numbers the user chose C. Spar v2: web 2.0 mm and caps
  8x3 mm unchanged, cap outer faces on the airfoil envelope; depth
  35.99 mm (root) / 17.99 mm (tip); mass 249.95 g (+29.77 g); spar +
  ribs 357.46 g (Stage 2: 327.9 g).
- **Skin:** closed NACA 0012 skin, Al 2024-T3, **shell elements (S6)**,
  candidates 0.5 / 1.0 / 1.5 mm. Solid tets were rejected: 0.5-1.5 mm
  needs 2-3 elements through the thickness and would exceed the
  ~230,000-element ceiling (a job near that size crashed the WSL host in
  Project 01).
- **Skin convention:** shell nodes on the airfoil envelope, material
  outward (`*SHELL SECTION, OFFSET=-0.5`, outward normals), so the outer
  mold line is the nominal NACA 0012 plus t.
- **Connection:** shared nodes with the solid on the rib rims and the
  spar-cap strips (conformal mesh), no `*TIE` / `*EQUATION`.
- **Load case (replaces the placeholder):** the 24.5 N placeholder was
  incompatible with the skin: at 1.0 mm the skin alone weighs 12.5 N and
  at 1.5 mm structure weight would cancel nearly all of the lift (net
  about 2 N). User: "whatever standard UAV has". Adopted limit load
  factor **n = 3.8** (the cap of 14 CFR 25.337(b): 2.1 + 24000/(W+10000)
  with a ceiling of 3.8; for a 44 lb aircraft the formula gives 4.49 so
  the cap applies; 3.8 is also quoted for normal-category aircraft in
  secondary sources. I could not retrieve the STANAG 4671 / CS-23 text;
  verify before citing). **Aircraft mass 20 kg is my assumption**, chosen
  so the wing pair is about 10 / 16 / 23% of the mass at t = 0.5 / 1.0 /
  1.5 mm. Half-wing lift 3.8 x 20 x 9.81 / 2 = **372.78 N** (elliptical,
  on the spar-cap face as in Stages 1-2), gravity x3.8 on solid and skin.
- **Requirements at limit load:** yield / von Mises >= 1.5 (yield
  **assumed 290 MPa, uncited**) and linear buckling factor >= 1.5.
- **Thickness lock (user):** t = 1.0 mm, the lightest candidate meeting
  both requirements. The 0.75 mm candidate suggested by interpolation
  (factor 1.5 at about 0.65-0.7 mm) was not run.

### Geometry (scripts/build_stage3_solid.py, build_stage3_skin_geom.py)

- Wing envelope: ruled loft of the root and tip NACA 0012 polylines.
  Because the outline scales linearly with chord about the spar axis,
  every facet is planar and each rib is an exact slice. Spar v2 and the
  rib slabs are built oversized and trimmed to the envelope by a
  boolean intersection, so cap tops, rib rims and the skin inner surface
  coincide exactly.
- Validation: volume vs. independent numerical integration (spar
  89,910.4455 mm^3, +0.00001%; spar + ribs 128,582.1216 mm^3), centroid
  z identical to the numerical value, 1 volume, no slivers. **Bug:**
  OCC `getBoundingBox` is loose after the boolean (reported the cutter
  extents); the bounding box is now taken from the nodes of a coarse
  surface mesh.
- Skin imprint: the 78 lateral envelope faces are fragmented with the
  solid into 834 pieces (446 shared with the solid, 388 free); skin area
  459,622.6823 mm^2 equals the independent quad-sum; shared area
  20,011.3512 vs. independent 20,011.3506 mm^2; the solid volume is
  unchanged.
- Known geometry limitations: the polyline trailing-edge closure is
  asymmetric (upper surface ends at the true TE point, the lower
  surface one point earlier, a slanted closing facet; inherited from
  Stage 2, flagged, not changed). The constant-width cap edge (x = -4 mm)
  crosses a scaled facet line near z = 742 mm, producing four wedge
  pieces (about 300 mm^2 in all) whose apex angle (0.18 degrees) bounds
  the skin triangle quality: inherent to the geometry, accepted.

### Skin connection and buckling pipeline (verification sub-stages)

- **`*TIE` is impossible for shells.** The first attempt failed with a
  Python error (`exportStep` does not exist; `WriteStep` does), then
  `gen3dnor` (a copy-paste mismatch: elements declared S6 but written
  with 8 nodes), then "nonpositive jacobian determinant" on every shell
  element even for a single element. Root cause, from the official
  manual and not from guessing: `*TIE` can only be used with 3D
  elements. Forum posts claiming otherwise were not relied on.
- **Shared-node flush test** (100x10x10 mm block, 1 mm shell plate on
  top, tip load 10 N): composite-beam reference 0.041112 mm. The first
  result was 0.021321 mm (+48% too stiff) with no solver error: shell
  elements built on the solid's faces had inconsistent winding, so with
  `OFFSET=-0.5` about half of them put material into the block;
  CalculiX's `*INFO ... opposite normals` messages were the symptom.
  144 of 264 elements had to be flipped; afterwards 0.042998 mm,
  -4.59% vs. the reference (attributed, not proven, to shear
  deformation of this stubby beam). **The code for this test and its
  reference calculation was run as inline code and was never saved as a
  script; it cannot be reproduced from the repository.** In the
  production pipeline the shell winding is set by an explicit outward
  test plus a directed-edge consistency check (0 flips needed, 0
  inconsistent edges).
- **Buckling pipeline verification:** flat plate 200x30x1.0 mm, simply
  supported, closed-form k = 4.0095 (m = 7), sigma_cr = 300.58 MPa:
  S8R meshes 20x4 / 40x8 / 80x16 give factors 0.9801 / 0.9778 / 0.9777
  (-1.99 / -2.22 / -2.23%, converged; the ~2% residual is attributed to
  fully restrained in-plane edges, not proven; accepted by the user as
  verified for screening). A curved root-bay panel (placeholder load,
  idealized edges) gave 2.103 / 2.098 at 8 / 5 mm; it was a pipeline
  stepping stone and is superseded by the full assembly.
  `buckle_verify_extract.py` always printed "deviation from the
  closed-form target of 1.0", meaningless outside the flat-plate case;
  fixed (the target is now an optional argument).

### Mesh and input generation

Pipeline (interpreter split: Netgen / Gmsh-OCC scripts run with
`/usr/bin/python3`, scipy-based builders in the `ccx223` conda env):
Netgen meshes the imprinted solid + skin faces in one conformal mesh
(solid, boundary and shared faces at the solid size, free skin faces at
the skin size; the global maxh is the skin size because per-face sizes
can only refine); free-face nodes carry no shared IDs and are welded by
coordinate (1e-6 mm, the solid node ID survives); the welded skin is
checked (no degenerate or duplicate triangles, edge-manifold, boundary
edges only at the two end planes, area = geometry); winding is set
outward and checked by directed edges; nodes are renumbered contiguously;
one midside-node table serves the C3D10 tets and the S6 triangles; the
`.inp` carries a `** LOADCASE LIFT_N= NG=` header that the extractors
read.

| Mesh (solid / skin) | C3D10 | S6 | Nodes | Equation estimate | Static solve | Buckle solve |
|---|---|---|---|---|---|---|
| 2.5 / 8 mm | 73,209 | 40,080 | 212,645 | 880k (75%) | 73 s, 5.7 GB | — |
| 4 / 8 mm | 44,899 | 34,698 | 150,322 | 660k (56%) | 50 s, 4.35 GB | 125 s, 4.7 GB |
| 4 / 6 mm | 44,868 | 44,528 | 169,941 | 778k (66%) | 72 s, 5.4 GB | 185 s, 5.8 GB |

The 8 mm skin already needs 5.7 GB in the static solve (8 GB watchdog
limit; WSL host crash history), so a finer skin was not attempted.
CalculiX's `.frd` holds only the solid original node IDs plus its own
expanded shell nodes (IDs above the input node count): skin-only
original IDs are absent, so extractors read the `.frd` node coordinates.

### Debugging and retractions (chronological)

- Netgen meshes free faces with their own node copies (probe: 21,756
  free-skin nodes, 0 shared IDs, 3,299 coincident with solid nodes,
  8,303 duplicate groups among free faces): fixed by the coordinate weld.
- First viewer images were identical (camera set after the render) and
  were deleted; fixed with one fresh plotter per image.
- **Reaction forces:** my first hypothesis (the "missing" 4 N sits on
  expanded root nodes) was wrong. `FORC` is nonzero at every node, and
  the sums over the original and over the expanded groups are each zero;
  free nodes carry about 25.5 N of constraint forces for 24.5 N applied
  (the same x2 factor for solid-only loads); `*NODE PRINT, TOTALS=ONLY`
  also failed to close (lift only 25.96 N vs. 24.5 N with spurious X and
  Z components; full model 7.54 N vs. 8.46 N). **No reaction-based
  equilibrium check is interpretable for expanded shells.**
- **Retracted claims:** (1) I called the fully-effective-skin beam a
  "bound"; it is only an approximate reference (vertical-offset ring, no
  Poisson or rib effects). (2) I predicted a limit-load tip deflection of
  about 3.4 mm by scaling the net force; lift scales x15.2 and gravity
  x3.8, and the correct superposition gives the 4.204 mm that was found.
  (3) I described the filler option as adding "a little mass".
  (4) I called the earlier tip agreement with the beam a clean
  confirmation; the deflection profile showed compensating inboard and
  outboard differences. (5) A stale note about an unresolved git
  divergence was wrong: origin/main contained nothing missing locally.
- Script mix-up: the buckling extractor text was saved into
  `make_buckle_inp.py`; fixed by rename and recreate. A Python run with
  the wrong interpreter (no numpy) failed harmlessly.
- Local-size routes for the trailing-edge study: a Netgen mesh-size file
  failed with "Problem in Surface mesh generation" and hung (now run
  under `timeout`); vertex sizes had no effect; **edge sizes on the 34
  edges near the trailing edge at rib 800 work** (`TE_MODE=edge`).

### Verification results

**Load application (t = 1.0 mm, 2.5 mm solid / 8 mm skin).** Tip deflection
of three n = 1 variants vs. the fully-effective-skin beam: lift only
0.33767 vs. 0.33867 mm (-0.29%), + solid gravity 0.28437 vs. 0.28511
(-0.26%), full 0.09193 vs. 0.09287 (-1.01%; absolute offset about 1 um in
each case); effect of solid gravity FEA / beam 0.995, of skin gravity
1.001. Limit load: 4.2040 mm = 15.215 x 0.33767 - 3.8 x 0.24575.

**Tip deflection vs. beam, by thickness (limit load):** 0.5 mm 8.359 vs.
8.292 mm (+0.82%); 1.0 mm 4.204 vs. 4.219 (-0.35%); 1.5 mm 2.523 vs.
2.603 (-3.07%); control t = 0.05 mm 27.603 vs. 27.561 (+0.15%, skin stress
ratio 0.995 at z = 300, 0.993 at z = 700). The deflection profile at 1.0 mm
differs from the beam by +23% at z = 100, +8% at 300, +3.6% at 500 and
-0.35% at the tip (extra inboard slope about 1.5e-4 rad).

**Skin stress / beam, expanded nodes (sigma_zz / -M y / I_tot):** z = 100:
1.006 (all thicknesses roughly 1.0); z = 300: 0.968 / 0.933 / 0.891;
z = 500: - / 0.886 / -; z = 700: 0.879 / 0.739 / 0.578 at t = 0.5 / 1.0 /
1.5 mm. The pattern is unchanged without gravity and the ratio is lowest
on the spar strip; the absolute offset is about 0.6 MPa for t >= 0.5 mm.
**Cause not found** (candidates: shear lag or shear deformation, skin
modeling).

**Mesh convergence.** Solid 2.5 -> 4 mm (t = 1.0, skin 8): tip -0.02%,
skin p99 -0.7%. Skin 8 -> 6 mm: tip +0.002%; buckling factors at
t = 1.0: 3.1388 -> 3.1226 (-0.5%; all six modes -0.5 to -1.1%), at
t = 0.5: 0.9122 -> 0.9325 (+2.2%).

**Buckling at limit load (8 mm skin):** 0.5 mm 0.912 (FAIL); 1.0 mm 3.139;
1.5 mm 8.611 (modes 8.611, 8.786, 12.15, 13.00, 13.95, 14.90). All modes:
upper skin (100%), aft of the spar (x 40-166 mm), out-of-plane, root bays
(0-200 and 200-400 mm), 0% near the edges or rib planes. Plate
plausibility (t = 1.0): buckle width 92.6 mm, static sigma_zz 10.36 MPa x
3.1388 = 32.5 MPa vs. simply supported plate k = 4.02 -> 31.6 MPa
(ratio 1.027).

**Stress.** Skin von Mises (CalculiX expanded nodes; the interpretation is
verified only near the root): p99 36.1 / 18.8 / 11.9 MPa; peak at the
trailing-edge vertex at the ribs 300 / 238 / 169 MPa; maximum at least 3 mm
from the trailing edge 39.4 / 34.2 / 25.3 MPa (assumed-yield factors of
safety 7.4 / 8.5 / 11.5).

**Trailing-edge refinement (t = 1.0, rib 800).** Local median edge near the
vertex 2.12 / 1.23 / 0.54 / 0.28 mm: skin peak 231.8 / 298.7 / 444.1 /
535.6 MPa, solid-only peak 82.9 / 123.8 / 271.2 / 473.3 MPa (no
convergence); rings 1-3 mm 90.3 / 98.4 / 56.3 / 59.8 (still drifting), 3-6
mm 29.3 / 31.7 / 35.4 / 35.7, 6-15 mm 15.1 / 18.1 / 19.1 / 19.3, 15-40 mm
4.9 / 5.3 / 5.7 / 5.9; unrefined ribs moved by at most 2.5%.

All tables are regenerated by `scripts/collect_stage3_results.py`
(`results/stage3_*.csv`, `stage3_summary.json`).

### Findings

- The skin carries essentially all of the bending (skin effectiveness
  fraction 1.00): at limit load a spar alone would deflect 30.5 mm against
  4.2 mm with the 1.0 mm skin.
- Skin buckling governs the thickness: 0.5 mm buckles below limit load
  (0.91), 1.0 mm has a factor of 3.1, 1.5 mm 8.6. The buckle is local to
  the aft panel between spar and trailing edge (width about 93 mm), so the
  200 mm rib spacing is not the controlling length (plate k is nearly
  constant for panels longer than they are wide). Denser ribs should
  help only below about 93 mm spacing (a square-panel estimate gives k of
  about 5.8 at 50 mm and about 24 at 20 mm) — an estimate, not computed.
- With 200 mm spacing the lightest candidate meeting buckling >= 1.5 and
  stress (factor of safety >= 1.5 outside the singularity) is 1.0 mm.
- The factor scales roughly as 1 / load; 0.5 mm would pass for about
  12 kg or less, 1.0 mm up to about 40 kg (approximate).

### Figures

figures/stage3_mesh_overview.png and stage3_mesh_root.png (mesh classes),
stage3_buckling_mode1.png (critical mode, plan view of the upper skin,
the white stripes are gaps between node rows of the point cloud),
stage3_thickness_summary.png, stage3_te_refinement.png (all from
scripts/plot_stage3_summary.py and render_stage3_figures.py).
scripts/view_skin_mesh.py opens a rotatable view (WSLg).

### Limitations and open items

- Unsourced: aircraft mass (assumed), yield (assumed, uncited), the
  primary load-factor source.
- The load is still applied to the spar-cap face (Stage 4).
- Skin rotations at shared nodes are uncoupled from the solid (hinged
  skin); not quantified.
- Linear buckling only; no imperfections or post-buckling.
- Trailing-edge geometry idealized; the singular peak is excluded and the
  real trailing-edge stress is unknown; refinement at one rib only (other
  ribs inferred).
- Wedge-piece sliver triangles (0.18 degrees) accepted; their effect was
  not isolated (a removal test was considered, not run).
- Unexplained skin stress / curvature deficit relative to the beam.
- No reaction-based equilibrium for the shell model.
- Skin-size convergence is two levels; a finer skin did not fit in memory.
- Yield and the 0.75 mm candidate not run; early flush-test code lost.

### Reproducibility (commands as run; `ccx223` = conda env with scipy)
/usr/bin/python3 scripts/build_stage3_solid.py spar # and: ribs
/usr/bin/python3 scripts/build_stage3_skin_geom.py # imprint + face list
/usr/bin/python3 scripts/mesh_stage3_solid.py ribs 2.5 # no-skin baselines
/usr/bin/python3 scripts/mesh_stage3_skin_raw.py 4.0 8.0 # raw Netgen mesh
python3 scripts/build_stage3_skin_inp.py mesh/stage3_skin/raw_s4_k8.npz 1.0 372.78 3.8 # ccx223
/usr/bin/python3 scripts/run_ccx_watchdog.py simulation/stage3_skin_s4_k8_t1_n3p8 8
/usr/bin/python3 scripts/extract_stage3_skin.py simulation/stage3_skin_s4_k8_t1_n3p8
python3 scripts/make_buckle_inp.py <static.inp> <buckle.inp> # then watchdog
python3 scripts/extract_stage3_buckle.py <buckle job>
/usr/bin/python3 scripts/analyze_buckle_modes.py <buckle job>
/usr/bin/python3 scripts/check_buckle_vs_plate.py <buckle job> <static job> <factor>
/usr/bin/python3 scripts/verify_skin_stress.py <job> # also check_deflection_profile.py,
# check_skin_load_variants.py
TE_MODE=edge /usr/bin/python3 scripts/mesh_stage3_skin_raw.py 4.0 8.0 <te_h> 800
/usr/bin/python3 scripts/te_hotspot_report.py <job>
/usr/bin/python3 scripts/collect_stage3_results.py
python3 scripts/plot_stage3_summary.py; python3 scripts/render_stage3_figures.py

Variant inputs for the load check were made by deleting the `*DLOAD` /
`ESKIN, GRAV` lines from the n = 1 input with `sed`. Meshes, solver
files and viewer files are not in git (regenerated by the scripts; see
`.gitignore`).

## Stage 4: lift distributed over the skin (analytic placeholder pressure)

Stages 1-3 applied the half-wing lift as nodal forces on the spar's top-cap
face. Stage 4 applies the same total lift to the skin. **The pressure
distribution is an analytic placeholder, not an aerodynamic prediction**
(like the 24.5 N load of Stages 1-2); CFD-derived pressure belongs to
Project 04.

### Decisions (locked by the user)

- Chordwise distribution: **B**, thin-airfoil flat-plate loading
  dCp proportional to sqrt((1 - xi)/xi), xi from the leading edge. Centre of
  pressure at 25% chord, 5% chord ahead of the spar axis (30%).
- Split: **antisymmetric**, suction on the upper skin and pressure on the
  lower skin, each carrying half of the net lift.
- Policy: if the critical buckling factor at t = 1.0 mm fell below 1.5
  under the distributed load, report it and ask the user; do not reopen the
  thickness lock. (Not triggered.)
- Unchanged: geometry (spar v2, six ribs at 200 mm), t = 1.0 mm, material,
  limit-load case (n = 3.8, 372.78 N per half-wing, gravity x 3.8; the 20 kg
  aircraft mass and the 290 MPa yield remain ASSUMPTIONS), meshes (solid 4 mm,
  skin 8 mm; skin 6 mm for convergence), fixed root, requirements
  (buckling factor >= 1.5, yield factor of safety >= 1.5).
- Implementation choices made by the assistant (not engineering decisions):
  consistent nodal forces instead of `*DLOAD ... P` (the sign convention of
  the CalculiX shell pressure was not confirmed), 64 sub-triangles per
  element, exact antiderivatives for the singular leading edge, and
  normalisation of the net vertical force.

### Input format (read before writing the builder)

The production `.inp` (13.07 MB): `*NODE` `id, x, y, z` (150,322 nodes),
`*ELEMENT, TYPE=C3D10` (44,899), `*ELEMENT, TYPE=S6` (34,698, ids from
44900; rows are three corner nodes followed by the midside nodes of edges
(a,b), (b,c), (c,a)), then `*DLOAD` (two `GRAV` lines), `*CLOAD` (3,349 lines,
all dof 2, sum 372.78 N), `*NODE FILE`, `*EL FILE`, `*END STEP`. The
builder asserts on this structure.

### Load builder (scripts/apply_stage4_pressure_load.py)

- Net load per unit planform area p_net = q(z) f(xi) / c(z), with
  q(z) = q0 sqrt(1 - (z/1000)^2), q0 = L / (1000 pi/4), f = (2/pi)
  sqrt((1-xi)/xi), c(z) = 300 - 0.15 z. Each surface carries a normal pressure
  p_net / 2; because n_y dA is the planform area element, the vertical
  component per planform area is exactly p_net / 2 on each surface.
- Upper/lower is decided by the sign of the element centroid y. Force is
  +n (outward) on the upper skin and -n on the lower skin.
- Each triangle is split into 64 sub-triangles. The mean of f over the xi
  range of a sub-triangle comes from the exact antiderivative
  G(xi) = (2/pi)(theta + sin theta cos theta), theta = asin sqrt(xi); the
  mean of q over its z range from the exact integral of sqrt(1 - s^2). The
  force is distributed to the six S6 nodes with the quadratic shape
  functions at the sub-triangle centroid (a consistent nodal load).
- Net vertical force normalised to LIFT_N (the scale factor is printed and
  gated). The old `*CLOAD` block is replaced (not duplicated); modes: `full`
  (pressure + existing gravity), `lift` (no `*DLOAD`), `grav` (no `*CLOAD`).

### Pre-solve gates (all passed on the first run)

| Check | 4/8 mm | 4/6 mm |
|---|---|---|
| raw net Fy before scaling [N] | 372.767289 | 372.770428 |
| scale factor (gate 1%) | 1.000034 | 1.000026 |
| net Fx / Fz after scaling [N] | -0.0381 / -0.0038 | -0.0387 / -0.0038 |
| root moment sum(Fy z) vs 158,212.746 N mm (gate 0.5%) | -0.0001% | +0.0002% |
| sum(Fy x) about the spar axis vs 25%-chord value (gate 2%) | +0.041% | +0.031% |
| chordwise centre of pressure per 100 mm bin (gate +-0.01) | 0.2493-0.2506 | 0.2495-0.2504 |
| outward-normal orientation test | 0 failures | 0 failures |
| loaded nodes / on shell-only nodes | 69,840 / 95% (354.80 N) | 89,500 / 95% (354.79 N) |

The values in the written `.inp` files were re-read independently
(scripts/collect_stage4_results.py): net Fy 372.780000 N, the same moment and
torque errors. The residual sum(-Fx y) = 425 N mm (about 10% of the
sum(Fy x) torque) is the horizontal component of the normal pressure, which
is real for a pressure acting normal to the skin; it is not part of the
hand twist estimate.

**Downward-force elements (inherited trailing-edge closure).** 46 skin
elements (41 on the 4/6 mesh; the 4/6 case was not examined separately) have
a force vector pointing downward. scripts/diag_stage4_negy_elements.py (4/8
mesh): all lie at xi = 0.9937-0.9982, share the normal (0.196, -0.980,
0.020), have a combined area of 177.3 mm^2 (0.04% of the skin) and a signed
planform area of -173.8 mm^2. They are the slanted closing facet of the
Stage 2 `airfoil_outline()` (upper surface to the true trailing-edge point,
lower surface one point earlier), which faces down but whose centroid has
y > 0, so it is treated as upper skin. The loading goes to zero at the
trailing edge; the load there was not computed (estimated far below
0.01 N). Accepted, documented; the rule was not changed.

### Solves (CalculiX 2.23, watchdog limit 8 GB; all exit code 0, 0 warnings)

| Job | Wall time | Peak RSS |
|---|---|---|
| 4/8 static, full | 92.8 s | 4.35 GB |
| 4/8 static, lift only / gravity only | 91.1 s / 88.8 s | 4.35 / 4.36 GB |
| 4/8 buckling (6 modes) | 214.1 s | 4.73 GB |
| 4/6 static, full | 137.8 s | 5.53 GB |
| 4/6 buckling (6 modes) | 323.5 s | 5.96 GB |

Both 4/6 runs took about 1.8 times as long as the Stage 3 runs on the same
mesh (72 s and 185 s) at similar memory; the cause was not investigated.

### Results

**Tip deflection (spar axis) and superposition.** 4/8: full 4.197839 mm
(Stage 3: 4.202998, -0.123%); lift only 5.131732 mm; gravity only
-0.933894 mm; sum 4.197838 mm. The lift-only value is within 0.12% of the
Stage 3 value inferred from the spar-cap-loaded model (15.215 x 0.33767 =
5.138 mm); gravity-only matches Stage 3 (-0.9338 mm). 4/6: 4.197943 mm
(Stage 3: 4.203088, -0.122%). The exact superposition checks the linear
solutions; the agreement of the lift-only value shows that the load on
shell-only nodes is transmitted in vertical bending. It does not verify the
chordwise distribution or the antisymmetric split (vertical bending depends
only on the net force and the root moment).

**Stress** (scripts/compare_stage4_stress.py; Stage 3 definitions unchanged,
design stress = maximum von Mises over expanded skin nodes with dTE >= 3 mm
and z > 10 mm; expanded-node interpretation verified near the root only):

| | S3 4/8 | S4 4/8 | S3 4/6 | S4 4/6 |
|---|---|---|---|---|
| design stress [MPa] | 32.853 | 32.095 | 32.843 | 32.087 |
| skin p99 [MPa] | 18.693 | 18.510 | 18.336 | 18.230 |
| skin maximum (TE vertex) [MPa] | 231.8 | 229.5 | 231.4 | 229.2 |
| solid-only maximum [MPa] | 82.9 | 82.2 | 83.7 | 82.9 |

Mean absolute change of von Mises over the design nodes 0.68 MPa (4/8) and
0.61 MPa (4/6); largest local changes +4.08 MPa at (-59.1, 9.8, 399.3) and
-3.65 MPa at (63.3, -11.7, 398.5) (4/8), next to rib 400. On the 4/6 mesh the increase is the same (+4.08 MPa
at the same position); the largest decrease is -3.63 MPa at (68.2, 11.3,
398.3), a different node on the opposite surface. A rough estimate (pressure about 0.0015 MPa on a panel about
90 mm wide and 1 mm thick gives a few MPa of local panel bending) is of the
same order; it was not checked further.

**Where the design maximum sits** (scripts/diag_stage4_design_nodes.py).
In both load cases and on both meshes it is the tip leading-edge node
(x = -45, y = 0, z = 1000), the nose of the free skin edge at the tip rib,
with six of the top ten nodes lying within 0.4 mm of it; the next-highest nodes are at
dTE = 3.1-3.7 mm beside ribs 1000, 800, 600 and 400, i.e. just outside the
3 mm trailing-edge exclusion band. Excluding z >= 995 mm gives 29.33 (Stage 3)
and 29.75 MPa (Stage 4); also excluding dTE < 6 mm gives 24.76 and 24.38 MPa.
These cutoffs were chosen after seeing the data (sensitivity rows, not a
rule). Whether the tip-nose value is converged was not tested. The
conclusion does not depend on the definition: Stage 3 and Stage 4 agree
within about 2%, the factor of safety against the assumed 290 MPa yield is
8.8 (Stage 3) and 9.0 (Stage 4) with the official definition, and stress
does not govern at 1.0 mm. The Stage 3 value (about 33 MPa) is therefore a
local number, not the bending field (p99 about 18.5 MPa); Stage 3 did not
look at where it sits.

**Buckling** (`*BUCKLE`, 6 modes; factor multiplies the whole applied load):

| | 4/8 | 4/6 |
|---|---|---|
| Stage 3 critical factor | 3.1388 | 3.1226 |
| Stage 4 modes 1-6 | 3.1938, 3.3658, 4.3522, 4.6241, 4.7899, 5.0939 | 3.1775, 3.3453, 4.3237, 4.5750, 4.7573, 5.0658 |
| Stage 4 / Stage 3, critical | +1.75% | +1.76% |

Requirement 1.5: PASS (+113% / +112%). All six factors fall 0.55-1.06% from
the 8 mm to the 6 mm skin. Modes 1-5 are local out-of-plane buckles of the
upper skin aft of the spar in the 0-200 mm bay (mode 1 peaks at
(87.0, 13.3, 58.8) on the 6 mm mesh against (89.6, 13.1, 60.2) on the 8 mm
mesh and (91.5, 12.1, 132.3) in Stage 3); mode 6 lies in the 200-400 mm
bay. All hot nodes are on the upper skin, none within 20 mm of the leading or
trailing edge or within 3 mm of a rib plane. Note that in the output of
analyze_buckle_modes.py "mode 1" is the static preload and modes 2-7 are the
buckling modes 1-6. The cause of the +1.75% was not isolated; von Mises on
the upper aft panel (mean -2.2%, p95 +7.4%, maximum +2.8% from Stage 3 to
Stage 4) is not the compressive stress that governs the buckle, and no
attribution is made. The shell-node nodal loads are fixed-direction (not
follower) loads in the preload.

**Twist** (scripts/extract_stage4_twist.py: least-squares fit of Uy against x
over the solid original nodes in each rib plane; theta = -slope, nose-up
positive; the ribs are not rigid, rms fit residual 0.0034-0.0046 mm for the
full case, at most 0.0073 mm over all cases shown):

| z [mm] | 200 | 400 | 600 | 800 | 1000 |
|---|---|---|---|---|---|
| Stage 4 full [1e-4 rad] | 1.444 | 2.694 | 3.545 | 3.934 | 3.987 |
| Stage 4 lift only | 1.109 | 2.019 | 2.530 | 2.607 | 2.503 |
| Stage 4 gravity only | 0.335 | 0.675 | 1.015 | 1.327 | 1.485 |
| Stage 3 (spar-cap lift) | 0.862 | 1.463 | 1.750 | 1.718 | 1.592 |
| hand estimate, lift only | 0.890 | 1.711 | 2.398 | 2.858 | 2.999 |

Root control 4e-8 rad in every case; lift-only plus gravity-only equals the
full case (z = 200: 1.109 + 0.335 = 1.444; z = 1000: 2.503 + 1.485 = 3.988).
4/6 full case: 3.976e-4 rad at the tip (-0.28%). Tip twist 0.0229 degrees;
a tip leading edge 45 mm ahead of the spar rises about 0.018 mm, 0.4% of the
bending deflection: twist is not a design driver here. The gravity
contribution (37% of the tip total) was not in my prediction; presumably the
centre of gravity of skin and solid lies aft of the spar axis (not computed).

The hand estimate (scripts/estimate_stage4_twist_bredt.py) is a single closed
cell, t = 1.0 mm, ribs and web ignored, G = E/(2(1+nu)) = 27,481 MPa, enclosed
area 0.08171 c^2, perimeter 2.0395 c, GJ 9.71e9 (root) and 1.21e9 N mm^2
(tip), shear centre ASSUMED at the spar axis (torque arm 0.05 c). It was run
after the FEA values were seen, so it is not a blind prediction. It is
linear in the assumed arm (0.35 c would double it, 0.25 c gives zero).
Hand/FEA (lift only): 0.80, 0.85, 0.95, 1.10, 1.20 at z = 200 ... 1000; the
FEA twists faster inboard and levels off outboard; the cause (ribs, web,
shear lag, shear centre) was not tested. This is a plausibility check, not a
verification.

**Skin-mesh convergence (8 to 6 mm, two levels):** tip deflection +0.003%,
tip twist -0.28%, design stress -0.02%, p99 -1.5%, buckling factor -0.51%
(Stage 3, same refinement: -0.5%). No asymptotic-convergence claim.

### Predictions written down before the runs, and what happened

| Prediction | Outcome |
|---|---|
| tip deflection within +-3% of Stage 3 | -0.12% (hit) |
| design stress within +-10% | -2.3% (hit) |
| buckling 2.8-3.4 | 3.1938 (hit) |
| twist from the lift 2-4e-4 rad | 2.50e-4 (hit); gravity share not predicted |
| 4/6 tip deflection within 0.1% | +0.0025% (hit) |
| 4/6 buckling 3.16-3.18 | 3.1775 (hit) |
| 4/6 twist within a few % | -0.28% (hit) |
| both cumulative load-distribution differences well under 1% | spanwise 0.17% (hit); chordwise 1.36% (MISS) |
| 4/6 static solve about 75 s | 137.8 s (miss) |
| 4/6 buckling solve 3-5 minutes | 5.4 minutes (slightly over) |

For the chordwise miss, the hypothesis (not tested) is that nodal loads sit at
element nodes about 8 mm apart (xi of about 0.027 at the root) while the
cumulative loading is steep at the leading edge, where about 22% of the lift
lies within the first 0.03 of chord; the centre of pressure and the torque are
unaffected (0.2493-0.2506; +0.04%). The leading-edge pressure peak is
resolved no better than the skin mesh allows.

### Retractions and corrections (chronological)

1. I told the user to expect Bredt-Batho agreement with the FEA twist to
   about 20-30%. That had no basis without a known shear centre (the arm is
   only 5% chord and the estimate is linear in it); retracted. The
   documented statement is plausibility only.
2. Stage 3 reported a "design skin stress" of about 33 MPa without examining
   where the maximum lies. Stage 4 found it at the tip leading-edge node and
   just outside the trailing-edge band; the wording in the README was
   qualified.
3. Predictions that missed (above): the gravity share of the twist, the
   chordwise cumulative difference, and the solve times.
4. My first README insertion script dropped the `## Limitations` heading (the
   replacement text omitted it) and contained two wording errors (the twist
   range was a hand/FEA ratio mislabelled as a percentage difference; "factor of
   safety 9 or more" was 8.8 in Stage 3). Found with `grep -n '^#'` and
   corrected with guarded patches. A stray token also appeared in a pasted
   patch script (removed). The plot script was patched once because the
   twist figure's legend covered a curve and its title misstated the fit
   residual (0.003-0.005 mm; the Stage 3 curve reaches 0.0072 mm).
5. The 46-element diagnosis was first run as an inline here-document and was
   saved as scripts/diag_stage4_negy_elements.py afterwards (re-run and
   compared with the first output) so that it regenerates from the repository.
6. The first version of this section dropped the `## Status Log` heading
   (the same mistake as in the README, repeated), said that the largest local
   stress changes were identical in both meshes (only the increase is; the
   decrease is at a different node) and that ranks 5-10 of the design nodes
   cluster at the tip nose (six of the top ten do). Found on re-reading the
   text against the outputs and corrected.

### Findings

- Distributing the same lift over the skin as an analytic placeholder
  pressure changes the tip deflection by -0.12%, the design stress by -2.3%
  and the critical buckling factor by +1.75% (3.14 to 3.19); the 1.0 mm lock
  is unchanged and the policy for a factor below 1.5 was not triggered.
- The wing twists nose-up by 0.023 degrees at the tip (lift 2.50e-4 rad,
  gravity 1.49e-4 rad); small, not a design driver.
- The stress result does not depend on the definition used, and stress does
  not govern at 1.0 mm (factor of safety about 9 against the assumed yield).

### Figures

figures/stage4_load_distribution.png (applied nodal load against the analytic
span and chord distributions; cumulative differences 0.17% and 1.36%) and
figures/stage4_twist_span.png (from scripts/plot_stage4_figures.py; the dashed
curve is the hand estimate). Tables: results/stage4_{static_summary,twist,
buckling_summary,load_check}.csv and results/stage4_summary.json
(scripts/collect_stage4_results.py).

### Limitations and open items

- The applied pressure is an analytic placeholder: the antisymmetric split,
  the thin-airfoil chordwise shape and the elliptical span are assumptions.
  The 20 kg aircraft mass and the 290 MPa yield remain unsourced.
- Shell-node `*CLOAD` was verified through vertical bending only (lift-only
  tip deflection, exact superposition), not for the chordwise or
  antisymmetric detail; reaction forces cannot check equilibrium for
  expanded shells.
- The nodal load follows the chordwise distribution to 1.4% (maximum
  cumulative difference), limited near the leading edge by the 8 mm skin
  mesh; the 4/6 mesh figure was not produced.
- The twist reference is a hand estimate with an assumed shear centre; the
  FEA/hand differences by station are unexplained.
- The design stress is set by local locations (tip leading-edge node; band
  beside the trailing-edge exclusion); convergence of the tip-nose value was
  not tested.
- 46 (4/8) and 41 (4/6) elements of the inherited trailing-edge closure are
  treated as upper skin with a downward force; the 4/6 set was not examined.
- Buckling is linear with fixed-direction loads; the cause of the +1.75% change
  was not isolated; two mesh levels only.
- The unexplained skin-vs-beam stress/curvature deficit of Stage 3 was not
  re-examined and is still open. The slow 4/6 solves were not investigated.

### Reproducibility (commands as run; `ccx223` = conda env with scipy)
/usr/bin/python3 scripts/apply_stage4_pressure_load.py simulation/stage3_skin_s4_k8_t1_n3p8.inp simulation/stage4_s4_k8_t1_n3p8.inp full
# also: ..._lift.inp with mode lift, ..._grav.inp with mode grav; for 4/6 use stage3_skin_s4_k6_t1_n3p8.inp -> stage4_s4_k6_t1_n3p8.inp
/usr/bin/python3 scripts/run_ccx_watchdog.py simulation/stage4_s4_k8_t1_n3p8 8 # also _lift, _grav, and the k6 job
/usr/bin/python3 scripts/extract_stage3_skin.py simulation/stage4_s4_k8_t1_n3p8
python3 scripts/make_buckle_inp.py simulation/stage4_s4_k8_t1_n3p8.inp simulation/stage4_s4_k8_t1_n3p8_buckle.inp 6 # ccx223
/usr/bin/python3 scripts/run_ccx_watchdog.py simulation/stage4_s4_k8_t1_n3p8_buckle 8
/usr/bin/python3 scripts/extract_stage3_buckle.py simulation/stage4_s4_k8_t1_n3p8_buckle
/usr/bin/python3 scripts/analyze_buckle_modes.py simulation/stage4_s4_k8_t1_n3p8_buckle
/usr/bin/python3 scripts/extract_stage4_twist.py simulation/stage4_s4_k8_t1_n3p8 simulation/stage4_s4_k8_t1_n3p8_lift simulation/stage4_s4_k8_t1_n3p8_grav simulation/stage3_skin_s4_k8_t1_n3p8
/usr/bin/python3 scripts/estimate_stage4_twist_bredt.py
/usr/bin/python3 scripts/compare_stage4_stress.py simulation/stage3_skin_s4_k8_t1_n3p8 simulation/stage4_s4_k8_t1_n3p8
/usr/bin/python3 scripts/diag_stage4_design_nodes.py simulation/stage3_skin_s4_k8_t1_n3p8 simulation/stage4_s4_k8_t1_n3p8
/usr/bin/python3 scripts/diag_stage4_negy_elements.py
/usr/bin/python3 scripts/collect_stage4_results.py
python3 scripts/plot_stage4_figures.py # ccx223


## Stage 5: rib spacing (reduced scope: 200 mm against 100 mm)

Locked question: how does uniform rib spacing (100, 50, 20 mm against the
200 mm baseline) affect tip deflection, skin buckling resistance and margin?
Outcome: only 200 mm and 100 mm could be solved at the validated mesh on this
machine (11 GB RAM). **The 50 mm and 20 mm cases were not run.**

### Decisions (locked by the user)

- D1 = A: skin t = 1.0 mm at every spacing (the thickness lock); the 0.5 and
  0.75 mm candidates were not studied.
- D2: the Stage 4 analytic placeholder pressure is the baseline load; the
  200 mm baseline is the Stage 4 result (critical buckling factor 3.1938).
- D3: rib model unchanged (solid, 1.5 mm thick, no lightening holes).
- Spacings uniform, root and tip ribs included: 11 / 21 / 51 ribs.
- Scope (2026-10-07): option A, drop 50 and 20 mm now and extend after a RAM
  upgrade. Presented and not chosen: B (coarser solid mesh, would need the
  200 mm baseline re-verified), C (ribs as shells, contradicts D3), D (dense
  ribs over part of the span only), E (iterative solver), F (linear elements).
- Assistant implementation choices: the `RIB_SPACING` environment variable
  (default 200 = unchanged behaviour and file names), the `_rs<N>` suffix on
  all outputs, `.gitignore` rules for the regenerable `_rs*` geometry, and the
  `NETGEN_KW` environment variable for extra `GenerateMesh` keywords.

### Parametrised pipeline (guarded patches, all anchors checked)

- build_stage3_solid.py: `RIB_STATIONS` from `RIB_SPACING` (must divide the
  1000 mm span), output name with the suffix.
- build_stage3_skin_geom.py: input/output names with the suffix; the expected
  solid volume at other spacings comes from the builder's own numerical
  integration (relative tolerance 1e-7) instead of the 6-rib constant, and a
  sidecar `stage3_solid_volume_rs<N>.txt` carries it to the input builder.
- mesh_stage3_skin_raw.py: names with the suffix; the face-count check compares
  the classified faces with the imprint JSON (388/446/72 is kept at 200 mm);
  `NETGEN_KW` (e.g. `closeedgefac=None`) adds keywords and a tag suffix.
- build_stage3_skin_inp.py: the suffix is taken from the raw-mesh file name and
  the sidecar volume is read.
- extract_stage3_skin.py: refuses to run if the job name's `_rs<N>` disagrees
  with `RIB_SPACING` (its beam reference depends on the rib stations).
- Regression check at 200 mm: solid volume 128,582.1216 mm^3, 834 skin pieces,
  tracked `stage3_skin_faces.json` byte-identical; the rebuilt tracked solid
  STEP differs only in its header timestamp (restored with `git checkout`).

### Geometry (builder and imprint; all validations PASS)

| Spacing | Ribs | Solid volume [mm^3] | Mass | Skin pieces (shared / free) | Solid boundary faces / non-skin | Imprint time |
|---|---|---|---|---|---|---|
| 200 mm | 6 | 128,582.1216 | 357.46 g | 834 (446 / 388) | 518 / 72 | |
| 100 mm | 11 | 160,188.1389 | 445.32 g | 1590 (814 / 776) | 956 / 142 | |
| 50 mm | 21 | 223,469.4464 | 621.25 g | 3098 (1548 / 1550) | 1830 / 282 | 2 m 53 s |
| 20 mm | 51 | 413,368.7873 | 1149.17 g | 7630 (3754 / 3876) | 4456 / 702 | 13 m 26 s |

The skin area after the imprint is 459,622.6823 mm^2 in all cases (the
envelope does not depend on the ribs) and the shared areas match the
independent rim-plus-cap estimates (e.g. 23,333.3512 against 23,333.3506 mm^2
at 100 mm). At 20 mm the imprint has four sliver pieces: two free pieces of
0.0037 mm^2 and two of 0.7228 mm^2 (scripts/diag_rs_slivers.py). They are
remnants of the inherited wedge piece (the cap-edge line x = -4 crossing a
scaled facet line near z = 742.3, apex angle about 0.18 degrees) cut by the
ribs at z = 720 and 740; the other spacings have no piece below 1.14 mm^2.

### Meshes (solid 4 mm / skin 8 mm, Netgen; counts from the input builder)

| Spacing | C3D10 | S6 | Nodes | Elements (share of 230k) | Equation estimate (share of 1.18M) |
|---|---|---|---|---|---|
| 200 mm | 44,899 | 34,698 | 150,322 | 79,597 (34.6%) | 660,486 (56.0%) |
| 100 mm | 68,050 | 46,028 | 214,545 | 114,078 (49.6%) | 921,135 (78.1%) |
| 50 mm | 112,518 | 65,800 | 334,992 | 178,318 (77.5%) | 1,401,108 (**118.7%**, over) |
| 20 mm | not meshed; about 246,000 tets extrapolated (6, 11, 21 ribs: +4630 and +4447 tets per added rib) | | | | |

The equation estimate is 3 per node plus 3 more per skin node (an estimate). A
rough extrapolation from the Stage 3 solves (880k equations needed 5.7 GB)
puts the 50 mm static solve at about 9 GB, above the 8 GB watchdog limit and
close to the 11 GB host, which was crashed at about 1.18 million equations in
Project 01. The 50 mm builder run stops at the ceiling check ('size exceeds
ceiling'); the mesh itself is built and checks out (skin boundary edges only at
the root and tip ends, tet volume +0.0000%, 0 tet-boundary faces without a
surface triangle).

### Netgen investigation at 50 mm (chronological)

1. The first mesh run at 4/8 (timeout 2400 s) ended with exit 124 after 41 min
   26 s, using only 634 MB (the 100 mm mesh takes 54 s and 380 MB). The timeout
   was far too generous; the log held only a Netgen line ' ERROR: Problem in
   Surface mesh generation' (Python's buffered output was lost on the kill).
2. An unbuffered 300 s re-run printed 'Faces: 3380 (free 1550, shared 1548 +
   non-skin 282)' as predicted and no 'GenerateMesh OK'; my `grep -v` filter
   hid the ERROR line (it begins with a space). A coarse trial at solid 6 /
   skin 10 mm failed the same way: not a mesh-size effect.
3. scripts/diag_netgen_faces.py (modes map, verbose, face, edges, single, topo)
   verified that Netgen's face number k is `shape.faces[k-1]` (1732 of 1732
   indices on the 100 mm mesh). Message levels 5 and 9 print identical logs.
   scripts/parse_netgen_face_log.py: the ERROR came while face 2603 was in
   progress ('give up with qualclass 201', 'number of frontlines = 9'), Netgen
   continued, and it spun on face 2986 (started twice, never finished; a front
   of more than 2,100 lines on a 110 mm^2 face), so faces 2987 to 3380 never
   started.
4. **Misreading, retracted:** I first read the 'faces = N trials = N elements =
   N' lines as the volume mesher; they are the 2D surface mesher's (the first
   one is at log line 14, interleaved with the face lines).
5. The faces are ordinary free skin pieces (2603: a quadrilateral 6.8 x 48.7 mm
   at the trailing edge; 2986: a strip 2.3 x 48.5 mm, the root-end bay of the
   wedge family; 2987: the wedge triangle with a 0.13209 mm base). Each of them
   alone meshes in 0.03 to 0.05 s (12, 17 and 28 triangles). Topology is
   conforming (no edge used by one face; 2 or 3 faces per edge at 100 and 50
   mm). The shortest edge, 0.13209 mm (the wedge base), is the same at 100 and
   50 mm. **API slip:** my first `edges` mode used an unverified `e.start.p`
   and printed None; it was rewritten with `e.center` and `e.mass`.
6. The `GenerateMesh` documentation lists the OCC option
   `closeedgefac: Optional[float] = 2.` ('if None it is disabled'), on by
   default. With `closeedgefac=None` all 3380 faces are meshed with no ERROR;
   with 1.0 the run still fails (2812 faces done, one ERROR). The full 50 mm
   mesh with the option takes 88.4 s (2 min 39 s wall, 605 MB).
7. With `closeedgefac=None` the raw meshes at **100 mm and 200 mm are
   array-for-array identical** to the default ones (scripts/compare_raw_meshes.py:
   P, tets, tris, tri_idx, in_tet_bnd), so no earlier result needed re-meshing;
   the option only matters for the 50 mm geometry. Why the default refinement
   fails there was not determined.

### Tet quality (scripts/check_tet_quality.py; mean ratio q, 1 = regular tet)

| Mesh | Tets | q min | q p1 | q median | tets q < 0.05 / < 0.10 | Max edge ratio |
|---|---|---|---|---|---|---|
| 200 mm (Stage 3/4 mesh) | 44,899 | 0.0311 | 0.483 | 0.818 | 6 / 35 | 46.2 |
| 100 mm | 68,050 | 0.0301 | 0.529 | 0.822 | 7 / 32 | 63.4 |
| 50 mm (not solved) | 112,518 | 0.0306 | 0.545 | 0.829 | 4 / 28 | 23.0 |

No tet below 0.02. The worst tets in all three meshes are at x about -3.8 to
-4.0 mm, |y| about 10.5 to 10.75 mm, z about 745 to 749 mm (just beyond the
wedge apex), away from the reported hot spots. The smallest-volume tets
(0.000499 mm^3 at 100 and 50 mm, near the 0.13 mm wedge base) are not among
the worst by q. This also characterises the Stage 3 and 4 production mesh,
which had only been checked for positive volume. The effect of the poorest tets
was not tested.

### Solves (CalculiX 2.23, watchdog limit 8 GB; all exit code 0, 0 warnings)

| Job (100 mm, 4/8) | Wall time | Peak RSS |
|---|---|---|
| static, full | 134.2 s | 6.19 GB |
| static, lift only / gravity only | 139.3 s / 136.3 s | 6.19 / 6.17 GB |
| buckling (6 modes) | 325.6 s | 6.71 GB |

Pre-solve load gates at 100 mm: scale factor 1.000055, root moment +0.0002%,
sum(Fy x) +0.073%, centre of pressure 0.2499 to 0.2501 in every 100 mm bin,
94.1% of the load on shell-only nodes, 88 elements with a downward force (46 at
200 mm; the ratio 1.91 equals (2 x 11 - 1) / (2 x 6 - 1), consistent with the
same trailing-edge closure facet cut into more strips; not separately
diagnosed).

### Results (scripts/collect_stage5_results.py; results/stage5_comparison.csv)

| Quantity | 200 mm | 100 mm | Change |
|---|---|---|---|
| Solid mass / skin mass / total | 0.3575 / 1.2778 / 1.6352 kg | 0.4453 / 1.2778 / 1.7231 kg | +5.37% total |
| Tip deflection full / lift only / gravity only [mm] | 4.197839 / 5.131732 / -0.933894 | 4.046027 / 4.994565 / -0.948538 | -3.62% / -2.67% / +1.57% |
| Beam with rib-section inertia, full [mm] | 4.194180 | 4.129242 | |
| FEA against that beam | +0.09% | -2.02% | |
| Tip twist (nose-up) | 3.987e-4 rad = 0.0228 deg | 4.127e-4 rad = 0.0236 deg | +3.5% |
| Skin p99 von Mises | 18.510 MPa | 18.592 MPa | +0.4% |
| Skin TE-vertex maximum (excluded zone) | 229.5 MPa | 248.2 MPa (rib 700) | |
| Design stress, dTE >= 3 mm (official) | 32.095 MPa | 37.718 MPa | +17.5% |
| Design stress, dTE >= 15 mm | 32.095 MPa | 30.719 MPa | |
| Variant z < 990 and dTE >= 6 (post hoc) | 24.380 MPa | 25.143 MPa | +3.1% |
| Critical buckling factor | 3.1938 | 3.9123 | +22.5% |
| Buckling modes 1 to 6 | 3.1938, 3.3658, 4.3522, 4.6241, 4.7899, 5.0939 | 3.9123, 4.3297, 5.2570, 5.4994, 5.6508, 5.9790 | +17 to +29% |

- Superposition holds exactly (100 mm: 4.994565 - 0.948538 = 4.046027 mm; twist
  at z = 200: 1.211e-4 + 0.325e-4 = 1.536e-4 rad).
- The extractor beam reference (no rib inertia) moved from 4.218995 to 4.177824
  mm (-0.98%); the FEA moved -3.62%.
- scripts/beam_with_rib_inertia.py adds, over each rib's 1.5 mm, the second
  moment of the filled airfoil outline (I_solid: root 550,634, mid 174,224, tip
  34,415 mm^4; I_ring 115,081, 49,646, 15,372; I_spar at the root 17,577). Its
  no-rib variant reproduces `beam_tip` exactly. With rib inertia the beam
  changes by -0.587% (200 mm) and -1.162% (100 mm), i.e. -0.58% between the two
  spacings, against -2.67% in the FEA lift response. At 200 mm the rib-aware
  beam agrees with the FEA to +0.09% (full) and +0.17% (lift only); the earlier
  -0.50% was mostly the rib sections. At 100 mm the FEA is 2.0% stiffer, in both
  the lift and the gravity responses, so it is a structural stiffness effect.
  Interpretation (not tested): closely spaced ribs suppress section distortion
  and shear lag in the skin over the whole span.
- Buckling modes at 100 mm: all upper skin aft of the spar, none within 20 mm of
  the leading or trailing edge or within 3 mm of a rib plane. Mode 1 peaks at
  (79.4, 14.5, 49.2) with hot nodes at z 29 to 73 (one lobe in the first bay);
  mode 2 at z = 149.7 (hot region z 125 to 172); all six modes lie at z < 200
  (at 200 mm mode 6 was in the 200 to 400 mm bay). The bay labels printed by
  analyze_buckle_modes.py are fixed to the 200 mm stations and are wrong at 100
  mm. In the plan-view figure the 200 mm mode shows two bright spots in the root
  bay (peaks near z = 60 and 132); the cause of the +22.5% was not isolated.
- Design stress: the top design nodes at 100 mm are at 3.0 to 3.5 mm from the
  trailing-edge vertex at the ribs (z = 700.8, 899.2, 900.8, 699.2, 499.2, ...),
  plus the tip leading-edge node (30.719 MPa). The dTE >= 3 mm metric is
  mesh-sensitive (the Stage 3 refinement study showed the 3 to 6 mm ring rising
  from 29.3 to 35.7 MPa; the 100 mm mesh is finer there, free-skin median edge
  3.98 mm against 5.62 mm) and is not comparable across spacings. The
  factor of safety against the assumed 290 MPa is at least 7.7 (official) or
  about 11.5 (variant). The "|z-rib|" column of diag_stage4_design_nodes.py is
  tied to the 200 mm stations and is meaningless at 100 mm.

### Predictions written down before the runs, and what happened

| Prediction | Outcome |
|---|---|
| 100 mm solid about 160,000 mm^3 / about 445 g | 160,188 / 445.32 g (hit) |
| 100 mm piece count about 1,600 | 1,590 (hit by luck, see retraction 1) |
| 50 mm solid about 223,500 / 621 g; pieces about 3,100 | 223,469 / 621.25 g; 3,098 (hit) |
| 20 mm solid about 413,000 / 1.15 kg; pieces about 7,620 | 413,369 / 1.149 kg; 7,630 (hit) |
| 20 mm imprint takes 7 to 20 minutes | 13 min 26 s |
| 100 mm tets 55k to 90k | 68,050 (hit) |
| 50 mm tets about 114k; equations 1.3 to 1.45 M | 112,518; 1.40 M (hit) |
| 100 mm static 120 to 250 s, 6.0 to 6.3 GB | 134 s, 6.19 GB (hit) |
| 100 mm buckling about 6.6 GB, about 315 s | 6.71 GB, 326 s (hit) |
| Weight of the 100 mm solid at n = 3.8: 17.2 N | 16.60 N (arithmetic error) |
| Tip deflection about 4.15 mm (-1%) | 4.046 mm (-3.6%) (**miss**) |
| Buckling factor 3.1 to 3.5 | 3.91 (**miss**) |
| Skin p99 and design stress within a few % | p99 +0.4% (hit); official design stress +17.5% (**miss**) |
| Rib-section inertia explains -0.5 to -1% | -0.58% between the spacings (hit, low end) |
| Mode peaks near z = 60 and 49 mm | 60.2 and 49.2 mm (hit) |
| Tet quality: median 0.7 to 0.8; worst tets near x -4, z 700 | median 0.82 to 0.83; worst at z 745 to 749 (partial miss) |

### Retractions and corrections (chronological)

1. The 100 mm piece-count 'pattern' used 78 x 21 = 1,614; the product is 1,638.
   The two-point linear structure found afterwards (free = 77.6 (n - 1), shared =
   73.6 (n - 1) + 78) is only an extrapolation device, accurate to 0.13% at 50
   and 20 mm.
2. I said the builder took about four minutes; that came from file timestamps
   that include human time. Measured: 19.7 s (50 mm), 82.5 s (20 mm).
3. The 2,400 s Netgen timeout was too generous (41 minutes lost); later runs used
   short limits, `python3 -u` and size-guarded logs.
4. My `grep -v '^[[:space:]]'` filter hid the ERROR line of the unbuffered run.
5. The volume-mesher misreading (above).
6. The unverified `e.start.p` API in the first `edges` mode (above).
7. The 17.2 N weight prediction was an arithmetic error (16.60 N).
8. The Stage 3 plate-theory expectation that ribs help only below about 93 mm
   spacing, and my prediction for this run, were wrong at 100 mm.
9. The README sentence that a 6 mm skin at 100 mm 'would exceed the memory
   ceiling' was unsupported: the extrapolation gives about 1.09 M equations
   (below the ceiling), about 7.7 GB static and about 8.2 GB for buckling, at or
   above the 8 GB watchdog limit; it was not tried. Corrected in the README.
10. The official design-stress metric (dTE >= 3 mm) is mesh-sensitive and not
    comparable across spacings; the +17.5% is not a physical stress increase.
11. My guess for the worst-tet location (z about 700) was off; the apex is at z
    about 742 and the worst tets are at 745 to 749.

### Findings

- Halving the rib spacing adds 5.4% mass and raises the critical buckling factor
  by 22.5% (3.19 to 3.91; all six modes rise) and lowers the tip deflection by
  3.6%, leaves the field stress unchanged (p99 +0.4%) and raises the twist
  slightly (+3.5%). At t = 1.0 mm the buckling requirement was already met by
  +109%, so this margin is not needed at this thickness; thinner skins with
  denser ribs were not studied.
- About a fifth of the extra bending stiffness is the rib sections themselves;
  about 2% remains unexplained (suspected section-distortion restraint).
- The Stage 3 and 4 tip deflection agrees with a rib-aware beam to 0.1 to 0.2%.
- With two spacings there is no trend and no optimum can be named.

### Figures

figures/stage5_comparison.png (scripts/plot_stage5_figures.py) and
figures/stage5_buckling_mode1.png (scripts/render_stage5_modes.py, plan view of
buckling mode 1 for both spacings with rib lines).

### Limitations and open items

- Only two spacings; 50 and 20 mm not run (RAM). 50 mm needs about 9 GB (extrapolated)
  and `NETGEN_KW="closeedgefac=None"`; 20 mm needs a mesh strategy (its element
  count extrapolates above the ceiling) and attention to the sliver pieces.
- No skin-size or solid-size convergence at 100 mm (a 6 mm skin is estimated at
  about 7.7 GB static / 8.2 GB buckling, not tried).
- About 2% of the 100 mm bending stiffness, the small twist increase and the
  +22.5% buckling change are not traced to a mechanism.
- The dTE >= 3 mm design-stress metric is mesh-sensitive.
- Tet-quality tail at the wedge apex in all meshes; effect not tested. The
  inherited trailing-edge closure facet and the wedge pieces remain.
- The load is the analytic placeholder; 20 kg and 290 MPa are assumptions; the
  ribs are solid without lightening holes.
- Why Netgen's default close-edge refinement fails on the 50 mm geometry is not
  known; why the 100 mm solves ran slower than in Stage 3 is not known either.

### Reproducibility (commands as run; `ccx223` = conda env with scipy)
RIB_SPACING=100 /usr/bin/python3 scripts/build_stage3_solid.py ribs
RIB_SPACING=100 /usr/bin/python3 scripts/build_stage3_skin_geom.py
RIB_SPACING=100 NETGEN_KW="closeedgefac=None" /usr/bin/python3 -u scripts/mesh_stage3_skin_raw.py 4.0 8.0 # identical to the default mesh at 100 and 200 mm
python3 scripts/build_stage3_skin_inp.py mesh/stage3_skin/raw_s4_k8_rs100.npz 1.0 372.78 3.8 # ccx223
/usr/bin/python3 scripts/apply_stage4_pressure_load.py simulation/stage3_skin_s4_k8_rs100_t1_n3p8.inp simulation/stage4_s4_k8_rs100_t1_n3p8.inp full # also lift, grav
/usr/bin/python3 scripts/run_ccx_watchdog.py simulation/stage4_s4_k8_rs100_t1_n3p8 8 # also _lift, _grav
python3 scripts/make_buckle_inp.py simulation/stage4_s4_k8_rs100_t1_n3p8.inp simulation/stage4_s4_k8_rs100_t1_n3p8_buckle.inp 6 # ccx223
/usr/bin/python3 scripts/run_ccx_watchdog.py simulation/stage4_s4_k8_rs100_t1_n3p8_buckle 8
RIB_SPACING=100 /usr/bin/python3 scripts/extract_stage3_skin.py simulation/stage4_s4_k8_rs100_t1_n3p8 # also _lift, _grav; refuses without RIB_SPACING
RIB_SPACING=100 /usr/bin/python3 scripts/extract_stage3_buckle.py simulation/stage4_s4_k8_rs100_t1_n3p8_buckle
/usr/bin/python3 scripts/extract_stage4_twist.py simulation/stage4_s4_k8_rs100_t1_n3p8 simulation/stage4_s4_k8_rs100_t1_n3p8_lift simulation/stage4_s4_k8_rs100_t1_n3p8_grav simulation/stage4_s4_k8_t1_n3p8
/usr/bin/python3 scripts/compare_stage4_stress.py simulation/stage4_s4_k8_t1_n3p8 simulation/stage4_s4_k8_rs100_t1_n3p8
/usr/bin/python3 scripts/diag_stage4_design_nodes.py simulation/stage4_s4_k8_rs100_t1_n3p8
RIB_SPACING=100 /usr/bin/python3 scripts/beam_with_rib_inertia.py # and without RIB_SPACING for 200 mm
/usr/bin/python3 scripts/check_tet_quality.py mesh/stage3_skin/raw_s4_k8.npz mesh/stage3_skin/raw_s4_k8_rs100.npz mesh/stage3_skin/raw_s4_k8_rs50_closeedgefacnone.npz
/usr/bin/python3 scripts/compare_raw_meshes.py mesh/stage3_skin/raw_s4_k8.npz mesh/stage3_skin/raw_s4_k8_closeedgefacnone.npz
/usr/bin/python3 scripts/collect_stage5_results.py
python3 scripts/plot_stage5_figures.py; python3 scripts/render_stage5_modes.py # ccx223
# 50 mm (not solved): RIB_SPACING=50 on the build scripts and NETGEN_KW="closeedgefac=None" on the mesh script;
# Netgen diagnostics: scripts/diag_netgen_faces.py, parse_netgen_face_log.py, diag_rs_slivers.py, diag_short_edges.py
# extension after a RAM upgrade: repeat the 100 mm chain with RIB_SPACING=50 (then 20, after a mesh strategy for it)

## Status Log

- Specification defined and approved
- Stage 1 (spar alone): complete and verified; committed locally
  (1f64068). Bugs: density units (1000x), equilibrium sign. Addenda: the
  root-reaction residual is a reporting convention; stress verified.
- Stage 2 (+ ribs, placeholder spacing): complete and verified
  (documentation in this entry); ribs contribute weight, not bending
  stiffness
- Stage 3 (spar v2 + skin, thickness locked at 1.0 mm): complete and documented in this
  entry; committed locally (8e27d78)
- Stage 4 (load distributed over the skin as an analytic placeholder pressure): complete and
  documented in this entry; committed locally (see git log)
- Stage 5 / Phase A (rib spacing): complete in a reduced scope (200 mm against 100 mm; 50 and 20 mm
  not run: RAM limit); documented in this entry; committed locally (see git log)
- Phase B (rib placement search): not started
- Git: Stage 1 (1f64068), Stage 2 (7b67903) and Stage 3 (8e27d78) committed locally; Stage 4 and Stage 5 committed locally (see git log);
  nothing pushed
