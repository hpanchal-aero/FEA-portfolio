# Project 03 — UAV Wing Structural Analysis

A staged FEA study of a tapered UAV wing structure under combined
self-weight and aerodynamic loading, built with CalculiX, Netgen, and
Python. Full development history lives in
[`development-log.md`](development-log.md); this file is the summary.

## Research Question

How does rib configuration (spacing, then non-uniform placement at a
fixed rib count) affect a tapered UAV wing's tip deflection, skin
buckling resistance, and structural margin under a combined
self-weight + elliptical aerodynamic load, built up from a verified
spar-to-skin structural model with independently selected skin
thickness?

## Status

**Stages 1 (spar alone) and 2 (spar + ribs) complete and verified.**
Stages 3-6 not started.

## Approach

**Geometry:** NACA 0012 airfoil, tapered planform (root chord 300 mm,
tip chord 150 mm), 1 m half-span.

**Material:** Al 2024-T3 (E = 73,100 MPa, nu = 0.33, rho = 2780 kg/m^3
— standard literature values, not yet cited to a primary source).

**Load cases:** self-weight (gravity) and an elliptical spanwise
aerodynamic load, applied together as one combined/superposed static
case. Constraint: von Mises stress <= yield / 1.5. A linear eigenvalue
buckling analysis will be added once skin exists (Stage 3+), since
ribs' main structural role is resisting skin buckling.

**Build sequence** (progressive verification, Project-01-style):
1. **Spar alone — complete.** Tapered I-beam spar (root height 28.8 mm,
   tip height 14.4 mm, 30% chord position), verified against tapered
   Euler-Bernoulli beam theory.
2. **+ ribs at a placeholder spacing — complete.** Six true-airfoil ribs
   (1.5 mm thick, 200 mm spacing) fused to the spar.
3. + skin (full airfoil, skin structurally participating) — skin
   thickness selection study runs here — not started.
4. Load application refinement (distributing the elliptical load over
   the completed skin) — not started.
5. Rib spacing study (100/50/20 mm spacings) — not started.
6. Rib placement study (non-uniform positions at the winning rib
   count, screened for self-weight+aero structural margin) — not
   started.

Spar position and dimensions are held fixed for the entire project (a
documented scope limitation); only skin thickness and rib
configuration are varied.

## Findings — Stage 1 (spar alone)

Isolated tapered I-beam spar (30% chord, root height 28.8 mm tapering
to 14.4 mm at the tip, 2 mm web, 8x3 mm caps), root fixed, under a
**placeholder** combined load: self-weight + a 24.5 N total elliptical
lift (assuming a 5 kg aircraft at load factor n=1 — not a sourced
design requirement, to be replaced in Stage 4).

| | 2.0 mm mesh | 1.5 mm mesh |
|---|---|---|
| Tip deflection vs. tapered beam theory | -0.16% | -0.21% |
| Flange bending stress vs. beam theory | within 0.14% | within 0.22% |
| Peak von Mises (root corner) | 16.60 MPa | 18.24 MPa |

Tip deflection (3.84 mm) and flange stress agree with beam theory to
better than 0.25% and are converged between the two mesh levels. Peak
von Mises (at the root flange-web corner, a geometric discontinuity)
shows a real but modest mesh sensitivity (+9.9% from 2.0 to 1.5 mm),
consistent with the discontinuity-driven peak-stress behaviour
documented in Project 02; it is not treated as a design value. Both
values remain far below the ~290-345 MPa yield of Al 2024-T3. A
prismatic-root sanity bound (4.19 mm) is only a same-order check: its
simplifications partly cancel (see the development log).

## Findings — Stage 2 (spar + ribs)

Six NACA 0012 airfoil-shaped ribs (1.5 mm thick, at z = 0, 200, 400,
600, 800, 1000 mm) fused to the Stage 1 spar, same placeholder load.

| | 2.0 mm mesh | 2.5 mm mesh |
|---|---|---|
| Spar-axis tip deflection | 3.679 mm | 3.679 mm |
| vs. beam reference (ribs as point weights) | -0.12% | -0.11% |
| Flange bending stress vs. beam theory | within 0.06% | within 0.06% |
| Peak von Mises (root corner, unconverged) | 16.96 MPa | 16.16 MPa |

![Deformed shape, Uy](figures/contours/stage2_h2p0_uy_side.png)

*Uy, deformation exaggerated x50, black = undeformed.*

![von Mises](figures/contours/stage2_h2p0_vonmises.png)

*von Mises, colour scale capped at the 99th percentile (11.6 MPa).*

**The ribs add weight, not bending stiffness.** Tip deflection fell
4.2% relative to the spar alone (3.841 to 3.679 mm), and the ribs'
added self-weight alone predicts a 4.25% fall: the FEA matches a beam
model in which the ribs are point weights with zero bending stiffness
to within about 0.1%. Von Mises in the ribs is near zero. The ribs cost
+108 g (+49% mass) and, being off the spar axis, add a small twist
(about 0.1 degrees at the tip). Their structural value has to come from
the skin (buckling resistance, shape retention), which is what Stage 3
onward examines.

The path to this result included a load-application error (the
nodal load distribution followed the locally refined mesh around the
ribs, shifting load toward the root by 4% on the finer mesh) and an
earlier, incorrect attribution of the deflection change to rib
stiffness. Both are documented in the development log.

## Limitations

- The applied load (24.5 N total lift) is an explicit placeholder, not
  a sourced aircraft weight or lift requirement — to be replaced in
  Stage 4 once the full structure exists and a real load case is
  defined.
- The load is applied directly to the spar's top-cap outer face (a
  Stage 1-2 approximation, since no skin exists yet to distribute
  aerodynamic pressure over) — this will change once skin is added.
- Peak von Mises stress at the root is not mesh-converged and is not
  used as a design value.
- Convergence is demonstrated between two mesh levels; a finer level
  was excluded on memory grounds.
- The ribs' twist effect is not checked against torsion theory.
- The finding that ribs add no bending stiffness applies to this load
  case (spanwise bending, no skin).
- Al 2024-T3 material properties (E, yield, density) are standard
  literature values, not yet verified against a primary/certified
  source.
- Spar position (30% chord) and dimensions are fixed for the entire
  project — not varied or optimized (a deliberate scope decision, see
  development log).

## Conclusion

Stages 1-2 establish a verified FEA pipeline for this project: the
Netgen/CalculiX/Python toolchain, a tapered-beam numerical reference
for deflection and bending stress, and mesh/solve/extract scripts that
match that reference to within 0.25% on deflection and stress. Adding
ribs to the spar showed they contribute weight and twist but no bending
stiffness, which sets up the skin stages: the ribs' justification has
to be resisting skin buckling, and that is what the eigenvalue analysis
in Stage 3 will test. A real density-unit bug (1000x error, present
but never triggered in Projects 01-02) was found and fixed in Stage 1;
see the development log.
