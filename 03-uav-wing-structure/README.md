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

**Stage 1 (spar alone) complete and verified.** Stages 2-6 not started.

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
2. + ribs at a placeholder spacing (structure only) — not started.
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
| Equilibrium (root reaction vs. applied) | -0.14% | -0.09% |
| Tip deflection vs. tapered beam theory | -0.16% | -0.21% |
| Peak von Mises | 16.60 MPa | 18.24 MPa |

Tip deflection (3.84 mm) and equilibrium are essentially converged
between the two mesh levels. Peak von Mises (at the root flange-web
corner, a geometric discontinuity) shows a real but modest mesh
sensitivity (+9.9% from 2.0 to 1.5 mm) — consistent with the same
discontinuity-driven, slower-converging peak-stress pattern already
documented in Project 02; both values remain far below the ~290-345
MPa yield of Al 2024-T3. A conservative prismatic-root sanity bound
(4.19 mm, ignoring taper and self-weight) correctly overestimates
stiffness relative to the actual tapered result, as expected.

## Limitations

- The applied load (24.5 N total lift) is an explicit placeholder, not
  a sourced aircraft weight or lift requirement — to be replaced in
  Stage 4 once the full structure exists and a real load case is
  defined.
- The load is applied directly to the spar's top-cap outer face (a
  Stage-1 approximation, since no skin exists yet to distribute
  aerodynamic pressure over) — this will change once skin is added.
- Peak von Mises stress at the root is not mesh-converged (see
  Findings above) and is not used as a design value at this stage.
- Al 2024-T3 material properties (E, yield, density) are standard
  literature values, not yet verified against a primary/certified
  source.
- Spar position (30% chord) and dimensions are fixed for the entire
  project — not varied or optimized (a deliberate scope decision, see
  development log).

## Conclusion

Stage 1 establishes a verified FEA pipeline for this project: the
Netgen/CalculiX/Python toolchain, the tapered-beam numerical
verification method, and the mesh/solve/extract scripts all produce
results matching closed-form-adjacent theory to within 0.2% on the
governing (deflection, equilibrium) quantities. A real density-unit
bug (1000x error, present but never triggered in Projects 01-02) was
found and fixed here — see development log. Ready to proceed to Stage
2 (adding ribs).
