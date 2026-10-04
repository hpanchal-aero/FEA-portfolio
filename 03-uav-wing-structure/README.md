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

**Stages 1 (spar alone), 2 (spar + ribs) and 3 (spar + ribs + skin,
skin thickness selected) complete.** Skin thickness is locked at
**1.0 mm**. Stages 4-6 not started.

## Approach

**Geometry:** NACA 0012 airfoil, tapered planform (root chord 300 mm,
tip chord 150 mm), 1 m half-span, spar at 30% chord, six ribs at
200 mm spacing (1.5 mm thick, true airfoil outline).

**Material:** Al 2024-T3 for spar, ribs and skin (E = 73,100 MPa,
nu = 0.33, rho = 2780 kg/m^3 — standard literature values, not yet
cited to a primary source).

**Spar:** Stages 1-2 used a spar whose depth was 80% of the local
airfoil thickness. Stage 3 found that this leaves a 3.6 mm (root) /
1.8 mm (tip) gap between the spar caps and the airfoil outline, so no
skin could touch the spar. The spar was therefore redesigned ("spar
v2": depth equal to the airfoil depth, same 2 mm web and 8x3 mm caps,
+30 g). Stage 1-2 results describe the old spar and are kept as the
verified baseline of the pipeline.

**Load cases:**
- Stages 1-2: self-weight + a 24.5 N elliptical lift (5 kg aircraft,
  n = 1) — an explicit placeholder.
- Stage 3: **limit load**, load factor n = 3.8, an assumed 20 kg
  aircraft, so a half-wing elliptical lift of 372.78 N, with
  self-weight x 3.8 on solid and skin. The load factor follows the 3.8
  cap of the transport-category formula (14 CFR 25.337(b)); the primary
  STANAG 4671 / CS-23 text was not retrieved. **The aircraft mass is an
  assumption and not from any standard.**
- Requirements at limit load: von Mises <= yield / 1.5 (yield assumed
  290 MPa, uncited) and linear buckling factor >= 1.5.

**Skin model:** closed airfoil skin of S6 shell elements (0.5 / 1.0 /
1.5 mm candidates) sharing nodes with the solid ribs and spar caps
along their contact lines (conformal Netgen mesh, no constraint
cards). CalculiX's `*TIE` does not support shell elements.

**Build sequence:**
1. Spar alone — complete.
2. + ribs at 200 mm spacing — complete.
3. **+ skin, skin thickness selection and buckling — complete
   (t = 1.0 mm locked).**
4. Load application refinement (distributing the elliptical load over
   the skin instead of the spar cap) — not started.
5. Rib spacing study (100 / 50 / 20 mm) — not started.
6. Rib placement study (non-uniform positions at a fixed count) — not
   started.

Spar position and dimensions are held fixed for the project (a
documented scope limitation); only skin thickness and rib
configuration vary.

## Findings — Stage 1 (spar alone, old spar, 24.5 N placeholder load)

| | 2.0 mm mesh | 1.5 mm mesh |
|---|---|---|
| Tip deflection vs. tapered beam theory | -0.16% | -0.21% |
| Flange bending stress vs. beam theory | within 0.14% | within 0.22% |
| Peak von Mises (root corner) | 16.60 MPa | 18.24 MPa |

Tip deflection (3.84 mm) and flange stress agree with beam theory to
better than 0.25% and are converged between the two mesh levels. The
peak von Mises (root flange-web corner) is not converged and is not a
design value. A prismatic-root sanity estimate (4.19 mm) is a
same-order check only: its simplifications partly cancel (see the
development log).

## Findings — Stage 2 (spar + ribs, old spar, 24.5 N placeholder load)

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

**The ribs add weight, not bending stiffness:** the tip deflection fell
4.2% (3.841 to 3.679 mm), and the ribs' added self-weight alone
predicts 4.25%. The ribs cost +108 g (+49% mass) and add a small twist.
Their structural value has to come from the skin.

## Findings — Stage 3 (spar v2 + ribs + skin, limit load)

![Mesh overview](figures/stage3_mesh_overview.png)

*Production mesh (t = 1.0 mm): free skin (blue), skin sharing nodes with
the solid (orange: rib rims and spar-cap strips), other solid faces
(grey).*

![Root close-up](figures/stage3_mesh_root.png)

**Skin thickness study** (limit load n = 3.8, lift 372.78 N per
half-wing, 200 mm rib spacing, aircraft mass 20 kg assumed, yield
290 MPa assumed):

| Skin t | Skin mass | Tip deflection | Buckling factor (>= 1.5) | Design skin stress* | Result |
|---|---|---|---|---|---|
| 0.5 mm | 0.64 kg | 8.36 mm | 0.91 (6 mm skin mesh: 0.93) | 39.4 MPa | **fails** buckling |
| **1.0 mm** | 1.28 kg | 4.20 mm | **3.14** (6 mm: 3.12) | 34.2 MPa | **passes; locked** |
| 1.5 mm | 1.92 kg | 2.52 mm | 8.61 | 25.3 MPa | passes, 50% heavier |

*Maximum von Mises over skin nodes at least 3 mm from the trailing edge
and more than 10 mm from the fixed root (see the singularity finding
below); factors of safety against the assumed yield are about 7.4 / 8.5
/ 11.5.

![Thickness study](figures/stage3_thickness_summary.png)

![Critical buckling mode](figures/stage3_buckling_mode1.png)

*Critical mode at t = 1.0 mm (factor 3.14): a single local buckle of the
upper skin aft of the spar in the root bay. All six computed modes, at
every thickness, are upper-skin, aft-of-spar, out-of-plane buckles in
the two root bays, none at the edges or at rib planes.*

**Verification and convergence (what was and was not shown):**
- Load application: results at the n = 1 and n = 3.8 load cases
  superpose exactly (lift x15.215 minus gravity x3.8 reproduces the
  limit-load tip deflection to four digits), and the effect of solid and
  skin gravity matches the beam model to 0.5% and 0.1%.
- Tip deflection vs. a tapered beam with a fully effective skin ring
  (an approximate reference, not a bound): +0.82% / -0.35% / -3.07% at
  0.5 / 1.0 / 1.5 mm; with a nearly absent skin (0.05 mm control) the
  model matches the spar-only beam to +0.15%.
- Mesh: solid 2.5 -> 4 mm changes the tip deflection by -0.02%; skin
  8 -> 6 mm changes the t = 1.0 mm buckling factor by -0.5% (and the
  0.5 mm factor by +2.2%), and the skin stress percentiles by under 2%.
- Buckling plausibility: the critical stress of the first mode (32.5 MPa)
  agrees with simply supported flat-plate theory for the buckle's own
  width (31.6 MPa) to 2.7% — a plausibility check, not a verification.

**Trailing-edge stress peak is a singularity of the model.** The
largest skin stress (238 MPa at t = 1.0 mm) sits at the trailing-edge
vertex at the ribs. Local refinement at one rib (median edge 2.1 to
0.28 mm) made the peak grow without bound (232 -> 536 MPa skin,
83 -> 473 MPa solid), while the stress farther than 3 mm from the
trailing edge converged (3-6 mm ring: 29.3 -> 35.7 MPa; 6-15 mm:
15.1 -> 19.3 MPa) and unrefined control ribs did not move. This peak is
excluded from the stress requirement; the true trailing-edge stress of a
real finite-thickness trailing edge is not predicted.

![Trailing-edge refinement](figures/stage3_te_refinement.png)

**Sensitivity to the assumed aircraft mass.** The analysis is linear, so
the buckling factor scales roughly as 1 / load: 0.5 mm would pass for an
aircraft of about 12 kg or less, 1.0 mm up to about 40 kg. The locked
thickness is valid for the 20 kg assumption only.

**Unresolved:** away from the root the skin stress and bending curvature
in the FEA are lower than the beam reference (stress ratio at
z = 700 mm: 0.88 / 0.74 / 0.58 at 0.5 / 1.0 / 1.5 mm; the tip deflection
agreement within 3% is partly a compensation with extra inboard
flexibility). A skin-free control agrees, so the cause is skin-related
and not found. Reaction-force output is not interpretable for CalculiX
expanded shells, so global equilibrium is checked indirectly (see the
development log).

## Limitations

- **Unsourced inputs:** aircraft mass (20 kg, my assumption), yield
  (290 MPa, uncited), and the load factor's primary source (STANAG 4671
  / CS-23 not retrieved). The thickness lock depends on the first two.
- The elliptical load is still applied to the spar's top-cap face
  (Stage 4 will distribute it over the skin).
- Shell rotations at shared nodes are not coupled to the solid (solid
  nodes have no rotational DOF): the skin is effectively hinged to the
  caps and rib rims. Expected to be conservative for buckling; not
  quantified.
- Linear eigenvalue buckling only: no imperfections or post-buckling.
- Trailing-edge geometry is idealized (asymmetric polyline closure, sharp
  junction at the ribs); the trailing-edge stress is a mesh-dependent
  singularity and is excluded. Refinement was done at one rib.
- Four thin wedge-shaped skin pieces (0.07% of the skin area) force
  triangle angles down to 0.18 degrees; accepted, effect not isolated.
- The skin-vs-beam stress and curvature discrepancy above is unexplained.
- Global equilibrium by reaction forces could not be evaluated for the
  shell model.
- Convergence is shown between two mesh levels per mesh parameter.
- Stages 1-2 results describe the old (80%-depth) spar and the 24.5 N
  placeholder load.
- Al 2024-T3 properties are literature values, not verified against a
  primary source. Spar position and dimensions are fixed project-wide.

## Conclusion

Stages 1-3 establish a verified solid-and-shell FEA pipeline and a
skin-thickness selection under a defined limit-load case. With 200 mm
rib spacing the lightest of the three candidate skins that meets both
the buckling requirement (factor 3.1 at 1.0 mm against 1.5) and the
stress requirement (factor of safety about 8 outside the trailing-edge
singularity) is 1.0 mm; 0.5 mm buckles below limit load. Skin buckling
is local (the aft panel between spar and trailing edge, about 93 mm
wide), so ribs closer than about that width, not the 200 mm baseline
spacing, are what the rib studies of Stages 5-6 will probe. The result
is conditional on the stated assumptions and carries the documented
unresolved items. See the development log for the full record,
including every error and retraction.
