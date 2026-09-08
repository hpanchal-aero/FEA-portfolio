# Stage 3 (Shoulder Fillet) — Baseline Sanity Check: CLOSED

Status: FROZEN. This record documents the verified baseline mesh/
solve for Stage 3. Geometry, mesh, sizing field, material, BCs,
loading, and solver setup for this baseline are not to be modified.
Any further refinement (convergence sweep) must be recorded as new,
separate mesh levels — this baseline itself stays fixed.

## Geometry (locked spec, unchanged)
- 60 x 40 mm cantilever
- Root section (x=0-30mm): depth D = 6 mm
- Tip section (x=30-60mm): depth d = 4 mm
- Symmetric opposite shoulder fillets, r = 0.95 mm
  (originally targeted r=1.0mm / r/d=0.25 exactly; reduced to 0.95mm
  after r=1.0mm failed the OCC fillet solver (h=r degeneracy) and
  r=0.99mm produced a mesh-quality defect from a 0.01mm sliver face -
  see build_fillet_case.py docstring for full construction history)
- Transition centered at x = 30 mm

## Baseline mesh (this specific level - the "fine" level, 0.075mm)
- 261,724 nodes / 156,219 C3D10 elements
- 784,701 equations (direct spooles solver)
- Max aspect ratio: 12.89
- Mean / median aspect ratio: 1.96 / 1.82
- Zero inverted (nonpositive-volume) elements
- Zero elements with AR > 20

## Loading and boundary conditions
- Material: Al 7075-T6, E=71,700 MPa, nu=0.33
- Root (x=0): full encastre, `*BOUNDARY root, 1, 3`
- Applied load: Fz = -235.44 N at reference node, via `*DISTRIBUTING`
  coupling to the 62-face tip surface at x=60

## Equilibrium check — PASSED
- Root Sum Fz = +235.4397 N (applied -235.44 N; opposite sign, equal
  magnitude, as required)
- Relative equilibrium error: ~0.0001%
- Root Sum Fx = +0.0002 N, Sum Fy = -0.0002 N (both ~0.0001% of
  applied load - no spurious lateral/axial reaction)

## Displacement
- Peak |U| = 0.4120 mm (node at tip edge, x=60/y=20/z=2)
- Dominant component Uz = -0.4112 mm; small Ux=0.0253mm (Poisson/
  bending coupling); Uy ~ 0
- Physically sane: smaller than Stage 1's ~1.05mm tip deflection
  under the same load/span, consistent with this geometry's greater
  average thickness (6mm root vs. Stage 1's uniform 4mm)

## Stress vs. Peterson reference
- sigma_nom (locked, 6M/(t d^2) at x=30) = 66.22 MPa
- Peterson Kt (r/d=0.25 target -> 0.2375 actual, D/d=1.5, h/r target
  1.0 -> 1.0526 actual, given r=0.95mm not the originally locked
  r=1.0mm) = 1.429 (unadjusted reference value, per the original
  specification gate - not recomputed at the as-built r/d)
- Predicted peak: 94.63 MPa

| | sigma_xx | Observed Kt | vs. Peterson (94.63 MPa) |
|---|---|---|---|
| Top fillet | +95.99 MPa | 1.450 | +1.44% |
| Bottom fillet | -94.98 MPa | 1.434 | +0.37% |

## Sanity check status: 8/8 CLOSED
1. Element type/count — C3D10, 156,219 — CLOSED
2. Node count — 261,724 — CLOSED
3. Material properties — Al 7075-T6 confirmed in deck — CLOSED
4. Root boundary condition — full encastre on 157-node root NSET — CLOSED
5. Load magnitude/direction/location — 235.44N, -z, distributed
   coupling — CLOSED
6. Reaction force / equilibrium — CLOSED (this record; originally
   flagged as a gap, closed via `analysis_fillet_rf.inp` rerun)
7. Displacement/stress scale — sane, consistent with Stage 1 — CLOSED
8. Rigid-body-motion check — no singular-pivot solver error; now also
   confirmed explicitly via exact equilibrium — CLOSED

## OPEN ITEM — not resolved, not a sanity-check failure
The mesh-convergence trend between the medium (0.15mm, baseline
solve) and this fine (0.075mm) level showed:
- Bottom-surface sigma_xx: 94.97 -> 94.98 MPa (+0.01%) — converged
- Top-surface sigma_xx: 93.82 -> 95.99 MPa (+2.31%) — NOT yet under
  the project's 2% convergence criterion, and still trending upward
  in the same direction as the coarse->medium step

This is an OPEN convergence item, explicitly not resolved by this
baseline verification and explicitly not a sanity-check failure —
the baseline itself is physically and numerically sound (equilibrium,
displacement scale, and Peterson agreement are all excellent). Two
attempted further refinements (0.06mm, 0.05mm) were mesh-quality-
verified but rejected as solver candidates on resource grounds
(extrapolated ~1.2-1.7M equations, too close to or above the
project's documented ~1.18M-equation direct-solver OOM point). A
mesh-growth diagnostic attributed ~70-73% of the excess growth to
the immediate fillet-surface refinement band itself and ~27-29% to
the fixed 1.5mm transition width, with far-field essentially
unaffected — recorded for reference if further refinement is
attempted later, but no further mesh work has been done as of this
record.

**This item remains open for a future session/decision:** whether to
accept the fine-level top-surface result with the 2.3% residual
documented as an honest limitation, or to pursue a further-refined
mesh (e.g. via a reduced TRANSITION_WIDTH, per the diagnostic
finding) before calling Stage 3's convergence study complete.
