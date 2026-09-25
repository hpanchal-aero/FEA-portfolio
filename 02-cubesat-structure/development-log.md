# Project 02 — 1U CubeSat Primary Structure — Development Log

## Specification (locked)

**Research question**: How does CubeSat structural configuration
affect mass, stiffness, and stress margin under representative
launch loading?

**Material**: Al 6061-T6 (E=68,900 MPa, nu=0.33, rho=2700 kg/m^3,
yield ~276 MPa). These are commonly quoted typical values; no source
has been cited yet. The 276 MPa yield is a typical value, and design
allowables (e.g. MMPDS) are typically lower. To be sourced before any
final margin claim.

**Loading**: 12g quasi-static axial launch acceleration. Explicitly a
**preliminary design assumption, not a sourced launch-vehicle
requirement** — same status as Project 01's 12g load assumption,
carried forward for methodological consistency across the portfolio.

**Structural configurations to compare**:
1. Configuration A — solid-panel frame (flat panels, corner rails, no
   lightening features)
2. Configuration B — ribbed/lightened frame (same envelope, panels
   with lightening pockets or rib pattern)

**Constraints**:
- Max von Mises stress <= yield/FoS, FoS=1.5 -> allowable ~184 MPa
- Mass checked against the 1U CDS budget (1.33 kg) as a design
  target, not necessarily a hard pass/fail constraint. 1.33 kg is
  widely quoted in secondary sources; it has not been verified against
  the primary CDS text (Rev 14.1 restructured its requirements into
  guidelines). The Configuration A structure mass (257.31 g) is the
  primary structure only, compared against a whole-spacecraft budget.
- No displacement limit set — CubeSat panel structures are typically
  stiffness/frequency-governed rather than displacement-governed;
  open to revision if the FEA results suggest otherwise

**Explicitly out of scope**: modal/frequency response, random
vibration loading — deferred to Project 07 (CubeSat modal analysis)
per the original portfolio sequencing.

**Verification strategy** (amended after sub-stage 2a; the original
text said "single panel as plate bending"): no single closed-form
reference exists for the full assembled structure, unlike Project 01's
progressive bracket stages. Individual members are verified against
uniaxial (axial-bar) theory, sigma = P/A and delta = PL/(AE), because
the load is in-plane and no transverse load exists, so plate bending
theory does not apply. Closed-form plate-buckling (panel) and Euler
(rail) checks are added for the isolated members; no eigenvalue
buckling FEA is performed and no buckling analysis of the assembly is
done. The full assembly is treated as an internal-consistency study
(equilibrium checks, mesh levels, a boundary-condition diagnostic)
rather than claimed as closed-form-verified — the same honest framing
Project 01 used for its own Stage 5 full-bracket integration.

**Toolchain**: geometry is built with the Gmsh OCC kernel (geometry
only; no Gmsh meshing is used), meshing with Netgen (straight-sided
C3D10 by manual midside-node promotion), CalculiX 2.23 (ccx223 conda
environment), Python (NumPy, PyVista, Matplotlib). Same locked
toolchain as Project 01's Stage 5 onward.

**Directory structure** (actual): `scripts/`, `simulation/`, `mesh/`
(geometry and mesh inputs in `mesh/<case>_study/`), `results/`,
`figures/`. A `geometry/` directory exists but is unused.

## Status Log

- [x] Specification defined and approved
- [x] Baseline geometry (Configuration A) — validated
- [x] Baseline mesh
- [x] Baseline solve + sanity check (equilibrium, uniaxial comparison)
- [~] Mesh convergence study — two levels (3.0, 2.0 mm); means stable, peak not converged (BC-sensitive); further refinement deliberately skipped
- [x] Feature-level analytical verification (panel 2a, rail 2b)
- [ ] Configuration B (ribbed/lightened) — begins after Configuration A is finished; early geometry candidates exist (see below)
- [ ] Comparative analysis (A vs. B: mass, stress, stiffness)
- [~] Documentation and figures — this update covers Configuration A
- [~] Git: Configuration A stage committed locally; push deferred until the project is complete

## Load-sharing assumption

1.33 kg x 12 g x 9.81 m/s^2 = 156.6 N, applied as uniform pressure
0.1643 MPa over the 953 mm^2 top face (0.1643 x 953 = 156.578 N).
Rails and panels share a rigid base and a common 100 mm length, so
equal axial strain gives equal axial stress. This is a bounding
assumption, not a sourced launch load. Applying the same pressure to
the whole combined top face follows from that assumption and was
stated explicitly in the input-file header.

## Configuration A geometry

100 mm envelope (x, y in [-50, 50], z in [0, 100]); four
8.5 x 8.5 x 100 mm rails at the corners with outer faces flush with
the envelope; four 2 mm x 83 mm x 100 mm panels between the rails,
outer faces flush; top and bottom panels omitted; all eight bodies
fused into one solid (real rails are bolted to panels at discrete
fasteners; this idealization cannot capture bolted-joint contact or
local joint stress).

Validation (build script `scripts/build_config_a_base_frame.py`, all
reproduced in the regenerated run): rail and panel pre-fuse dimension
checks (the script aborts on any mismatch, and the run completed);
rail/panel overlap 0.000000 mm^3; fuse gives exactly 1 volume; volume
95,300.0000 mm^3 against the analytic 95,300 (rails 28,900 + panels
66,400); bounding box exact; 18 surfaces, no slivers; z = 50 section
area 953.0000 mm^2 against 953; model cleanliness 1 volume / 18
surfaces; mass 257.31 g at 2700 kg/m^3. Outputs:
`mesh/config_a_study/config_a_base_frame.brep` and
`config_a_clean.step`.

## Sub-stages 2a and 2b (feature-level verification)

Isolated 83 x 2 x 100 mm panel (2a) and 8.5 x 8.5 x 100 mm rail (2b);
base fixed, uniform 0.1643 MPa on the top face, 1.5 mm C3D10 mesh.
Reference: uniaxial theory (delta = 0.00023846 mm), not plate bending.
Results were re-extracted from the stored `.frd` files by
`scripts/extract_isolated_member.py` on 2026-09-21. Mean von Mises uses
the slab 25 < z < 75 mm.

| | 2a panel | 2b rail |
|---|---|---|
| Nodes / elements / equations | 81,735 / 47,888 / 245,205 | 18,755 / 11,398 / 56,265 |
| Applied force | 27.2738 N | 11.8707 N |
| Bottom reaction (Z) | 27.273802 N (+0.0000%) | 11.870676 N (+0.0000%) |
| Top reaction (Z) | -27.273829 N | -11.870680 N |
| Top Uz mean vs -0.00023846 mm | -0.00023537 mm (magnitude 1.29% below) | -0.00023767 mm (0.33% below) |
| Mean von Mises, slab | 0.163828 MPa (-0.288%) | 0.164298 MPa (-0.001%) |
| Peak von Mises | 0.345026 MPa (2.11x mean) | 0.259284 MPa (1.58x mean) |
| Peak location | panel edges/corners at z = 0 (top 10 all z < 3 mm) | rail cross-section corners at z = 0 (top 10 all z < 3 mm) |
| Closed-form check | plate, k = 4, b = 83, t = 2: sigma_cr 147.70 MPa, P_cr 24,517.89 N, applied/critical 0.0011 | Euler, K = 2: I = 435.0052 mm^4, P_cr 7,395.26 N, applied/critical 0.0016 |

The peaks are located at the fixed base (Saint-Venant boundary
effect), checked by listing the highest-stress nodes and their
coordinates. The buckling values are closed-form checks on the isolated
members (simply supported plate assumed; fixed-free column assumed);
they say nothing about buckling of the assembly.

**Superseded earlier figures.** Earlier session values gave 2a/2b
displacement differences of 1.10% / 0.14% and mean von Mises of
0.161993 / 0.163368 MPa (about -1.4% / -0.55%). The displacement
percentages match a comparison against the rounded theory value
0.000238; against the exact 0.00023846 they are 1.29% / 0.33%. The
earlier mean-stress region definition is not known. Both are replaced
by the re-extracted values above. Node and element counts, equilibrium,
top displacement and peak stress otherwise reproduce the earlier
figures.

**Sign convention note.** `extract_isolated_member.py` reports the
top-displacement difference as (FEA + theory)/theory, positive when the
FEA displacement magnitude is smaller than theory, while
`extract_config_a.py` reports (FEA - (-theory))/(-theory), negative in
the same situation. The two describe the same physical effect (the FEA
structure is slightly stiffer than theory). All tables in this project
state the effect as "magnitude below theory".

## Full-assembly meshing: failures and resolution

Three meshing attempts on the fused assembly failed with
self-intersection warnings and "Meshing failed" (uniform 3.0 mm;
uniform 0.8 mm; face-specific sizing). The first two were diagnosed as
a thin-wall sizing problem (3.0 mm exceeding the 2 mm panel
thickness), and the third assumed panel and rail faces stay
identifiable after the boolean fuse.

**Root cause, found later: those diagnoses were wrong.** The saved
BRep was a compound holding the valid 18-face solid plus one stray
953 mm^2 face at z = 50, left in the model by the sectional-area check
(the intersect operation deleted its inputs but kept its result).
Netgen meshed that face as part of the geometry. After extracting the
single solid, uniform 3.0, 2.0 and 1.5 mm all meshed with zero
intersection warnings, so the thin wall was not the blocker.

The observation that OCC merges coplanar touching faces (the four
outer sides came out as single 100 x 100 mm faces, so panel and rail
faces cannot be classified by bounding box) is true, but it was not
the cause of the failures. The face-classification approach solved a
problem that did not exist.

Fix: the build script now removes the section residue, asserts 1
volume and 18 surfaces before saving, and writes both BRep and STEP.
The regenerated geometry meshes identically (31,013 / 83,854 /
155,900 tets at 3.0 / 2.0 / 1.5 mm). The 1.5 mm mesh was generated as
a feasibility test only and was not solved. The defective BRep is kept
locally as `config_a_base_frame_DEFECTIVE.brep` and is not intended for
commit. Two obsolete scripts were deleted: `mesh_and_solve_config_a.py`
(the face-classification attempt) and `build_and_solve_2a_panel.py`
(an earlier Gmsh-based 2a attempt, superseded by the Netgen script).

## Extraction error: reaction sum

The first extraction reported a Z equilibrium error of -100%. Cause: an
error in the extraction script. It summed the CalculiX RF output over
all nodes, which is zero by construction (the loaded top nodes carry
the opposite sign of the base reaction). The reaction of the fixed base
is the sum over the constrained bottom nodes. The corrected script
reports bottom, top and all-node sums; the solve was not at fault.

## Convergence stage results (fixed base)

| | 3.0 mm | 2.0 mm |
|---|---|---|
| Tets / C3D10 nodes / equations | 31,013 / 62,076 / 186,228 | 83,854 / 158,069 / 474,207 |
| Min / mean tet quality | 0.2945 / 0.8811 | 0.2407 / 0.8827 |
| Bottom Z reaction vs 156.5779 N | +0.0000% (156.577904 N) | -0.0000% (156.577890 N) |
| Top reaction (Z) | -156.578024 N | -156.577888 N |
| Top Uz mean vs -0.00023846 mm | -0.00023439 mm (magnitude 1.71% below) | -0.00023453 mm (1.65% below) |
| Mean von Mises, panel region | 0.164240 MPa | 0.164240 MPa |
| Mean von Mises, rail region | 0.160468 MPa (-2.33%) | 0.160611 MPa (-2.25%) |
| Peak von Mises (location) | 0.368 MPa at (-50, -50, 0) | 0.422 MPa at (-50, 50, 0) |
| Peak / slab mean | 2.25x | 2.59x |
| Peak RSS / CalculiX wall time | 0.73 GB / 13 s | 2.18 GB / 36 s |

Tet quality is q = 12(3V)^(2/3) / sum(l^2), which is 1 for a regular
tet; no tet was below 0.1. Peak RSS and wall time were read from the
terminal output during the runs and are not stored in a file. The peak
sits at the four outer base corners in both runs, with the ten highest
nodes within about 1.5 mm of the base.

The region means are the stable quantities. The "all nodes in slab"
mean moved from -0.42% to -0.61% against theory only because the rail
share of nodes in the slab rose from 16.8% to 26.1% with refinement
(a node-count-weighted mean), so it is not used.

Mean stress and displacement are stable between the two meshes. Two
levels show stability but give no convergence order. The peak rose 15%
and is reported as **not converged and BC-sensitive**. The 1.5 mm and
2.5 mm levels were deliberately not run.

## Roller-base diagnostic (3.0 mm, diagnostic only)

Base restraint changed to Uz = 0 on all bottom nodes plus the minimum
rigid-body constraints (node 1: Ux, Uy; node 5: Uy), via
`scripts/make_roller_bc_case.py`. Results: top Uz -0.00023846 mm,
+0.0002% against theory; all three region means and the maximum nodal
stress at 0.164298 MPa (peak = 1.00x mean); base reaction +0.0001%
(156.578018 N); maximum lateral displacement 7.87e-5 mm against
nu * epsilon * 100 mm = 7.87e-5 mm (valid here because node 1 pins
x = -50, so the far side moves across the full 100 mm width).

The 1.7% displacement offset, the 2.25x stress peak and the 15%
peak growth on refinement are therefore attributed to the lateral
restraint of the fixed base (this is a model effect, not a mesh error).
This run also confirms that the load, material, C3D10 promotion and
connectivity across the fused rail-panel junction reproduce a uniform
stress state exactly. The locked boundary condition remains the fixed
base.

Not explained: (a) why the rail region reads about 2% below the panel
region under the fixed base while the isolated rail (2b) shows no such
deficit; (b) the area-weighted isolated-member displacement offsets
(1.29% panel, 0.33% rail, about 1.0%) are below the assembly's
1.65-1.71%; the meshes differ (1.5 mm against 3.0 / 2.0 mm), so no
attribution is made. Neither was tested.

## Limitation: axial-only loading

Under the locked axial 12g case, stress in each member of a prismatic
configuration is set by total cross-section area and stiffness by
EA/L, which the roller run reproduced exactly. This load case cannot
distinguish arrangements of equal area in such a configuration, so for
Configuration A and prismatic variants the comparison is primarily a
mass-versus-area trade. The specification is unchanged. Configurations
with non-prismatic panels (pockets) are not covered by this statement,
because their stiffness and local stress depend on how material is
distributed along the load path.

**Proposed extension (not adopted):** lateral or combined
axial + lateral quasi-static loading, where the arrangement of panels
and rails affects stiffness and stress. Modal analysis remains
deferred to Project 07.

## Configuration B status

Not begun as a study; it starts after Configuration A is finished.
Early geometry candidates exist from 2026-09-19
(`scripts/build_config_b_panels.py`, `mesh/config_b_study/`,
`results/config_b_geometry_validation.csv`): isolated 83 x 100 x 2 mm
panels with a through-cut pocket and a 5 mm retained frame, four pocket
shapes (rectangle, circle, cross, 2 x 2 grid) at five removal levels
(10-30% of panel area) with equal removed area per level. All 20
passed geometry validation (removed area within 1e-12 %, no slivers).
This pocket-shape and removal-level design has not been adopted for the
study and will be revisited when Configuration B starts. Removal
levels are fractions of panel area; removing 30% of panel area removes
about 20.9% of the total structure volume (panels are 66,400 of
95,300 mm^3). None of these geometries has been meshed or solved, and
they should be re-checked as single solids in Netgen before use, given
the stray-face issue found in Configuration A.

## Open items and unverified sources

- Al 6061-T6 property values and the 276 MPa yield: source to be cited.
- CDS 1U mass (1.33 kg): not verified against the current CDS text.
- Rail-region mean-stress deficit in the assembly: unexplained.
- Convergence order for mean quantities: not established (two levels).
- Repository hygiene before commit: the parent `.gitignore` is
  Project 01's; Project 02 large `.inp` files, `_DEFECTIVE.brep` and
  the duplicate STEP need decisions.

## Reproducibility

```bash
cd 02-cubesat-structure
# Configuration A geometry (Gmsh OCC), validated; writes BRep + STEP
/usr/bin/python3 scripts/build_config_a_base_frame.py
# Mesh feasibility test (mesh generation only)
/usr/bin/python3 scripts/test_config_a_meshing.py
# Mesh + write CalculiX input for a uniform size, e.g. 3.0 mm or 2.0 mm
/usr/bin/python3 scripts/mesh_config_a.py 3.0
# Solve under an RSS-polling memory watchdog (8 GB limit)
/usr/bin/python3 scripts/run_ccx_watchdog.py simulation/config_a_h3p0 8
# Extract and check results (equilibrium, displacement, stress)
/usr/bin/python3 scripts/extract_config_a.py simulation/config_a_h3p0
# Roller-base diagnostic
/usr/bin/python3 scripts/make_roller_bc_case.py simulation/config_a_h3p0.inp simulation/config_a_h3p0_rollerBC.inp
/usr/bin/python3 scripts/run_ccx_watchdog.py simulation/config_a_h3p0_rollerBC 8
/usr/bin/python3 scripts/extract_config_a.py simulation/config_a_h3p0_rollerBC
# Mesh figures (PNG + STL of the outer surface)
/usr/bin/python3 scripts/view_config_a_mesh.py simulation/config_a_h3p0.inp
# Sub-stages 2a / 2b: mesh + .inp, solve, re-extract
/usr/bin/python3 scripts/mesh_and_solve_2a_panel.py
/usr/bin/python3 scripts/mesh_and_solve_2b_rail.py
/usr/bin/python3 scripts/run_ccx_watchdog.py simulation/stage2a_panel_solve 8
/usr/bin/python3 scripts/run_ccx_watchdog.py simulation/stage2b_rail_solve 8
/usr/bin/python3 scripts/extract_isolated_member.py a
/usr/bin/python3 scripts/extract_isolated_member.py b
```

Netgen scripts run under the system Python (`/usr/bin/python3`); the
solver is CalculiX 2.23 from the `ccx223` conda environment
(`/home/harsh/miniconda3/envs/ccx223/bin/ccx`). The 2a/2b solves were
originally run by direct `ccx -i` calls, and no wrapper for those runs
was saved. The watchdog commands for 2a/2b in this list are the
reproducible equivalent and have not themselves been re-run.

## Configuration B: design and load application

Same rail layout as Configuration A. Each of the 4 panels gets one
through-cut pocket, centered at z=50, with a 5 mm perimeter frame
retained (matching the isolated-panel design already generated at
project kickoff, `scripts/build_config_b_panels.py`). Four shapes
(rect, circle, cross, 2x2 grid) x three removal levels (10/20/30% of
panel area, chosen from the originally planned 5, per the decision to
start with 3 and add more only if trends needed it — they did not).

**Load application mechanism, verified before use.** Config A's
uniform-pressure `*DSLOAD` assumes equal stress across members, valid
only for prismatic geometry. For non-prismatic pocketed panels this
would silently impose an assumption rather than let it emerge, so
Config B instead ties every NTOP node's Uz to a reference node
(`*EQUATION`, off-structure at (0,0,150), Ux/Uy fixed to avoid
singularity) carrying the total 156.5779 N via `*CLOAD`
(`scripts/make_equation_load_case.py`, `scripts/mesh_config_b.py`).
Verified on Config A's 3.0mm mesh, both boundary conditions:

| | Roller base | Fixed base |
|---|---|---|
| Bottom reaction vs 156.5779 N | +0.0001% | +0.0000% |
| Top Uz vs theory | +0.0002% | -1.6500% (pressure case: -1.71%) |
| All region means | 0.164298 MPa (matches pressure case exactly) | panel +0.22%, rail -0.39% vs pressure case |
| Peak VM | 0.164298 MPa (1.00x mean) | 0.368391 MPa (pressure case: 0.368141) |

The roller case reproduces the pressure-loaded roller result exactly
(uniform stress emerges from tied displacement, not imposed). The
fixed-base case shows small (0.2-0.4%) load redistribution toward the
stiffer rails versus uniform pressure — physically explainable, well
inside the ~2% offset already attributed to the base BC. Mechanism
approved for Config B use.

## Configuration B: geometry

Build script `scripts/build_config_b_assembly.py`. Pocket dimension
formulas reused from `scripts/build_config_b_panels.py` (rect: solves
for w,h at fixed aspect ratio matching the panel's interior W_INT x
H_INT; circle: r = sqrt(A/pi); cross: two full-length arms solved for
arm width; grid: 4 separated square blocks, side = sqrt(A/4)).

**Sectional-check bug found and fixed.** The original SECTION_Z=10mm
(reused from Config A reasoning) assumed every shape's pocket stays
clear of z=10 at every level. True for rect/circle/grid (verified:
max z half-extent at 30% is 27.7/28.15/28.95mm, all < 40mm). **False
for cross**: its vertical arm always spans the full interior height
(H_INT/2 = 45mm half-extent) regardless of removal level, so z=10 sits
inside the arm at every level. First cross run (10%) failed the
sectional check (measured 910.9 vs expected 953.0 mm^2, -4.42%).
Fixed by moving SECTION_Z to 2mm (provably inside the always-uncut
5mm frame margin for every shape by construction) and adding an
explicit per-shape cutter z-extent assertion so this class of mistake
cannot recur silently. Cross re-run passed after the fix; rect,
circle, grid were confirmed safe at the original z=10 by the same
z-extent check (not re-run, since their z=10 result was already
correct for all levels used).

**All 12 geometries validated**: 1 solid each, cleanliness check
(total faces = solid boundary faces) passing for every case, per-panel
removed-area error 0.0000% at every level, total volume exact
(88,660 / 82,020 / 75,380 mm^3 at 10/20/30%, i.e. 95,300 -
4 x level% x 16,600), bounding box exact, section area 953.0000 mm^2
at z=2 for every case. Surface counts per shape (constant across
levels): rect 34, circle 22, cross 66, grid 82. Masses: 239.38 /
221.45 / 203.53 g (7.0% / 13.9% / 20.9% of Config A's 257.31 g removed
at 10/20/30% panel-area removal — not equal to the panel-area
percentage, since panels are only 66,400 of the 95,300 mm^3 total).

## Configuration B: meshing feasibility and mesh sizing

Feasibility test (mesh generation only, no `.inp`) at 3.0mm on all 12
`.step` files: all pass, zero intersection warnings, 21,590-28,313
tets. Solved at 3.0mm: all 12 pass, equilibrium <1e-4% throughout,
peak RSS not individually captured for this batch (logging mistake:
the summary command grepped "peak RSS" against the wrong file — that
string is only printed by run_ccx_watchdog.py to stdout, never written
into the CalculiX .ccx.log; fixed for later batches by piping full
stdout to a per-case log file). What is known: all 12 exited 0 with no
watchdog kill, so all stayed under the 8GB limit.

**A second logging mistake in the same batch**: the extraction
commands were run through `tail -12`, which truncated
`extract_config_a.py`'s Equilibrium and Overall-axial-stiffness
sections from the captured output for 10 of the 12 cases (only the
standalone cross-10% test had full output). Recovered by reading the
already-written `results/*_configB_results.json` files directly
(read-only, no re-solve needed) rather than re-running anything.

## Configuration B: mesh convergence

3.0mm -> 2.0mm, all 12: stiffness stable everywhere (0.03-0.14% shift).
Net-section stress: rect/circle stable (0.06-2.0% shift); cross/grid
not stable (cross 1.5-8.2%, grid 10.8-20.3%, both worsening with
removal level for grid specifically).

**1.5mm was attempted for cross/grid** (the two non-converged shapes)
and rejected before solving: feasibility test showed 53-69% of the
element ceiling and 57-69% of the equation ceiling — well above the
informal 25-30% comfort margin and in the same risk range as Config
A's projected-but-never-solved 1.5mm case, given the documented WSL
crash history. **1.8mm chosen instead** as a safer intermediate step
(34-40% element ceiling, matching Config A's already-solved 2.0mm case
in relative terms). Solved for all 6 (cross/grid x 3 levels): peak RSS
1.80-2.36 GB, all well under the 8GB limit, wall times 68-103s (versus
~20-30s at 2.0mm).

**Grid's growing 3.0->2.0mm percentage difference (10.8% to 20.3%
across levels) was initially read as a possible real or worsening
trend, or a flaw in the corner-exclusion metric.** The third data
point (1.8mm) resolved this: the 2.0->1.8mm step shrank sharply for
every cross/grid case (to 0.03-2.5%, an order of magnitude smaller
than the prior step), confirming genuine convergence rather than
divergence. This is recorded as an explicit caution: two mesh levels
were not enough to distinguish "converging slowly from a bad starting
point" from "diverging" for this metric, and a third level changed the
conclusion. Stiffness was already flat across all three levels
(0.01-0.16% total movement 3.0->1.8mm).

**Converged/reported mesh level per shape**: rect and circle at 2.0mm
(already stable there); cross and grid at 1.8mm (needed the third
level). This is a mixed-mesh-level comparison across shapes, stated
explicitly rather than left implicit, since "same mesh size" and
"converged mesh size" were different requirements for different
shapes here.

## Configuration B: net-section stress method

Not an integrated force/area value. Mean axial stress SZZ (signed,
not von Mises, since von Mises at a re-entrant corner is dominated by
local concentration) over panel-region corner nodes at z = 50 +/-
1.5mm, reported as a full-section mean and with the highest/lowest
decile (by SZZ value) excluded as a corner-singularity-exclusion
proxy. Implemented in `scripts/extract_config_b.py`.

**Limitation, observed directly on the cross-10% case**: the
excl-decile trim only produced a 1.5% gap from the full mean, far
smaller than expected for a shape with severe corner concentration.
Inspection of the peak-VM node list showed the five worst corners all
lie within the z = 50 +/- 1.5mm section band itself, meaning at this
severe a throat essentially the entire cross-section is stress-
elevated, not just the corners — the decile-trim proxy does not
cleanly separate "corner singularity" from "bulk throat" behavior for
cross. Documented as a known limitation of the method, not adjusted to
force a larger gap (which would be tuning the metric to a preferred
answer).

## Configuration B: comparison and figures

`scripts/compare_config_b.py` reads the converged-level result JSON
per shape, mass from the build-script logs (`build_assembly_*.log`,
grepped — mass is not stored in any results JSON), and recomputes
Config A's mass from its already-validated volume for the baseline
row. Config A's "net-section stress" in the comparison is the
theoretical/FEA-confirmed uniform applied stress (-0.1643 MPa), NOT
computed by the same nodal-mean-with-exclusion method as Config B —
stated explicitly in the script and the output, not treated as
equivalent precision.

**A Python `glob.glob("...{cross,grid}...")` brace-expansion mistake**
was caught before it produced wrong output: Python's glob module does
not support shell-style brace expansion, so the pattern silently
matched nothing and the convergence comparison printed an empty table
with no error. Diagnosed by checking file existence directly, fixed by
using a list of shape names instead of a brace-expansion string.

Outputs: `results/config_b_comparison.csv` (13 rows: 12 designs + 1
Config A baseline), and 4 figures in `figures/` (mass, stiffness,
net-section stress, and stiffness/mass, each vs. removal level, all
four shapes plus the Config A baseline as a dashed reference line).
See README findings section for the results and interpretation.

## Configuration B: status log

- [x] Load-application mechanism designed and verified (roller + fixed base, vs Config A)
- [x] Geometry: 12 designs built, validated, sectional-check bug found and fixed
- [x] Meshing feasibility: all 12 at 3.0mm
- [x] Solve + extract: all 12 at 3.0mm and 2.0mm; cross/grid also at 1.8mm
- [x] Mesh convergence: rect/circle converged at 2.0mm; cross/grid converged at 1.8mm (after a 2-level false read reversed by the 3rd level)
- [x] Comparison table and figures
- [~] Documentation: this update
- [ ] Git commit and push (pending)

## Open items and unverified sources (updated)

- Al 6061-T6 property values and the 276 MPa yield: source to be cited.
- CDS 1U mass (1.33 kg): not verified against the current CDS text.
- Rail-region mean-stress deficit in Config A's assembly: unexplained.
- Why rect and grid track so closely on stiffness/efficiency across all three levels: observed, not investigated.
- Net-section stress corner-exclusion proxy does not work well for cross (gap collapses); documented as a method limitation, not fixed.
- Config B does not have its own roller-base diagnostic (Config A's was treated as validating the mechanism generally, on the prismatic case; not repeated on any pocketed geometry).
- Repository hygiene before commit: Config B's large `.inp` files, mesh logs, and any duplicate/scratch artifacts need the same `.gitignore` review Config A received.

## Configuration B: random-shape exploration (20 -> 223 -> 20 -> 5)

Extension beyond the 4 named shapes: at the highest-scoring removal level
(10%, selected by the combined score below, not by raw stiffness/mass
efficiency alone), explore whether an unconstrained random search over
pocket geometry can beat the best hand-picked named shape.

### Combined selection score (replaces pure stiffness/mass efficiency)

score = (1/3)*(stiffness / Config_A_stiffness)
      + (1/3)*(mass_removed_fraction, i.e. (Config_A_mass - mass)/Config_A_mass)
      + (1/3)*(1 - net_section_stress / 184 MPa allowable)

An earlier version used mass as an inverse ratio (Config_A_mass/mass)
rather than a removed-fraction; that version barely varied across removal
levels (a few percent span) and was rejected because it rewarded removing
LESS mass, the opposite of what a mass-optimization study should reward.
Under the corrected formula, 10% removal scores highest on average across
rect/circle/grid (0.632 vs 0.615 at 20%, 0.608 at 30%), so 10% (not the
originally-picked 20%) is the level used for this exploration.

### Random polygon generator (scripts/build_config_b_random.py)

N points (cycled 5 through 12 across attempts) placed uniformly at random
inside the panel's interior region (73x90mm), connected two ways -- convex
hull, and angular/star order around the centroid -- both guaranteed simple
(non-self-intersecting) by construction. Each candidate polygon is scaled
to the exact target area (830 mm^2, shoelace formula + uniform radial
scaling from centroid) and rejected if any edge falls under the 3mm
minimum feature size, or if it doesn't fit the frame margin, or if it
duplicates an already-accepted shape. Fixed seed (2) for reproducibility.
No cap on the number of valid candidates kept (a "keep the first 20" cap
was tried and removed on request, per the actual research question: how
many good candidates exist, not just find 20). Vertex coordinates for
every accepted candidate are saved to a per-candidate JSON alongside the
isolated-panel .brep, specifically so the full-assembly stage does not
depend on re-deriving them from the seeded random generator a second time.

Result: 223 valid candidate polygons from 150 attempts (0.1-1.6s runtime),
vertex counts spanning 3-12 (hull collapses higher point counts down to
their convex boundary; star-order keeps all N points). All 223 passed
isolated-panel validation (exact area, no slivers, correct bbox, clean
single-solid topology).

### Full-assembly geometry (scripts/build_config_b_random_assembly.py)

Adapts build_config_b_assembly.py's methodology (same rails, same pocket
applied to all 4 panels, same z=2mm section-margin check, same fuse +
cleanliness validation) to read arbitrary polygon vertices from the saved
JSON rather than a closed-form shape formula. All 223 candidates built and
validated in ~30s total; 223/223 passed (exact volume match to the
named-shapes' 10% target 88,660 mm^3, correct bbox, section area 953.0000
at z=2, no stray faces).

### Mesh feasibility (scripts/mesh_config_b_random_feasibility.py)

Linear-mesh-only (no .inp, no solve) feasibility pass on all 223 at all 3
converged mesh levels (3.0, 2.0, 1.8mm): 669/669 succeeded, zero failures,
46 minutes total wall time. Tet counts tightly clustered across all 223
(3.0mm: 27,644-28,835; 2.0mm: 76,356-78,331; 1.8mm: 93,149-101,589), all
well under both ceilings (worst case 44.2% of the element ceiling at
1.8mm). C3D10 node/equation estimates use a promotion ratio (2.049
C3D10-nodes-per-linear-tet) empirically measured on one actual matching
geometry (cross_10pct@3.0mm: 28,313 tets -> 58,000 promoted nodes); these
are estimates for screening purposes, not exact counts.

### Selection funnel: 223 -> 20 -> 5

Given 223 candidates, solving all of them at all 3 mesh levels was judged
disproportionate to what the exploration needed to answer. Funnel used:

1. **223 -> 20** by worst-case (1.8mm) estimated equation-ceiling %,
   ascending -- a resource-similarity filter, NOT a structural-performance
   filter (no candidate had been solved yet at this point). The spread
   was tight (48.5-52.9% across all 223), so this step is closer to a
   tie-breaker than a meaningful safety filter; all 223 were already
   comfortably solvable.
2. **20 solved at 3.0mm** (mesh_config_b.py extended with a --random
   <tag> mode alongside its original <shape> <level_pct> mode). All 20
   solved cleanly, equilibrium <1e-4% throughout, peak RSS tightly
   clustered 602-638 MB. Stiffness spread 472,702-563,179 N/mm (~19%
   range) -- confirming real structural differentiation despite the
   resource-similarity selection.
3. **Combined score computed for all 20 at 3.0mm**; 6 candidates beat the
   best named shape (rect, score 0.6343) -- all 6 were hull-derived,
   none were star-derived.
4. **Top 5 carried to 2.0mm and 1.8mm** (the same 3-level convergence
   treatment every other Config B design received). All 5 solved cleanly
   at both levels (2.0mm: RSS 1.88-1.93 GB; 1.8mm: RSS 2.37-2.41 GB).

### Convergence and final result

| tag | stiff 3.0 | stiff 2.0 | stiff 1.8 | net 3.0 | net 2.0 | net 1.8 |
|---|---|---|---|---|---|---|
| attempt16_hull | 563,179 | 562,961 | 562,908 | 0.2239 | 0.2210 | 0.2210 |
| attempt37_hull | 561,725 | 561,505 | 561,457 | 0.2172 | 0.2160 | 0.2155 |
| attempt54_hull | 560,729 | 560,514 | 560,472 | 0.2417 | 0.2351 | 0.2348 |
| attempt40_hull | 560,221 | 559,923 | 559,867 | 0.2171 | 0.2165 | 0.2167 |
| attempt15_hull | 560,239 | 559,961 | 559,901 | 0.2260 | 0.2255 | 0.2240 |

Stiffness fully converged for all 5 (2.0->1.8mm shift <=0.01%). Net-section
stress converged cleanly (2.0->1.8mm shift <=0.63% for all 5) -- unlike
grid (named shapes), none of these 5 needed the third level to REVERSE an
earlier misleading trend; they simply settled quickly and stayed settled.

**Final ranking at the converged (1.8mm) level, all 5 beat rect (score
0.6343):**

| rank | tag | stiffness | net stress | score | margin over rect |
|---|---|---|---|---|---|
| 1 | attempt16_hull | 562,908 | 0.2210 | 0.6372 | +0.46% |
| 2 | attempt37_hull | 561,457 | 0.2155 | 0.6365 | +0.35% |
| 3 | attempt54_hull | 560,472 | 0.2348 | 0.6360 | +0.27% |
| 4 | attempt15_hull | 559,901 | 0.2240 | 0.6357 | +0.22% |
| 5 | attempt40_hull | 559,867 | 0.2167 | 0.6357 | +0.22% |

All 5 winning polygons are 3-4 vertex (triangle/quad) shapes -- convex-hull
outputs from 5-point random seeds. No star-order (more irregular, more
concave) polygon appeared among the top 6 at any stage. Stress contours
(figures/contours/config_b_random_attempt{16,37,54,15,40}_hull_h1p8_vm_
contour.png) confirm the pattern visually: stress concentrates at the
polygon's more acute vertices (typically 2-3 hot spots per panel), and
peak von Mises does NOT track the combined-score ranking -- attempt15_hull
has the highest peak among the 5 winners (0.7775 MPa) despite ranking 4th
on the combined score, because the score is driven by net-section
(spatially-averaged) stress and stiffness, not the single worst local
point. This is stated explicitly so the contours are not mis-read as
implying peak-stress ordering matches the ranking.

### Honest scope of this result

This identifies the best design among 20 resource-screened candidates out
of 223 generated -- not a claim of finding the global optimum over the
full space of possible pocket shapes, nor even over all 223 generated
candidates (only 20 were ever solved). The margin over rect is real and
mesh-converged but modest (0.2-0.5%). The result does support a real,
physically sensible qualitative finding: smoother/more-convex pocket
boundaries (fewer, gentler-angled vertices) outperform highly irregular
ones at this removal level and load case, consistent with reduced local
stress concentration.

### Bugs found and fixed during this stage

1. `occ.extrude()` returns a flat list of (dim,tag) tuples, not a
   (list,tag) 2-tuple like fuse/cut/intersect -- copy-pasted the wrong
   unpacking pattern from the named-shape builder; caused an immediate
   crash on the first random-assembly build. Fixed by filtering the flat
   list for dim==3 entries directly.
2. `mesh_config_b.py`'s C3D10 promotion loop originally excluded the
   *EQUATION reference node from statistics correctly, but an earlier
   attempt to add the same exclusion to plot_stress_contour.py's printed
   min/max used the wrong node list (the 4-corner-only `conn` array used
   for surface rendering, not the full 10-node-per-element connectivity),
   which marked ~2/3 of all nodes (every midside node) as "unconnected."
   This silently corrupted the printed max von Mises for every
   *EQUATION-loaded case (e.g. cross_10pct_h1p8 read 1.2259 MPa instead
   of the correct 1.4462 MPa) while leaving the actual contour PNGs
   unaffected (they use the full node set). Fixed by parsing full
   10-node element connectivity separately for the connectivity mask.
   All 37 affected contours were re-generated and cross-checked against
   two independent prior extraction paths after the fix.
3. `main()`'s .inp header-writing line in the --random branch of
   mesh_config_b.py referenced `shape`/`level`, which only exist in the
   named-shape branch, causing an UnboundLocalError after mesh generation
   had already completed successfully. Fixed by using a `label_for_print`
   string built consistently in both branches.
4. A manual multi-step edit to development-log.md (deleting the file
   and pasting only new content, rather than appending) destroyed the
   full ~333-line Config A history; recovered via `git show HEAD:...`
   against the last commit (02eac59) and a verified merge, since the
   Config A content had already been committed. No data was permanently
   lost, but this is recorded as a reminder that `git`-committed content
   is the safety net for exactly this kind of accidental overwrite.

