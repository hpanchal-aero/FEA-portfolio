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
