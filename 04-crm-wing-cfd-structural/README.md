# Project 04: CFD-to-structure analysis of a transonic airliner wing (uCRM-9)

One-way and once-corrected aerodynamic loading of the NASA Common Research Model wing: steady transonic CFD in OpenFOAM, load mapping to a CalculiX shell model of the uCRM-9 wingbox, and a deflected-geometry correction step to test how much static aeroelastic deformation changes the aerodynamic load and the structural response.

Status: documentation in progress. Sections below are added as they are verified against the result files in `results/`.

## 1. Overview

A transonic cruise condition (Mach 0.85, lift coefficient 0.5 on the reference area of the full wing) is solved with a steady compressible OpenFOAM model with slip walls (Euler-type, see Section 9) on a snappyHexMesh mesh of the wing alone. The surface pressure is converted into force- and moment-conserving nodal loads on the uCRM-9 wingbox shell model, which is solved statically in CalculiX. The structural displacement is then fed back once: the CFD wing is deformed (translation and twist per spanwise section), the CFD is rerun and re-trimmed to the same lift, and the new loads are mapped back to the undeformed structural frame.

Headline result, with the caveats in the Limitations section (unconverged meshes, first-order upwind Euler loads, one correction step only): at equal lift (CL = 0.5) the corrected solution has a tip deflection of 1.907 m against 2.801 m for the one-way solution on the same mesh (-31.9%), and an outer-wing twist of -4.72 degrees against -6.17 degrees (-23.5%). The one-way load is therefore not an adequate load source for this wing at this level of deformation.

## 2. Research question

How much does static aeroelastic deformation change the aerodynamic loads and the structural response (tip deflection and twist) of a transonic airliner wing, and is a one-way CFD-to-structure load transfer adequate, or does it need at least one correction step with the deformed shape?

The adequacy criterion was defined relative to the structural model's own bias against the reference (decision D3, refined in D22): one correction step on one mesh level, compared against one-way and against the dataset's reference loads.

## 3. Motivation and engineering significance

High-aspect-ratio transport wings bend and twist measurably under cruise load. The twist changes the local angle of attack, which changes the load, which changes the deformation. A one-way transfer (rigid CFD loads onto the structure) ignores this loop and is the usual first step in preliminary structural sizing. Whether it is conservative or not depends on the case, and the size of the effect decides whether a coupled analysis is needed. This project quantifies that effect for a public benchmark wing, using open tools only and a workflow that can be rerun from the terminal.

## 4. Scope and what this is not

- Steady, static, one correction step. No full two-way iteration to convergence, no flutter, no dynamics.
- Wing alone, extended to the symmetry plane (decision D11). No body, tail, engines or fuselage interference. The dataset's own reference loads are for a different configuration (wing-body-tail), which is documented where they are used.
- Euler-type CFD: slip walls, no turbulence model and no boundary layer. The laminar viscosity is kept at its physical value only as numerical regularisation (a zero viscosity produced NaN in the first momentum solve). Shock position and load distribution are therefore a lower-fidelity load source than a viscous solution would give (decision D17).
- Stress is not evaluated: the structural results kept are displacements and section kinematics.
- This is a verification-oriented study, not a validation against wind-tunnel or flight data.

## 5. Workflow
uCRM-9 geometry --> snappyHexMesh (L8 / L9 / L10)
|
OpenFOAM foamRun, steady Euler, first-order upwind
|
surface pressure --> strip-wise force/moment-conserving nodal loads
|
CalculiX S4 shell model of the uCRM-9 wingbox --> displacement
|
section translation + twist --> deformed CFD wing --> mesh --> CFD (re-trimmed to CL 0.5)
|
loads mapped back to the undeformed frame --> CalculiX --> corrected response

All steps are scripted (`scripts/`), run from the terminal under WSL2/Ubuntu 24.04, and post-processed with Python (NumPy, SciPy, Pandas, Matplotlib, PyVista).

## 6. Data sources, licence and attribution

- Wing geometry and structural model: uCRM, the undeflected Common Research Model, by T. R. Brooks, G. K. W. Kenway and J. R. R. A. Martins. Dataset: Mendeley Data, DOI 10.17632/gpk4zn73xn.1, licence CC BY 4.0. Paper to cite: Brooks, Kenway and Martins, "Benchmark Aerostructural Models for the Study of Transonic Aircraft Wings", AIAA Journal 56(7), 2018, pp. 2840-2855, DOI 10.2514/1.J056603.
- Underlying aerodynamic configuration: NASA Common Research Model. Vassberg, DeHaan, Rivers and Wahls, "Development of a Common Research Model for Applied CFD Validation Studies", AIAA Paper 2008-6919, 2008.
- The dataset files are not stored in this repository. They are downloaded by script into `data/` (ignored by Git). The meshes, decks, nodal loads and section tables built from them are regenerated by the scripts and are also not committed.
- Numbers quoted from the dataset (the NASTRAN reference tip deflection, shell mass, reference aerodynamic load) and the figures that show the wing outline are derived from the CC BY 4.0 material and are attributed here. This project is not affiliated with, or endorsed by, the dataset authors or NASA.

## 7. Computational domain and operating condition

Geometry: the uCRM-9 wing, taken from the dataset's wing-body-tail geometry and rebuilt as a closed watertight surface (STL) with a closed tip and a flat closure on the symmetry plane (`scripts/build_wing_geometry.py`, `scripts/build_wing_stl.py`, `scripts/wing_edge_matching.py`). The wing extends over x = 23 to 48 m, y = -29.5 to 29.5 m, z = 2.4 to 4.9 m (metres throughout). The CFD solves the half-space y >= 0.

Far field: a box of 500 m x 250 m x 500 m (x from -215 to 285 m, y from 0 to 250 m, z from -246 to 254 m), i.e. about 240 m (about 34 reference chords) upstream, downstream, above and below the wing and 220 m (about 7.5 half-spans) beyond the tip. The plane y = 0 is a symmetry plane.

| Quantity | Value |
|---|---|
| Mach number | 0.85 |
| Static temperature | 216.65 K |
| Static pressure | 21,662.7 Pa |
| Density | 0.34833 kg/m3 |
| Speed of sound / speed | 295.07 m/s / 250.81 m/s |
| Dynamic pressure | 10,955.9 Pa |
| Reference area (half wing) | 191.87 m2 (full wing 383.74 m2) |
| Reference length (mean aerodynamic chord) | 7.01 m |
| Moment reference point | (33.678, 0, 4.520) m |
| Dynamic pressure x reference area | 2.102 MN |

Gas: perfect gas, gamma = 1.4, R = 287.058 J/(kg K), Cp = 1004.703 J/(kg K). The angle of attack is set per run in `0/U`.

## 8. CFD mesh

Mesher: snappyHexMesh (OpenFOAM 11) on a 50 x 25 x 50 background block of 10 m cells. The wing surface is refined to level (7, L) with feature-edge refinement at level L, distance-based refinement regions around the wing (levels 7, 6, 5 and 4 at 0.15, 0.4, 1.0 and 3.0 m), three cells between refinement levels and a feature angle of 30 degrees. The mesh is snapped to the surface. No prism layers are used (`addLayers false`): two layer variants reached only 40.1% and 37.1% wing-face coverage against the 90% gate, so the project moved to a slip-wall (Euler-type) model on a mesh without a boundary layer (decision D17).

The nominal cell edge at the finest surface level L is 10 m / 2^L: 39 mm (L8), 19.5 mm (L9), 9.8 mm (L10). The level sweeps are recorded in `results/cfd_snappy_level_sweep_novolq.txt` and `results/cfd_snappy_level_sweep_novolq_defl.txt` (script `scripts/snappy_level_sweep.sh`).

Production meshes and their `checkMesh` results (re-run on the run directories, recorded in `results/production_meshes.txt` by `scripts/record_production_meshes.sh`):

| Mesh | Used for | Cells | Max non-orthogonality | Max skewness | Max aspect ratio | Min cell volume (m3) | checkMesh |
|---|---|---|---|---|---|---|---|
| L8 | mesh study | 636,768 | 64.3 deg | 3.98 | 6.29 | 7.4e-6 | OK |
| L9 | mesh study, base mesh, alpha search | 1,103,898 | 65.1 deg | 3.98 | 5.45 | 9.8e-7 | OK |
| L10 | mesh study | 2,018,605 | 64.6 deg | 3.98 | 5.97 | 1.0e-7 | 1 check failed: 2 incorrectly oriented faces in the face-pyramid check |
| L9, deflected geometry | corrected-load runs | 1,116,145 | 65.2 deg | 3.49 | 5.28 | 9.6e-7 | OK |

The production meshes were generated with the snappyHexMesh minimum-volume check switched off (`minVol -1e30`) and with tighter skewness limits (internal 3.5, boundary 4) than the standard variants (`minVol 1e-13`, internal 4, boundary 20); the main quality-control non-orthogonality limit is 65 in both. The standard variants give a different snapped mesh for the same cell count (L9: 1,421,908 points against 1,392,709; L10: 2,641,163 against 2,561,980) and were not used for the runs. All minimum cell volumes are positive. The measured maximum non-orthogonality of the L9 and deflected meshes (65.1 and 65.2 degrees) is slightly above the 65 degree control; this is not explained here and the meshes pass the default `checkMesh` limit of 70 degrees. The two incorrectly oriented faces on L10 are at the blunt trailing-edge base and were accepted knowingly (the run was monitored and the defect is documented).

The meshes of the mesh study are not nested (volume ratios 1.20 and 1.22 between successive levels, surface-cell ratios 1.34 and 1.37 instead of a constant ratio).

Trailing-edge resolution. The blunt trailing edge was kept (decisions D13, D16). The base thickness of the wing surface, measured on the STL with `scripts/te_thickness.py`, is 12.6 to 17.5 mm for y = 0.5 to 29.2 m (17.5 mm at the root station, 14.5 mm at y = 3 m, 12.6 to 13.9 mm elsewhere). Against the nominal cell edge this is 0.32 to 0.45 cells across the base on L8, 0.65 to 0.90 on L9 and 1.3 to 1.8 on L10. These ratios are nominal (cell edge from the refinement level, not measured on the mesh). The base is therefore under-resolved on L8 and L9 and only marginally resolved on L10, which affects the pressure near the trailing edge and the mesh-study differences.

## 9. CFD model

Solver: OpenFOAM 11, `foamRun -solver fluid`, steady-state compressible flow with a PIMPLE pressure-based algorithm. Thermophysics: `hePsiThermo`, perfect gas, constant-Cp enthalpy, laminar transport with the physical dynamic viscosity (1.4216e-5 Pa s, Pr 0.71) and no turbulence model. The viscous terms are kept only because a zero viscosity gave NaN in the first momentum solve on this mesh; with slip walls there is no boundary layer and the flow is Euler-type at the wall.

| Boundary | Velocity | Pressure | Temperature |
|---|---|---|---|
| inlet, outlet, farfieldY, bottom, top | freestreamVelocity | freestreamPressure | inletOutlet |
| wing | slip | zeroGradient | zeroGradient |
| symmetry (y = 0) | symmetryPlane | symmetryPlane | symmetryPlane |

Schemes (`fvSchemes.upwind`): first-order bounded upwind for momentum and energy, linear Gauss gradients with a cell limiter for the velocity gradient, linear corrected Laplacian. First-order upwind is the only scheme that converged robustly on this case (linearUpwind hit the pressure floor from about iteration 7; TVD developed pressure-floor columns at the root and tip, decisions D18 and D19). It adds numerical diffusion and smears the shock, which is why the CFD is treated as a lower-fidelity load source.

Linear solvers and relaxation: pressure GAMG (DIC smoother, relTol 0.01), U and h PBiCGStab with DILU (relTol 0.1), one non-orthogonal corrector, relaxation p 0.3, rho 0.01, U and h 0.7. Residual targets p 1e-6, U 1e-5, h 1e-5.

Limiters (`fvConstraints`): `limitPressure` with a floor of 0.1 p_inf = 2,166.27 Pa and a ceiling of 2 p_inf = 43,325.4 Pa, and `limitTemperature` between 100 K and 400 K. Activity recorded in the run logs: in the L8, L9, L10 and deflected alpha 0.36 degree runs `limitPressure` printed messages only during the first few iterations (the start-up transient); in the deflected alpha 2.07 degree run it is silent from iteration 7 to 1419 and then acts at every iteration from 1420 to 1500, while the pre-limit minimum pressure falls from 2,123 Pa to -1,075 Pa. The loads of that run are taken from iteration 1300, before any limiter activity, but the minimum pressure and temperature of the saved fields already fall steadily from iteration 1300 onward (Section 12.2). No `limitTemperature` activity message appears in any of these logs (whether this limiter reports its own activity was not verified).

Force coefficients (`forceCoeffsCompressible`): forces on the `wing` patch with the reference pressure 21,662.7 Pa subtracted, so that the open surface carries no spurious uniform load. Forces use the actual density field; the dynamic pressure uses the free-stream density. Reference area 191.87 m2 (half wing), reference length 7.01 m. The absolute drag coefficient (about 0.05) is not used in any conclusion.

Run procedure: `scripts/cfd_run_upwind.sh <alpha_deg> <tag> [iterations] [processes]` copies the template case `simulation/cfd_v1/euler`, sets the angle of attack, switches to the upwind schemes, decomposes, and runs `foamRun` in parallel (default 1500 iterations, 8 processes). The case templates and the per-run files that differ from the template are in `simulation/cfd_v1/`.

## 10. Structural model (CalculiX)

Source model: the uCRM-9 coarse wingbox in Nastran format, with 25,055 CQUADR shell elements, 23,738 grid points, 242 shell property sets and 968 MAT2 material cards (four per property set: membrane, bending, transverse shear and membrane-bending coupling). The stiffeners are smeared into the shell properties and there are no composite laminates. The dataset provides a Nastran static solution (gravity of 9.81 m/s2 in z plus aerodynamic forces) and a modal solution (20 modes) as the reference.

CalculiX deck (`scripts/ucrm_to_ccx.py`, `scripts/build_structure_decks.sh`): S4 shell elements with one shell section and one orthotropic material (`*ELASTIC, TYPE=ENGINEERING CONSTANTS`) per property set, root constraints from the dataset's SPC set, and 74 T3D2 truss pairs (E = 1e9 Pa, unit area) that carry the concentrated masses. The deck has 23,886 nodes and 25,129 elements. Property set 1, for example, has E1 = 70.0 GPa, E2 = 31.97 GPa and nu12 = 0.3. The material data come from the uCRM-9 dataset.

Constitutive idealisation (documented assumption): the membrane matrix of each property set is inverted to E1, E2, nu12 and G12; the transverse shear moduli are taken from the shear material card. The bending matrix, the membrane-bending coupling matrix and the 12I/t3 factor are not represented. E3, nu13 and nu23 are placeholders required by CalculiX's three-dimensional constants.

Masses: the shell mass in the deck is 11,958.49 kg (dataset: 11,958.5 kg, ratio 1.0000); concentrated masses of the dataset (engine, leading- and trailing-edge and fuel masses), 22,170.0 kg, are tied to the structure with `*EQUATION` relations. The total mass is 34,128.5 kg and the gravity load is -334.80 kN.

Solver: CalculiX 2.23, SPOOLES direct solver, one linear static step (`*STATIC`, no `*NLGEOM`). The coarse decks solve in about 10 s (peak memory about 1.5 GiB for the corrected CL 0.5 deck). The medium dataset mesh cannot be solved in this build (the SPOOLES solver crashes and PARDISO is not linked, decision D8), so all structural results are on the coarse mesh and the structural discretisation error is not quantified. Stress is not evaluated; the kept results are displacements and section kinematics.

Verification against the dataset reference (structure only, with all dataset loads: aerodynamic, gravity, engine, leading/trailing-edge and fuel masses):

| Quantity | CalculiX | NASTRAN reference | Difference | Tolerance set beforehand | Result |
|---|---|---|---|---|---|
| Tip deflection, mean of the 21 tip-section nodes | 2.3368 m | 2.5239 m | -7.41% | 2% | not met (documented limitation, decision D10) |
| Tip deflection, maximum | 2.3643 m | 2.5239 m | -6.32% | | |
| Modes 1 to 5 | 1.675, 2.582, 3.867, 4.760, 6.338 Hz | 1.629, 2.632, 3.816, 4.684, 6.315 Hz | +2.83, -1.88, +1.34, +1.64, +0.36% | 5% | met (worst 2.83%) |

CalculiX modes 6, 8 and 9 (7.86, 9.84 and 10.01 Hz) are local rib-web panel modes of the two root-rib property sets and have no counterpart in the NASTRAN list; multiplying the elastic constants of those two property sets by 10 removed them and moved the other modes by less than 1%. The global CalculiX modes at 8.25 and 10.67 Hz are within +1.7% and +0.8% of the NASTRAN modes at 8.11 and 10.59 Hz. The modal criterion uses only modes 1 to 5.

Equilibrium and mass checks (`scripts/check_equilibrium.py`): the shell mass matches the dataset, the sum of the aerodynamic nodal forces is +967.01 kN, gravity is -334.80 kN, and the vertical force balance closes to a few newtons (-2.5 N in the coupled cases).

The static tip deflection is therefore 7.4% below the reference and this bias is carried into every comparison with the dataset. An analysis of the neglected terms (`scripts/section_stiffness.py`) gives a NASTRAN-to-deck axial stiffness ratio of 0.92 (median over the property sets, range 0.87 to 0.98) when the membrane-bending coupling is left free and 1.00 when it is restrained; the observed tip ratio of 0.926 to 0.937 agrees in sign with the free case. This is consistent with the neglected coupling being the cause but does not prove it (the analysis assumes plane sections, no shear lag and uniaxial stress). A layered-shell representation of the full stiffness matrices was investigated and stopped: CalculiX 2.23 allows composite sections only for S8R and S6 elements, and the composite model did not verify against laminate theory (decisions D9 and D10).

## 11. Load transfer and correction procedure

1. CFD surface force. The wing force is F = sum of (p - p_inf) S_f over the wing faces, with S_f the face area vector pointing out of the fluid domain: pressure only, relative to the free-stream pressure (`scripts/cfd_surface_loads.py`).
2. Load mapping (`scripts/map_cfd_loads.py`, decisions D20 and D21). Target nodes are the nodes the dataset loads (y >= 2.94 m) plus the free skin nodes of the root stub (y < 2.94 m, shell elements with |n_z| >= 0.7, not constrained); the dataset applies no load inboard of 2.94 m, but the CFD does. The baseline puts each CFD face force on its nearest target node. A correction per spanwise strip (1 m by default) then restores the exact strip moment about the strip centroid with a minimum-norm, force-only change, so the total force is unchanged. Loads are force-only nodal loads, as in the dataset.
3. Section kinematics (`scripts/section_kinematics.py`). For each spanwise bin (1 m by default) the vertical displacement is fitted as uz = a + b (x - xbar) + c (y - ybar) by least squares, together with the mean ux and uy. The twist is -b (negative means nose-down washout). A root point with zero motion is prepended and a reconstruction of the displacement from the section values is checked.
4. Deformed CFD wing (`scripts/deform_wing_stl.py`). Each vertex is displaced by the section translation and twist, interpolated linearly in y, with zero motion at the root and constant extrapolation beyond the outermost section (y = 29.24 m). The displacement field is first order: translation and twist per section, no rotation about the x-axis. Normals are recomputed from the deformed triangles.
5. CFD on the deformed wing. A new level-9 mesh is generated for the deformed surface (1,116,145 cells) and the same flow condition is solved, first at the original angle of attack (0.36 degrees) and then, as a re-trim, at 2.07 degrees (decisions D22 and D23).
6. Pullback (`--pullback` option of `map_cfd_loads.py`). The face centres of the deformed CFD wing are mapped back to the undeformed (jig) frame of the structural nodes by solving c0 + D(c0) = c. The force vectors are kept as computed, so the loads are applied in the jig frame (first order, linear structure).
7. Lift trim by interpolation (`scripts/blend_wing_vtk.py`, decision D24). The wing-patch pressure of the two deformed states is blended linearly, p = (1 - w) p1 + w p2 with w = (CL_target - Cl1) / (Cl2 - Cl1), where Cl is computed from the surface pressure with q = 10,955.9 Pa and Aref = 191.870 m2. This is exact for a linear structure and assumes that the aerodynamic load is linear in the angle of attack between 0.36 and 2.07 degrees.
8. One correction step only. The loop is not iterated to convergence.

## 12. CFD results

### 12.1 Base angle of attack

The lift of the undeflected wing was searched with first-order upwind runs on the L9 mesh (`scripts/alpha_cl_table.py`, `results/alpha_cl_upwind.txt`). At iteration 1500:

| Angle of attack | CL (iteration 1500) |
|---|---|
| 0.30 deg | 0.49347 |
| 0.36 deg | 0.49976 |
| 0.50 deg | 0.51440 |

The lift-curve slope between 0.30 and 0.50 degrees is 0.1046 per degree and the linear estimate of the angle of attack for CL = 0.5 is 0.362 degrees, so 0.36 degrees is the base angle of attack of all undeflected runs. The undeflected wing reaches CL = 0.4998 on the L9 mesh. The deflected wing at the same angle of attack carries CL = 0.341 (a loss of 31.8% of the lift), and the re-trimmed deflected run at 2.07 degrees reaches CL = 0.525 (decisions D12, D19, D23).

### 12.2 State of convergence

The runs are steady-state with fixed iteration counts (1500 for the L9 and deflected runs). The L8 and L10 runs stopped on the residual control at iterations 1161 and 1148. The forces are still drifting slowly at the end of the runs: for the undeflected wing at 0.36 degrees on L9, Cl is 0.50129 at iteration 1148, 0.50054 at 1300 and 0.49976 at 1500 (about -0.1% per 100 iterations); for the deflected wing at 2.07 degrees Cl is 0.52572, 0.52501 and 0.52425 at 1148, 1300 and 1500. The solutions are therefore not converged to a strict steady state, and the iteration used as the load source is stated for every result.

Limiter activity (`scripts/check_limiter_fields.py`, `results/limiter_fields.txt`): at the iterations used as load sources (L8 1161, L9 1500, L10 1148, deflected alpha 0.36 degrees 1500, deflected alpha 2.07 degrees 1300) no cell of the saved field is at a pressure or temperature limiter bound. The deflected alpha 2.07 degree run develops a growing suction region: its minimum pressure is 6,553 Pa at iteration 1300, 3,252 Pa at 1400 and 2,166 Pa (the floor, with 99 cells of the saved VTK field at the floor) at 1500, and its minimum temperature is 161 K, 142 K and 116 K. The other load-source states have minimum pressures of 5,518 to 13,310 Pa and minimum temperatures of 182 to 196 K. The loads of the alpha 2.07 degree run are taken at iteration 1300, which is free of clipping but already belongs to a field in which this region is growing.

### 12.3 Mesh sensitivity (L8, L9, L10, alpha 0.36 degrees, undeflected wing)

`scripts/mesh_study.py`, `results/mesh_study_p0p36.txt`, `figures/mesh_study_p0p36.png`. Force coefficients at the common iteration 1148 (common window 948 to 1148); the tip deflection is that of the one-way CalculiX run driven by the loads of the last iteration of each CFD run (1161, 1500 and 1148), so this column does not use the same iteration as the coefficients.

| Mesh | Cells | Wing faces | Cl | Cd | Cm | Tip uz (m) | Max uz (m) |
|---|---|---|---|---|---|---|---|
| L8 | 636,768 | 109,730 | 0.48012 | 0.04718 | -0.16789 | 2.5960 | 2.6294 |
| L9 | 1,103,898 | 198,337 | 0.50129 | 0.04819 | -0.18714 | 2.8010 | 2.8378 |
| L10 | 2,018,605 | 372,363 | 0.51551 | 0.04911 | -0.19866 | 2.9303 | 2.9691 |

Successive changes (L8 to L9, then L9 to L10): Cl +4.41% and +2.84%, Cd +2.15% and +1.91%, Cm magnitude +11.5% and +6.2%, tip deflection +7.9% and +4.6%. The values are still changing by several percent per refinement, so the mesh is not converged.

Grid convergence index (Celik et al.), refinement ratio from the cell counts (r = 1.2013 for L8 to L9 and 1.2229 for L9 to L10), fine pair L10 and L9:

| Quantity | Observed order | Richardson extrapolation | GCI, observed order | GCI, order 1 |
|---|---|---|---|---|
| Cl | 2.56 | 0.5366 | 5.11% | 15.47% |
| Cd | 0.92 | 0.0537 | 11.63% | 10.61% |
| Cm | 3.17 | -0.2116 | 8.14% | 32.52% |
| Tip deflection | 2.89 | 3.094 m | 6.99% | 24.75% |

The observed orders of Cl, Cm and the tip deflection (2.6 to 3.2) are well above the formal first order of the upwind scheme, which means that the three meshes are not in the asymptotic range (the meshes are not nested, the refinement ratios are small and the blunt trailing edge is under-resolved). The GCI with the observed order is therefore not reliable and the conservative value with order 1 is the better indication: 10 to 33% for the coefficients and about 25% for the tip deflection. With the refinement ratio defined from the wing-surface face counts instead (1.344 and 1.370) the observed orders are 1.5 (Cl), 0.5 (Cd), 1.9 (Cm) and 1.7 (tip), and the order-1 GCI is 6.4 to 19.6%, so the numerical uncertainty estimate itself depends on the definition of the refinement ratio. The Richardson-extrapolated values are given for information and are not used. All absolute values in this study carry a mesh uncertainty of the order of 5 to 25%; the comparisons between one-way and corrected loads (Section 13) are made on the same mesh.

## 13. Coupling results

`scripts/coupling_table.py`, `results/coupling_summary.csv`, `results/coupling_summary.md`. All cases use the same coarse CalculiX deck, gravity and point masses (gravity -334.8 kN); the vertical force closes to 2.5 N in every case.

| Case | Load source | CFD CL | Aerodynamic Fz (kN) | Tip uz (m) | Tip uz vs dataset load | Twist, outer section (deg) |
|---|---|---|---|---|---|---|
| Dataset loads | dataset reference loads (wing-body-tail) | - | 967.0 | 2.3368 | - | -4.618 |
| One-way, L9 | undeflected wing, alpha 0.36 deg, iteration 1500 | 0.4998 | 1051.2 | 2.8010 | +19.9% | -6.167 |
| One-way, L10 | undeflected wing, L10 mesh, iteration 1148 | 0.5155 | 1084.3 | 2.9303 | +25.4% | -6.489 |
| Corrected, same alpha | deflected wing, alpha 0.36 deg, iteration 1500 | 0.3409 | 717.2 | 1.1071 | -52.6% | -3.562 |
| Corrected, alpha 2.07 | deflected wing, alpha 2.07 deg, iteration 1300 | 0.5250 | 1106.8 | 2.0327 | -13.0% | -4.900 |
| Corrected, CL 0.5 | deflected wing, interpolated in lift between the two deflected runs | 0.5000 | 1053.8 | 1.9068 | -18.4% | -4.718 |

The twist of the outermost 1 m section (centred at y = 29.24 m) is the negative of the slope b of the plane fit uz = a + b (x - xbar) + c (y - ybar); negative means nose-down washout. It is not validated against the NASTRAN reference. The dataset row uses the dataset's own wing-body-tail reference loads; the other rows use wing-alone CFD loads on the same deck, so the comparison with the dataset row mixes the difference in loads with the structural bias of -7.41% (Section 10).

Findings, all on the L9 mesh and the coarse structural deck:

- The one-way wing-alone load at CL = 0.4998 is 8.7% larger than the dataset load (1051.2 against 967.0 kN) and gives a tip deflection 19.9% larger and an outer-section washout of -6.17 degrees against -4.62 degrees.
- Deforming the wing at the same angle of attack removes 31.8% of the lift (CL 0.4998 to 0.3409): the structure bends and twists nose-down, which reduces the local angle of attack of the outer wing. The corrected tip deflection is 60.5% smaller than the one-way value in that case, but this compares different total loads.
- At equal lift the comparison is meaningful. The corrected solution at CL = 0.5 carries 1053.8 kN against 1051.2 kN for the one-way solution (+0.25%), yet its tip deflection is 1.907 m against 2.801 m (-31.9%) and its outer-section twist is -4.718 degrees against -6.167 degrees (-23.5%). The total vertical force is the same to within 0.3%, so the difference comes from the spanwise and chordwise distribution of the load: the washout moves load inboard and reduces the outer-wing load that drives the tip deflection.
- The L8 to L10 mesh effect on the one-way tip deflection (+4.6% from L9 to L10) is much smaller than the effect of the correction (-31.9%), but the GCI of the tip deflection (Section 12.3) is of the order of 25%, so the absolute values should not be over-interpreted.
- With respect to the criterion of decision D3 and D22 (is one-way adequate?): in this model the one-way load overestimates the tip deflection and washout by about 47% and 31% relative to the corrected solution (2.801 / 1.907 and 6.167 / 4.718), so one correction step changes the structural response substantially and the one-way load is not an adequate load source for this wing at this level of deformation.

## 14. Verification and validation

Verification (was the mathematical model solved correctly?):

| Item | Evidence | Outcome |
|---|---|---|
| Free-stream condition | Mach number of the free-stream field and `flight_conditions` | Mach 0.85, a = 295.07 m/s, U = 250.81 m/s |
| CFD mesh quality | `checkMesh` on the four production meshes (`results/production_meshes.txt`) | positive cell volumes, max non-orthogonality 64 to 65 deg, three meshes pass, L10 has 2 incorrectly oriented faces |
| CFD mesh convergence | three levels, GCI (Section 12.3) | not converged, asymptotic range not reached, GCI (order 1) 10 to 33% |
| CFD iterative convergence | force histories (`figures/force_history_mesh_study.png`, `figures/force_history_coupling.png`) | forces drift by about 0.1% per 100 iterations at the end |
| Limiters | saved fields (`results/limiter_fields.txt`) | no cell at a limiter bound at the load-source iterations |
| Structure, mass and equilibrium | `check_equilibrium.py` | shell mass ratio 1.0000, force balance closes to 2.5 N |
| Structure against the Nastran reference | `post_static.py`, `post_modal.py` (Section 10) | static tip -7.41% (tolerance 2%, not met), modes 1 to 5 within 2.83% (tolerance 5%, met) |
| Load mapping | `Fz / (q Aref)` against the CFD CL | within 0.3% for every case (Section 13) |
| Displacement mapping for plots | `check_mapping_vs_dat.py` against the printed CalculiX displacements | maximum difference 4.1 mm and RMS 32 um (L9 case); 3.3 mm and 25 um (corrected CL 0.5 case) |

Validation (does the model represent reality?): none was performed. The CFD was not compared with wind-tunnel data or with published CRM solutions, and the structural model was compared only with another computational model (Nastran) for the same structure, which is code-to-code verification and not validation. The conclusions in this report are therefore conclusions about the computational model and not about the real aircraft.

## 15. Limitations

Aerodynamic model
- Euler-type flow with slip walls and first-order upwind schemes: no boundary layer, no shock-boundary-layer interaction, no separation, and a smeared shock. The loads are a lower-fidelity source. The absolute drag coefficient (about 0.05) is not used.
- The CFD mesh is not converged (tip deflection +7.9% and then +4.6% per refinement, Cl +4.41% and then +2.84%). The observed orders (2.6 to 3.2 for Cl, Cm and tip deflection) exceed the formal first order, so the asymptotic range is not reached; the three meshes are not nested (volume ratios 1.20 and 1.22, surface-cell ratios 1.34 and 1.37). Only three levels exist.
- The blunt trailing edge (12.6 to 17.5 mm) is under-resolved on L8 and L9 (0.3 to 0.9 cells across the base, nominal) and only marginally resolved on L10 (1.3 to 1.8). The L10 mesh has two incorrectly oriented faces at the trailing-edge base, accepted knowingly.
- The production meshes were generated with the snappyHexMesh minimum-volume check disabled (tighter skewness limits were used instead); all minimum cell volumes are positive. The measured maximum non-orthogonality of L9 and the deflected mesh (65.1 and 65.2 degrees) is slightly above the 65 degree control and is not explained.
- The runs stop at fixed iteration counts (or on the residual control at 1161 and 1148 iterations for L8 and L10) with the forces still drifting by about 0.1% per 100 iterations. The load-source iteration differs between runs (L8 1161, L9 1500, L10 1148, deflected alpha 0.36 degrees 1500, deflected alpha 2.07 degrees 1300); the L9 value of Cl differs by 0.3% between iterations 1148 and 1500. In the mesh-study table the tip deflections and the coefficients use different iterations (Section 12.3).
- The deflected alpha 2.07 degree run develops a growing suction region (Section 12.2). Its loads are taken at iteration 1300, before any limiter activity. The tip deflection computed with the iteration 1500 loads is 2.0304 m against 2.0327 m with the iteration 1300 loads (-0.11%), so the structural result is insensitive to that choice; the iteration 1500 field itself is clipped at the pressure floor.
- The deflected-wing loads at CL = 0.5 are interpolated linearly in lift between two CFD states (alpha 0.36 and 2.07 degrees, weight of about 0.86 on the second state), assuming that the aerodynamic load is linear in the angle of attack in this range.
- About 0.7% of the wall faces show scatter in the pressure coefficient. Its effect on the integrated force is estimated to be below 0.1% (rough bound).
- The CFD is for the wing alone extended to the symmetry plane, whereas the dataset loads are for the wing-body-tail configuration. At CL = 0.5 the wing-alone load is 8.7% larger than the dataset load.

Structural model
- Coarse structural mesh only: the medium mesh cannot be solved in this CalculiX build (decision D8), so there is no structural mesh study and the structural discretisation error is not quantified.
- The shells represent only the membrane and transverse-shear stiffness of the smeared stiffened panels; bending and membrane-bending coupling are neglected. The static tip deflection is 7.41% below the Nastran reference (tolerance 2%, not met), which is consistent in sign with the neglected coupling but not proven to be caused by it. A layered-shell representation was investigated and stopped (decisions D9 and D10).
- The analysis is linear static (no `NLGEOM`). The one-way tip deflection is 2.80 m, 9.5% of the semispan of 29.38 m (corrected: 1.91 m, 6.5%); the effect of large deflections on the structural response was not assessed.
- Stress is not evaluated and there is no allowable stress in the dataset, so no factor of safety is given.
- The twist of the outer section is derived from a plane fit and is not validated against the Nastran reference.
- The displacement mapping used for plots leaves 583 of 61,552 result nodes of the expanded shell model (0.9%) unassigned (more than 0.05 m from any shell node); their nature was not identified. They are not used, and the mapped displacements agree with the printed ones to 4.1 mm at most.

Coupling procedure
- One correction step only; the aeroelastic loop is not iterated to convergence.
- The deformation applied to the CFD wing is first order: translation and twist per 1 m section with linear interpolation in y, no rotation about the x-axis, constant extrapolation beyond y = 29.24 m.
- The corrected loads are pulled back to the undeflected frame by moving the face centres and keeping the force vectors as computed, which is a small-displacement approximation.
- The CFD geometry is the undeflected (jig) shape and the frames of the CFD surface and the structural deck agree to about 0.1 m.
- The comparison of one-way and corrected loads with the dataset row mixes the difference in loads with the structural bias of -7.41%.
- Only one cruise condition (Mach 0.85, CL 0.5, 1 g) is analysed.

## 16. Conclusions

Within the assumptions above and on the L9 mesh:

1. The workflow from a public benchmark geometry through CFD, load mapping, a CalculiX shell analysis, a deformed-geometry CFD correction and a comparison table runs from the terminal with scripts and is reproducible.
2. Static deformation changes the aerodynamic load substantially. At the base angle of attack the deformed wing carries 31.8% less lift than the undeflected wing (CL 0.3409 against 0.4998). At equal lift (CL = 0.5) the corrected solution has a 31.9% smaller tip deflection (1.907 against 2.801 m) and a 23.5% smaller outer-section washout (-4.72 against -6.17 degrees), with the same total force to within 0.3%.
3. The one-way load is not an adequate load source for this wing: it corresponds to the undeflected shape, which the wing does not take in flight. One correction step moves the CFD toward the loaded shape and changes the structural response by far more than the mesh refinement from L9 to L10 does (+4.6% in tip deflection), although the mesh uncertainty of the absolute values (GCI of the order of 25% for the tip deflection) remains large.
4. The structural model is stiffer than the Nastran reference by 7.4% at the tip, and the CFD is a lower-fidelity Euler solution on meshes that are not converged. The absolute values (tip deflection, twist, loads) are therefore indicative, and the relative effect of the correction is the main result.
5. This is a verification study. No validation against experimental data was performed.

Natural next steps, none of which was done here: a viscous or higher-order CFD solution on a converged mesh sequence, a structural mesh study with a different solver or build, a geometrically nonlinear analysis, and iterating the correction to convergence.

## 17. Figures and videos

All figures are generated by scripts in `scripts/` and saved in `figures/`. The case-dependent figures are produced by `scripts/make_figures.sh`.

| File | Content | Script |
|---|---|---|
| `mesh_study_p0p36.png` | Cl, Cd, Cm and tip deflection against cell count for L8, L9 and L10, with the Richardson extrapolation and the GCI as error bar | `plot_mesh_study.py` |
| `force_history_mesh_study.png`, `force_history_coupling.png` | force-coefficient histories of the mesh-study runs and of the runs used in the coupling comparison (load-source iterations marked) | `plot_force_histories.py` |
| `contour_cp_surface_cl0p5.png` | Cp on the upper and lower wing surface: one-way L9 (CL 0.4998) and the corrected CL 0.5 pressure | `contour_cp_surface.py` |
| `contour_mach_slices.png` | Mach number on the planes y = 6, 14 and 22 m with the sonic line: undeflected alpha 0.36 degrees (iteration 1500), deflected alpha 0.36 degrees (1500), deflected alpha 2.07 degrees (1300) | `contour_mach_slices.py` |
| `contour_displacement_cl0p5.png` | vertical displacement of the shell model, planform and rear view at true scale: dataset loads, one-way L9, corrected CL 0.5 | `contour_displacement.py` |
| `oneway_coupling_p0p36.png` | spanwise load, deflection and twist: dataset loads, one-way L9, one-way L10 | `plot_oneway_results.py` |
| `coupling_correction_cl0p5.png` | spanwise load, deflection and twist: dataset loads, one-way L9, corrected CL 0.5 | `plot_oneway_results.py` |
| `cfd_cp_sections_<tag>.png` | chordwise Cp at several span stations of one CFD state; the tag is the case and the iteration: `u_p0p36_t1500` (undeflected L9), `u_p0p5_t1500` (undeflected, alpha 0.5 degrees), `u_L10_p0p36_t1148` (L10), `u_defl_a2p07_t1300` and `_t1500` (deflected, alpha 2.07 degrees), `b_m0p8_t1700` and `tvd600` (earlier TVD runs of the angle-of-attack search) | `plot_cp_sections.py <time> <case_dir> <tag>` |
| `cfd_section_overlay_b_m0p8_t1700.png` | wing sections with the cells at the pressure and temperature limiter bounds marked, earlier TVD run at alpha -0.8 degrees (evidence behind decision D18) | `plot_section_overlay.py` |

Videos (visualisations, not simulations):
- `flow_tracers_y06_small.mp4` (12 s, 1920 x 1080, 6.8 MB): tracer particles released upstream of the plane y = 6 m and advected through the steady CFD field, for the undeflected wing (alpha 0.36 degrees, iteration 1500) and the deflected wing (alpha 2.07 degrees, iteration 1300). The field is steady; the animation shows streamline "timelines" in 50-fold slow motion and is not a time-accurate solution. The full-quality version (42 MB) is not stored in the repository (`scripts/make_particle_video.py`).
- `wing_bending_cl0p5.mp4` (12 s): the shell model deforming from the undeflected shape to the loaded shape for the one-way L9 and the corrected CL 0.5 cases, at true scale (no exaggeration), with the undeflected mesh as a grey reference. It scales the linear static displacement with a load factor and is not a transient (`scripts/make_bending_video.py`).

Reading the figures:
- The corrected twist curve nearly coincides with the curve computed with the dataset loads (outer-section twist -4.72 degrees against -4.62 degrees) while the one-way twist is about 1.5 degrees more nose-down. The corrected load is higher inboard and lower outboard than the one-way load. The twist is not validated against the reference.
- The corrected-pressure Cp figure shows a small, localised suction pocket near y = 6 m (Cp minimum -1.28 against -1.18 for the one-way case), and the Mach slices show a local Mach number of 1.62 at y = 6 m for the deflected alpha 2.07 degree state. This is the growing suction region described in Section 12.2. Its cause was not established and it is part of the load source of the corrected solution (see Section 15).

## 18. Reproducibility

Environment: Windows with WSL2 and Ubuntu 24.04; OpenFOAM 11; CalculiX 2.23 (conda environment `ccx223`); Gmsh 4.12 with its Python API and the OpenCASCADE kernel; Python 3 with NumPy, SciPy, Pandas, Matplotlib, PyVista 0.48, imageio and ffmpeg. All commands run from the project directory.

Steps (script usage lines are in the script headers):
1. Data: `bash scripts/download_ucrm.sh` downloads and verifies the uCRM-9 dataset into `data/` (not committed).
2. Geometry: `python3 scripts/build_wing_geometry.py <wbt.igs>` and `python3 scripts/build_wing_stl.py` produce the closed wing surface `geometry/wing_full.stl` (not committed).
3. Structure-only decks and verification against the reference: `RUN=1 bash scripts/build_structure_decks.sh`; `scripts/post_static.py`, `scripts/post_modal.py` and `scripts/check_equilibrium.py` evaluate them (`results/structure_only_coarse.txt`, `results/static_equilibrium_coarse.txt`).
4. CFD meshes: the meshes are built from the base case `simulation/cfd_v1/feas` with `scripts/snappy_level_sweep.sh` (`QTIGHT`, `SEDX`, `TAG` and `BASECASE` select the variant). The production meshes are the `_novolq` variants (`feas_L8_novolq`, `feas_L9_novolq`, `feas_L10_novolq`, and `feas_L9_novolq_defl` for the deformed wing). The exact command lines were not saved; the settings of the meshes actually used and the mesh behind each run are recorded in `results/production_meshes.txt` (`scripts/record_production_meshes.sh`). The wing surface and edge mesh of the base case are copies of `geometry/wing_full.stl` and are not committed.
5. CFD runs: `bash scripts/cfd_run_upwind.sh <alpha_deg> <tag> [iterations] [processes]` (set `MESHSRC` to use another mesh), then `bash scripts/cfd_postprocess.sh <case_dir> <time>` to reconstruct and export VTK. The earlier two-stage angle-of-attack search used `scripts/cfd_run_alpha.sh`. Run files that differ from the template case are committed under `simulation/cfd_v1/runs/`.
6. Structure runs: `bash scripts/run_ccx.sh <dir> <job>`. Loads: `scripts/map_cfd_loads.py` (usage in its header; the option sets used for the production decks were not recorded in the decks). Correction step: `scripts/section_kinematics.py`, `scripts/deform_wing_stl.py`, a new mesh and CFD run for the deformed surface, `scripts/map_cfd_loads.py --pullback`, and `scripts/blend_wing_vtk.py` for the lift trim.
7. Tables and checks: `scripts/alpha_cl_table.py`, `scripts/mesh_study.py`, `scripts/coupling_table.py`, `scripts/check_limiter_fields.py`, `scripts/record_production_meshes.sh`.
8. Figures and videos: `scripts/plot_mesh_study.py`, `scripts/plot_force_histories.py`, `bash scripts/make_figures.sh`, `scripts/make_particle_video.py`, `scripts/make_bending_video.py`.

Not in the repository (regenerated by the scripts): the dataset, geometry files, meshes, CFD run folders (1.0 to 1.9 GB each), CalculiX decks and nodal loads, section tables and the full-quality particle video. Some intermediate steps (the exact option sets of the load mapping and of the mesh generation) were run interactively and are described here but not recorded as saved command lines.

## 19. References

1. Brooks, T. R., Kenway, G. K. W., and Martins, J. R. R. A., "Benchmark Aerostructural Models for the Study of Transonic Aircraft Wings", AIAA Journal 56(7), 2018, pp. 2840-2855. DOI 10.2514/1.J056603.
2. Brooks, T. R., Kenway, G. K. W., and Martins, J. R. R. A., "uCRM: undeflected Common Research Model", Mendeley Data, version 1, 2019, DOI 10.17632/gpk4zn73xn.1, licence CC BY 4.0.
3. Vassberg, J. C., DeHaan, M. A., Rivers, S. M., and Wahls, R. A., "Development of a Common Research Model for Applied CFD Validation Studies", AIAA Paper 2008-6919, 26th AIAA Applied Aerodynamics Conference, Honolulu, 2008.
4. Celik, I. B., Ghia, U., Roache, P. J., Freitas, C. J., Coleman, H., and Raad, P. E., "Procedure for Estimation and Reporting of Uncertainty Due to Discretization in CFD Applications", Journal of Fluids Engineering 130(7), 078001, 2008.
5. OpenFOAM 11 (OpenFOAM Foundation) and CalculiX 2.23 (G. Dhondt and K. Wittig), Gmsh 4.12 (C. Geuzaine and J.-F. Remacle).
