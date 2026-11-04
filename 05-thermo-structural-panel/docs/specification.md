# Project 05 specification: thermo-structural aerospace panel

Status: draft, written 2026-11-01 before any simulation. Tolerances below are set before the runs and will not be adjusted to fit results.

## Research question

How do a transient non-uniform heating pulse, panel thickness and edge restraint control the thermal stress and out-of-plane deformation (thermal buckling) of a thin aerospace skin panel?

## Locked decisions

| ID | Decision |
|---|---|
| P05-D1 | Flat monolithic rectangular plate, about 500 x 500 mm, 1 to 3 mm thick (stiffened and curved panels declined) |
| P05-D2 | Ti-6Al-4V, constant properties (`scripts/material.py`, cited typical vendor values) |
| P05-D3 | Baseline loading: transient heating pulse (surface heat flux ramps up, then soaks) |
| P05-D4 | Baseline edges: simply supported in a frame (out-of-plane displacement fixed, in-plane restrained) |
| P05-D5 | Thermal buckling: geometrically nonlinear analysis (NLGEOM) with a small imperfection seeded by the linear buckling mode |

Open decisions (to be asked before use): heat-pulse definition (peak flux, ramp and soak durations, spatial shape), back-face thermal condition, baseline thickness, final panel size.

## Analysis chain

1. Gmsh (Python API) builds a structured hexahedral solid mesh, quadratic elements, several elements through the thickness.
2. CalculiX 2.23 heat transfer: transient conduction with the surface flux.
3. CalculiX static step reading the temperature field as a thermal load, one-way (sequential). Linear for the verification cases, geometrically nonlinear with an imperfection for the main study.
4. Python extracts results, compares with closed-form solutions, runs the parameter sweep and produces PyVista and Matplotlib figures.

Units: SI (m, kg, s, N, Pa, W, J). Temperatures in degC, stress-free reference 20 degC.

## Verification plan (tolerances fixed in advance)

| ID | Case | Closed-form reference | Tolerance |
|---|---|---|---|
| V1 | 1D steady conduction through the thickness, fixed face temperatures | T(z) linear between the face temperatures | 0.5 % of the temperature difference |
| V2 | Uniform temperature rise, minimally constrained (free expansion) | edge displacement = alpha dT L, von Mises stress about zero | 1 % displacement; stress below 0.1 % of E alpha dT |
| V3 | Uniform temperature rise, in-plane restrained, below the buckling threshold (linear) | sigma_x = sigma_y = -E alpha dT / (1 - nu) | 1 % |
| V4 | Linear through-thickness gradient, minimally constrained | curvature kappa = alpha dT_g / h, stress about zero; corner-to-centre deflection = kappa L^2 / 4 | 2 % |
| V5 | 1D transient slab with a surface flux, insulated back face | series solution (Carslaw and Jaeger), reference form fixed when the case is built | 2 % of the surface temperature rise at the chosen times; time-step convergence reported |
| V6 | Linear buckling of the restrained plate under uniform heating | dT_cr = 2 pi^2 h^2 / (12 (1 + nu) alpha b^2) (square, simply supported, equal biaxial compression, thin plate; amended 2026-11-01, see Amendments) | 5 % |
| V7 | Nonlinear post-buckling response just above dT_cr | single-mode imperfect-plate curve, corrected in Amendment 3 (Amendment 2 criterion failed, on record) | 5 % (Amendment 3) |

Mesh convergence of temperature, stress and deflection on three meshes, with the grid convergence index (Celik et al. 2008), as in Project 04.

Validation: none claimed. All cases above are verification (closed-form solutions of the mathematical model). Validation would need published experimental or benchmark data for a comparable panel; if a suitable reference is found it will be added and the claim will be limited to what it supports.

## Planned limitations (to be stated in the README)

- Constant material properties: real conductivity and CTE rise with temperature, modulus and yield fall; results hold for a moderate temperature rise only.
- Property values are typical vendor values, not design allowables; stress is compared with typical room-temperature yield as a ratio, not a factor of safety, and this is unconservative when the plate is hot.
- Sequential one-way thermal-structural coupling; no radiation or convection from the faces unless added by an explicit decision.
- Idealised flat plate and idealised frame support; no stiffeners, curvature, rivet or fastener detail.
- Post-buckling results depend on the imperfection amplitude and shape.

## Amendments

### Amendment 1 (2026-11-01): reference formula for V6 corrected before the V6 run

The V6 row originally gave dT_cr = 4 pi^2 h^2 / (12 (1 + nu) alpha b^2). That coefficient belongs to a plate compressed in one direction only. With the in-plane restraint of this project (all four edges restrained) the uniform temperature rise produces equal compressive membrane stress in both directions, N = E alpha dT h / (1 - nu). A simply supported square plate under equal biaxial compression buckles in the one-half-wave mode (m = n = 1) at N_cr = 2 pi^2 D / b^2, with D = E h^3 / (12 (1 - nu^2)). Equating the two gives

dT_cr = 2 pi^2 h^2 / (12 (1 + nu) alpha b^2)

which is 0.575 K, 2.301 K and 5.177 K for h = 1, 2 and 3 mm at b = 500 mm (Ti-6Al-4V constant properties of scripts/material.py). The tolerance (5 %) is unchanged. The correction was made on review of the derivation, before any buckling analysis was run, so it is not a re-tuning to a result. Shear-deformation effects on a thin plate are of order (h/b)^2 = 1.6e-5 and are neglected in the reference.

### Amendment 2 (2026-11-01): definition of V7 (nonlinear post-buckling), fixed before any nonlinear run

Decisions taken by the project owner: the verification reference is the analytic single-mode curve below (a shell-element cross-check and a digitised published curve were declined); the imperfection is 0.1 h as the reference case with a sensitivity sweep at 0.01 h and 0.5 h.

Model. Plate 0.5 x 0.5 x 0.002 m, Ti-6Al-4V (constant properties of scripts/material.py), mesh 20 x 20 x 2 C3D20R (converged for the first buckling mode in V6c; one load level is repeated on 40 x 40 x 2). Supports: the membrane restraint of V3c and out-of-plane support on the mid-thickness edge lines. Load: uniform temperature rise dT applied as *TEMPERATURE in one geometrically nonlinear static step (*STATIC with NLGEOM), ramped linearly from 0 to 2.5 dT_cr, where dT_cr = 2.30101 K is the closed form of Amendment 1. Imperfection: every node is shifted in z by A sin(pi x/L) sin(pi y/L) (edges do not move; the finite-element first mode has this shape, V6: R^2 = 1.00000). Quantity compared: total centre deflection W = A + u_z at the centre mid-plane node.

Reference. With w = W sin(pi x/L) sin(pi y/L), initial imperfection A in the same shape, equal biaxial membrane compression N and immovable edges: equilibrium of the mode gives N = N_cr (1 - A/W), N_cr = 2 pi^2 D / L^2; zero mean in-plane edge displacement gives alpha dT - N (1 - nu)/(E h) = (pi^2 / (8 L^2)) (W^2 - A^2). With tau = dT/dT_cr, omega = W/h, a = A/h and c = 4/(3(1 + nu)):

omega^3 - [a^2 + c (tau - 1)] omega - c a = 0,

whose unique positive root is the reference. The derivation is the author's own and was not cross-checked against a published formula or table (none was found). It assumes one mode and a uniform membrane force (a Berger-type approximation), so the reference is approximate; this is why the tolerance is 10 %. The reference values are computed by scripts/ref_postbuckling.py and stored in results/v7_reference.txt before the first nonlinear run. The 3D finite-element linear threshold is 0.3 % below dT_cr (V6c, V6d); the reference keeps the closed-form dT_cr.

Acceptance criteria (fixed before the run, not to be changed after seeing results).
1. The reference case (A = 0.1 h) runs to tau = 2.5 without solver errors or warnings.
2. At tau = 1.2, 1.5 and 2.0, |W_FE / W_ref - 1| <= 10 % (A = 0.1 h).
3. Mesh check: at tau = 1.5 the 40 x 40 x 2 mesh changes W by at most 2 % relative to 20 x 20 x 2.
Reported without pass/fail: tau = 0.5, 0.8, 1.0 and 2.5 for the reference case, and the whole curves for A = 0.01 h and 0.5 h (imperfection sensitivity), with the corresponding reference curves.
A failed criterion is recorded as a failure and investigated; the criteria are not re-tuned.

Limitations of the verification. It verifies the nonlinear thermal load path, the restraint and the imperfection implementation against an approximate analytic relation; it is not validation against experiment. Properties are temperature independent.

### Amendment 3 (2026-11-02): V7 criterion of Amendment 2 failed; reference corrected (post hoc); blind checks

Result recorded. The first V7 run (A = 0.1 h, mesh 20 x 20 x 2, 50 increments, no warnings) gave W/h = 0.5449, 0.6985 and 0.9093 at tau = 1.2, 1.5 and 2.0, against 0.6118, 0.7980 and 1.0525 from the reference of Amendment 2: -10.9 %, -12.5 % and -13.6 %, outside the 10 % tolerance. V7 is recorded as failed against the Amendment 2 reference (results/v7_a0p1_amendment2_FAIL.txt); that reference and its tolerance are not changed.

Diagnosis. The finite-element mean edge force exceeds N_cr after buckling (N/N_cr = 1.02 to 1.30), which a uniform membrane force cannot produce. The Amendment 2 formula neglects the non-uniform membrane stress induced by the deflection.

Corrected reference. With w = W S and w0 = A S, S = sin(pi x/L) sin(pi y/L), the von Karman compatibility equation with the imperfection gives the particular membrane forces N_x = -kappa cos(2 pi y/L), N_y = -kappa cos(2 pi x/L), kappa = E h pi^2 (W^2 - A^2)/(8 L^2), in addition to a uniform compression N. This field has zero normal edge displacement at every edge point and zero edge shear traction, which are the in-plane edge conditions of the model (normal displacement restrained, tangential displacement free). Galerkin projection of the equilibrium equation on S gives N = N_cr (1 - A/W) + kappa/2, and zero mean edge displacement gives alpha dT - pi^2 (W^2 - A^2)/(8 L^2) = N (1 - nu)/(E h). Eliminating N:

omega^3 - [a^2 + c' (tau - 1)] omega - c' a = 0, with c' = 8/(3 (1 + nu)(3 - nu)) = 0.75094 for nu = 0.33,
N/N_cr = (1 - a/omega) + 3 (1 - nu^2)(omega^2 - a^2)/8,

where omega = W/h, a = A/h, tau = dT/dT_cr, and N is the mean edge force per unit width. The values are produced by scripts/ref_postbuckling.py corrected and stored in results/v7_reference_corrected.txt before the A = 0.01 h and 0.5 h runs.

Status of the corrected reference. It was derived after the first V7 result was known. Agreement of the reference case (A = 0.1 h) with it is therefore a post hoc consistency check and not an independent verification. It is a single-mode truncation and has not been compared with a publication (none was found). The A = 0.01 h and 0.5 h runs, whose predictions were stored before they were run, are the blind check.

Acceptance criteria (replace criterion 2 of Amendment 2; fixed before the sweep runs; not to be changed after seeing results).
1. Each of the A = 0.01 h, 0.1 h and 0.5 h runs reaches tau = 2.5 without solver errors or warnings.
2. Blind: for A = 0.01 h and A = 0.5 h, |W_FE / W_ref - 1| <= 5 % at tau = 1.2, 1.5 and 2.0 (corrected reference).
3. Post hoc: the same criterion for A = 0.1 h, labelled post hoc.
4. Mesh check (unchanged from Amendment 2): at tau = 1.5 and A = 0.1 h the 40 x 40 x 2 mesh changes W by at most 2 % relative to 20 x 20 x 2.
Reported without pass/fail: tau = 0.5, 0.8, 1.0 and 2.5, and the membrane force N/N_cr against the corrected relation.
V7 is reported as passed only if criteria 1 to 4 are met, together with the recorded failure against the Amendment 2 reference.

### Amendment 4 (2026-11-02): baseline loading definition (P05-D8)

Fixed by the project owner before any baseline run:
- Heating: Gaussian hot spot on the top face, q(x,y,t) = q0 f(t) exp(-r^2/(2 s^2)), centred at (L/2, L/2), s = L/6 = 83.33 mm, q0 = 30 kW/m^2.
- Time history: f(t) rises linearly from 0 to 1 over 5 s, then stays 1 for 25 s (total 30 s). No cool-down.
- Thermal boundaries: back face and four edges insulated. No convection, no radiation.
- Baseline panel: 500 x 500 mm, thickness 2 mm. Thicknesses 1 mm and 3 mm belong to the parametric study.
- Structural boundary conditions (P05-D4) and imperfection (P05-D7) are unchanged.
- Hand estimate, not a result: with no lateral conduction the peak rise is about q0 (27.5 s) / (rho c h) = 177 K. Conduction will lower it, so the verified value comes from the simulation.

### Amendment 5 (2026-11-02): V8a, thermal step verification (fixed before the run)

Case: baseline Gaussian pulse of Amendment 4 on the 20 x 20 x 2 and 40 x 40 x 2 C3D20R meshes (heat transfer, DIRECT, backward Euler), dt = 0.5 s (both meshes) and dt = 0.25 s (20 x 20).
1. Energy conservation (insulated plate): the volume-averaged temperature rise at t = 5, 10, ..., 30 s agrees within 0.05 % with sum_n P f(t_n) dt / (rho c V), where P is the total power printed by make_flux_top.py for that mesh and f(t_n) the ramp evaluated at the end of each increment (the exact algebraic result of the time scheme). The continuous-time value (rise 30.73 K at 30 s) is reported; the scheme adds dt/2 s of full-power equivalent for the ramp, so it is 0.9 % higher at dt = 0.5 s and 0.45 % higher at dt = 0.25 s.
2. Time step: the peak temperature rise at 30 s on the 20 x 20 mesh changes by less than 1 % of the rise between dt = 0.5 s and dt = 0.25 s.
3. Symmetry and through-thickness: the final field is symmetric under x <-> y and x -> L - x within 0.01 K, and the top-minus-bottom temperature difference at the plate centre at 30 s is within 10 % of q0 h / (2 k) = 4.48 K (quasi-steady conduction through the thickness).
Reported without pass/fail: peak temperature, 20 x 20 versus 40 x 40.
Withdrawn before any run: the continuous-time energy reference for criterion 1 and a "top minus bottom below 2 K" limit, both of which were wrong by hand calculation.

### Amendment 6 (2026-11-02): V8a attempt 1 and element choice for thermal steps

- Attempt 1 (thermal run on the C3D20R mesh): criteria 1 and 2 passed (energy deviation 2e-5 %, peak change with dt 0.43 %). Criterion 3 FAILED: top-minus-bottom at the plate centre was 1.27 K (20 x 20) and 3.92 K (40 x 40) against 4.48 K +/- 10 %. Along y = L/2 on the 20 x 20 mesh the difference was 1.3 K on element-corner columns and 1.8 to 3.6 K on mid-side columns (quasi-steady theory 2.2 to 4.4 K). Record: results/v8a_attempt1_criterion3_FAIL.txt.
- Diagnostic (not a verification): the same 20 x 20 x 2 nodes, flux and deck with full-integration C3D20 gave a smooth field (top minus bottom 4.20 K at the centre, 6 % below theory; peak at the centre node). Changing only the integration rule removes the pattern.
- Decision (implementation detail): thermal steps use C3D20 (full integration); structural steps keep C3D20R (verified in V3 to V7). The meshes have identical nodes and numbering; V8b checks the hand-over.
- V8a attempt 2: same cases and unchanged criteria 1 to 3, on C3D20 meshes v8a_n20_z2 and v8a_n40_z2 (job names v8a2_*). Criterion 3 is not blind: the diagnostic already gave 4.20 K on the 20 x 20 mesh.

### Amendment 7 (2026-11-02): V8b, hand-over of temperatures to the structural run (fixed before the run)

Thermal run on C3D20 mesh v8a_n20_z2 (dt = 0.5 s, results file .frd with every increment); linear static structural runs on the identical C3D20R mesh v6c_n20_z2 read it with *TEMPERATURE, FILE= and are held only against rigid-body motion (three corner nodes). Run A uses dt = 0.5 s, run B dt = 0.25 s.
1. Run A: nodal temperatures reported by the structural run at t = 5, 10, ..., 30 s (all nodes) agree with the thermal run to within 0.005 K (the .frd stores 5 significant digits).
2. Run B: at multiples of 0.5 s the centre-column temperatures agree with the thermal run to within 0.005 K.
3. Run B, other times from 0.75 s: report the maximum deviation from linear interpolation, previous frame and next frame of the thermal results, and name the rule CalculiX applies. No pass/fail; the rule decides how the baseline run may be stepped.
4. The thermal run reproduces V8a attempt 2 at the centre column at 5, 10, ..., 30 s within 0.005 K, and all three runs finish with zero warnings.
Not verified here: the structural response to the temperature field (the thermal expansion path was verified in V1 to V4).

### Amendment 8 (2026-11-02): V8b attempt 1 and deck change

- Attempt 1 (linear *STATIC step, two-value time line): the structural runs produced one temperature block each (t = 30 s) where 6 (run A) and 120 (run B) were expected, and a warning that the minimum increment (default 0) was changed. The script stopped before evaluating criteria. The one block available agreed with the thermal run to 5.0e-4 K over all nodes (hand-over syntax and node numbering confirmed at t = 30 s).
- Attempt 2: structural step is *STEP, NLGEOM with *STATIC dt, 30., 1.E-3, dt (fixed increments, as in the baseline). Criteria 1 to 4 of Amendment 7 are unchanged. The earlier "linear static" description in Amendment 7 is superseded by this line.

### Amendment 9 (2026-11-02): V8b attempt 2 and change of hand-over method

- Attempt 2 (NLGEOM step, fixed increments, *TEMPERATURE, FILE=): 0 warnings, 6 and 120 output blocks, thermal run reproduces V8a attempt 2 exactly, but criterion 1 FAILED (max |dT| 13.18 K) and criterion 2 FAILED (13.26 K). Record: results/v8b_attempt2_time_mapping_FAIL.txt.
- Diagnosis (scripts/diag_v8b.py, existing results only): structural temperatures equal T0 + (t/30)(T_last - T0) to 4.3e-4 K (run B, 120 times) and 5.0e-4 K (run A, all nodes). *TEMPERATURE, FILE= takes the last frame of the thermal results and the step ramps linearly to it; it does not follow a history inside one step.
- Attempt 3: one structural step (0.5 s, NLGEOM) per thermal frame, each with an explicit *TEMPERATURE field (all nodes) taken from the thermal .dat, so each step ramps linearly from the previous frame to the next. Run A uses one increment per step, run B two. Criteria 1 to 4 of Amendment 7 are unchanged; criterion 3 now includes every time that is not a multiple of 0.5 s (from 0.25 s).

Note on V8b attempt 3 (2026-11-02): run A produced 64 output blocks instead of 60. The structural step 8.0 to 8.5 s diverged at the 0.5 s increment and CalculiX cut it back (extra blocks at 8.125, 8.25, 8.4375, 8.75 s); all 60 frame times are present. The script assertion "exactly 60 blocks" was replaced by "every frame time present, extra blocks reported". Criteria 1 to 4 are unchanged. The cutback is a mechanical-solution event, not verified or interpreted in V8b; it is a flag for the baseline run.

V8b result (2026-11-02, attempt 3): PASS. Structural temperatures equal the thermal history at all 60 frames (run A, all nodes, 0.000 K; run B, centre column, 0.000 K; limit 0.005 K) and are linear in time between frames (max deviation 5e-5 K; previous- and next-frame rules deviate by 1.69 K). The thermal run reproduces V8a attempt 2 exactly; zero warnings. Record: results/v8b_transfer.txt. Attempts 1 and 2 (Amendments 8 and 9) remain on record as failures. Not verified in V8b: the mechanical response to the temperature field.

### Amendment 10 (2026-11-03): baseline run definition and acceptance (fixed before the run)

Baseline: 500 x 500 x 2 mm Ti-6Al-4V (P05-D2), simply supported frame with membrane restraint (P05-D4, as in V3c, V6, V7), seeded imperfection A = 0.1 h in the first-mode shape (P05-D7) directed toward +z, the heated face (the through-thickness temperature gradient also bows the panel toward its hot side); the -z sign is a sensitivity case for the parametric study. Gaussian pulse of Amendment 4. Thermal run on C3D20 (V8a), structural run on C3D20R, 60 NLGEOM steps of 0.5 s with one thermal frame per step (V8b), automatic cutback allowed (minimum increment 1e-3 s).
Acceptance of the run (code correctness only, no physical claim):
 B1. The run reaches t = 30 s with every frame time present in the output and zero warnings; cutbacks are reported, not penalised.
 B2. The centre-column temperatures reported by the structural run equal the thermal history at all 60 frames within 0.005 K.
 B3. Global equilibrium: each component of the total reaction over all nodes, at every printed increment, is at most 1e-3 of the largest magnitude of the x reaction on face x = 0 over the run.
Reported without pass/fail: centre and maximum deflection histories, edge membrane force history.
Elastic sanity (maximum von Mises stress against the 880 MPa yield stress) is evaluated in the next stage from the saved stresses. If the yield stress is reached, I stop and ask before continuing, because that changes the material model (P05-D2).

Baseline run result (2026-11-03): B1, B2 and B3 pass (no cutbacks, no warnings; B2 0.000 K; B3 ratio 7.9e-15). Note on B3: the sum of all nodal reactions vanishes identically when the element force vectors are self-equilibrated, so B3 is nearly a tautology and not evidence of convergence; the evidence of convergence is that CalculiX completed every step without cutback. The parametric study adds a balance of the membrane forces on opposite faces. Results (2 mm, imperfection 0.1 h toward +z): centre deflection W/h rises monotonically from 0.26 (2.5 s) to 4.16 (30 s, about 8.3 mm); the maximum deflection is at the plate centre at all reported times; |Rx| on face x = 0 is 1.04 kN at 2.5 s, 1.77 to 1.92 kN afterwards (1917.7 N at 30 s). Elastic stress check pending.

### Amendment 11 (2026-11-03): elastic stress check rule (fixed before the run)

Maximum von Mises stress over all integration points of the baseline run, at every printed increment, multiplied by 1.267 (= 1/0.789; the outermost integration plane of the two quadratic layers lies 79 % of the way from mid-thickness to the surface, so this bounds a bending-dominated surface stress). If the bound reaches 660 MPa (75 % of the 880 MPa yield stress, allowing for the lower yield stress at the roughly 200 C of the hot region), I stop and ask before continuing, because the constant-property elastic model (P05-D2) would then be in question. Below 660 MPa the elastic model stands for the baseline. This is a sanity check on the model, not a verification of the stress values.

Baseline stress check result (2026-11-03): maximum von Mises stress 89.6 MPa at t = 30 s (bound 113.6 MPa = 12.9 % of yield, 17 % of the 660 MPa stop bound), in the bottom layer next to the plate centre, equal-biaxial compression (min(sxx, syy) = -89.6 MPa), max |szz| 0.11 MPa. Stop rule of Amendment 11 not triggered; the elastic model stands for the 2 mm baseline. Record: results/p05_baseline_stress.txt. Hand estimate for scale (infinite plate, Gaussian hot spot, -E alpha dT0 / 2 with dT0 of about 175 K): 85 MPa; not a verification.

### Amendment 12 (2026-11-03): parametric study definition (P05-D9, P05-D10), fixed before generation

Cases (in-plane mesh 20 x 20, two C3D20R layers through the thickness for the structure, C3D20 for the thermal run, the Gaussian pulse of Amendment 4 with q0 = 30 kW/m2 for every thickness, imperfection 0.1 h toward +z unless stated):
 - 9 cases: h = 1, 2, 3 mm x edge condition {SS-restrained (baseline), SS in-plane free, clamped}.
 - 3 imperfection cases at h = 2 mm, SS-restrained: +0.01 h, +0.5 h, -0.1 h.
 - 1 mesh case: baseline at 40 x 40.
Edge conditions: SS-restrained as in V3c, V6, V7 and the baseline (frame equations plus uz = 0 on the mid-thickness edge lines). SS in-plane free: uz = 0 on the mid-thickness edge lines only, no frame equations; in-plane rigid-body modes removed by ux, uy at the plate-centre mid-plane node and uy at the mid-plane node at (L, W/2); with symmetric heating these supports should carry no reaction (checked). Clamped: normal in-plane displacement and uz fixed at every node of the four edge faces, tangential displacement free.
Gate: before any new case, the generalised case script reproduces the baseline (h = 2 mm, SS-restrained, +0.1 h): centre deflection and Rx histories equal to the baseline histories within 1e-6 relative.
Per case: peak temperature, deflection histories (centre and maximum), edge force, maximum von Mises stress with the 1.267 bound and the stop rule of Amendment 11 (if triggered for a case, I stop and ask), solver warnings, cutbacks. Limitation recorded in advance: at h = 1 mm the peak temperature is about 370 C, beyond the range in which constant material properties are justified (about 200 to 300 C).

### Amendment 13 (2026-11-03): per-case checks of the parametric study (fixed before generation)

Every case (scripts/run_case.py) must pass, in addition to the stop rule of Amendment 11:
 G1. Thermal run: zero warnings, and the volume-averaged temperature rise at t = 5, 10, ..., 30 s within 0.05 % of the discrete energy reference (the V8a criterion 1, repeated for each thickness and mesh).
 G2. Structural run: zero warnings, reaches t = 30 s, every 0.5 s frame time present in the output (cutbacks are reported, not penalised).
 G3. Hand-over: centre-column temperatures of the structural run equal the thermal history within 0.005 K at all 60 frames.
 G4. Restraint balance. SS-restrained and clamped: |Rx(face x = 0) + Rx(face x = L)| <= 1e-3 max|Rx(face x = 0)| at every printed increment (symmetric loading, equal and opposite edge forces). SS in-plane free: the reactions at the three rigid-body supports (centre node ux, uy; edge node uy) are at most 1e-4 E alpha h L (T_peak - 20 C) at every increment, since symmetric heating should not load them.
 G5. Stop rule: 1.267 x max von Mises < 660 MPa at every increment; otherwise stop and ask.
The sum of all nodal reactions is no longer used as a check (Baseline note).
Regression gate (Amendment 12) applies to the first use of the script: maximum difference of W_centre/h and of Rx(face x = 0) against results/p05_baseline_history.csv, relative to the maximum of the baseline values, at most 1e-6.

Regression gate result (2026-11-03): scripts/run_case.py --h 2 --bc rest --imp 0.1 --regress reproduced the baseline: relative difference 1.2e-7 in W_centre/h and 2.1e-8 in Rx(face x = 0) (limit 1e-6). These differences equal the rounding of the stored baseline CSV (6 and 4 decimals), so the two computations are identical within that precision. G1-G5 passed (energy deviation 1.77e-7, no cutbacks, face balance below the 6-digit print resolution, maximum von Mises 89.60 MPa, bound 113.57 MPa). Note: G3 in this run compared only the top-centre temperature with the stored thermal history because the thermal .dat is deleted after use; the full-field hand-over is the V8b result.

### Amendment 14 (2026-11-03): mesh-refinement check at 1 mm (fixed before the run)

Reason: the 1 mm case (20 x 20 in-plane elements, 2 through the thickness, in-plane element about 50 times the layer thickness) gave peak von Mises 62.89 MPa, lower than at 2 and 3 mm although the temperature rise is twice as large. Mesh sensitivity was checked at 2 mm in the earlier verification stages, not at 1 mm. The study matrix therefore gets a second mesh case: h = 1 mm, simply supported restrained, imperfection +0.1 h, 40 x 40 in-plane elements, same thermal and structural settings (matrix now 14 cases).
Criteria (comparison of the 20 x 20 result with the 40 x 40 result at t = 5, 10, 15, 20, 25 and 30 s; two meshes show sensitivity, not an asymptotic convergence order):
 M1. W_centre/h: relative difference <= 3 % at every listed time.
 M2. Peak von Mises: relative difference <= 10 % at t = 30 s.
 M3. Edge force |Rx(x = 0)|: relative difference <= 10 % at t = 30 s.
Consequence if M1-M3 pass: the 20 x 20 mesh is adequate for the 1 mm cases and the difference is reported as the mesh uncertainty. If any fails: the 1 mm cases are run on the 40 x 40 mesh (a further refinement check is then needed) and the failure stays on record.
All G1-G5 checks of Amendment 13 apply to the 40 x 40 run too.

Amendment 14 result (2026-11-03): 1 mm, simply supported restrained, imperfection +0.1 h; 20 x 20 against 40 x 40 in-plane elements (2 layers through the thickness in both). M1: relative difference in W_centre/h 0.133, 0.119, 0.115, 0.112, 0.109, 0.108 % at t = 5, 10, 15, 20, 25, 30 s (limit 3 %); M2: peak von Mises 62.885 against 63.442 MPa, 0.89 % (limit 10 %); M3: edge force 705.5 against 682.1 N, 3.3 % (limit 10 %). All three pass; G1-G5 pass for the 40 x 40 run (6 cutbacks, no warnings). The 20 x 20 mesh is adequate for the 1 mm cases; the differences above are the mesh uncertainty. Not covered: through-thickness refinement (2 layers kept). Thermal field: the peak temperature at 30 s is 369.18 C on 20 x 20 and 373.79 C on 40 x 40 (+4.6 K, 1.3 % of the rise; 3.9 % of the rise at 5 s at the top centre), so the reported peak temperatures carry a mesh uncertainty of a few percent while the structural response is insensitive. Heat input: 1302.19 W (20 x 20) against 1302.00 W (40 x 40).

Imperfection group result (2026-11-03): 2 mm, simply supported restrained, A/h = -0.1, +0.01, +0.1 (baseline), +0.5; all G1-G5 pass (stop rule not triggered). At 30 s: W_centre/h 4.1591, 4.1601, 4.1633, 4.2036 (spread 1.1 %); peak von Mises 91.93, 90.64, 89.60, 85.11 MPa; |Rx(x = 0)| 2069, 1985, 1918, 1649 N; W_centre/h at 5 s 0.862, 0.912, 0.957, 1.187. Cutbacks: 4, 3, 0, 0. The imperfection against the thermal bowing direction (-0.1 h) gives the highest stress and edge force and is the conservative bound. For -0.1 h, W_centre/h stays negative until 3.5 s (-0.281) and is +0.370 at 4.0 s: the panel snaps through to the hot-face side inside one 0.5 s increment, with no cutback. The quasi-static NLGEOM solution does not resolve the unstable path across the snap; the values before and after are valid equilibrium states, the path between is not resolved (limitation). The generator self-check was corrected to compare with |A| (patched; regenerated positive mesh byte-identical to the earlier one).

Thickness group result (2026-11-03): 3 mm and 1 mm, free and clamped, imperfection +0.1 h; all G1-G5 pass (stop rule not triggered). At 30 s (W_centre/h, peak von Mises, |Rx(x = 0)|): 3 mm free 1.6045, 86.79 MPa, 0 (largest rigid-body support reaction 3.3 N); 3 mm clamp 1.8154, 110.52 MPa, 18.2 kN; 1 mm free 12.589, 63.04 MPa, 0 (0.57 N); 1 mm clamp 11.648, 62.77 MPa, 2.10 kN. Cutbacks: 3 mm none, 1 mm free 1, 1 mm clamp 2. Edge-condition sensitivity grows with thickness: peak von Mises relative to the restrained case is +0.2 / -0.2 % (1 mm), -0.1 / +1.8 % (2 mm), -4.7 / +21.4 % (3 mm) for free / clamped; W_centre/h +4.0 / -3.7 %, -0.3 / -6.0 %, -21.7 / -11.4 %. The 3 mm clamped case has the largest stress of the study (110.52 MPa; bound 140.07 MPa, 15.9 % of yield); the location of its maximum is not yet known. The earlier statement of weak edge-condition sensitivity holds for 1-2 mm only.

### Amendment 15 (2026-11-03): mesh-refinement check at 2 mm (fixed before the run)

Same case and criteria as Amendment 14, applied to the baseline thickness: h = 2 mm, simply supported restrained, imperfection +0.1 h, 20 x 20 against 40 x 40 in-plane elements (2 layers through the thickness in both; reference results results/cases/h2_rest_a0p1_n20*). M1: relative difference in W_centre/h <= 3 % at t = 5, 10, 15, 20, 25, 30 s; M2: peak von Mises at 30 s <= 10 %; M3: |Rx(x = 0)| at 30 s <= 10 %. G1-G5 apply. Peak temperature differences are reported, not gated. Consequence if any fails: the 2 mm cases are rerun on the 40 x 40 mesh and the failure stays on record.

Amendment 15 result (2026-11-03): 2 mm, simply supported restrained, +0.1 h; 20 x 20 against 40 x 40 in-plane elements. M1: relative difference in W_centre/h 0.115, 0.121, 0.130, 0.132, 0.129, 0.127 % at t = 5, 10, 15, 20, 25, 30 s (limit 3 %); M2: peak von Mises 89.603 against 90.437 MPa, 0.93 % (limit 10 %); M3: edge force 1917.7 against 1881.8 N, 1.9 % (limit 10 %). All pass; G1-G5 pass for the 40 x 40 run (no cutbacks, no warnings). Peak temperature at 30 s 196.70 C (20 x 20) against 199.10 C (40 x 40): +2.4 K, 1.4 % of the rise (1 mm: +4.6 K, 1.3 %). Structural results are insensitive to the in-plane mesh; the peak temperature is not (about -1.4 % of the rise on the coarse mesh). Two meshes only: no convergence order is claimed.

Stress location (2026-11-03): reruns of the 3 mm restrained and 3 mm clamped cases with the output kept (scripts/run_case.py --keep-dat) reproduced the stored histories bit-identically. scripts/locate_max_stress.py: the largest von Mises stress of the last increment (t = 30 s, also the maximum over the pulse) is 110.517 MPa (clamped) and 91.042 MPa (restrained), in the four elements around the plate centre (centroid radius 0.018 m), as equibiaxial in-plane compression (sxx and syy about -110.5 and -91.1 MPa, all other components below 1 MPa). Maximum in the edge strip (within 0.05 m of an edge): 66.31 MPa (clamped), 34.23 MPa (restrained). The clamped panel's higher stress is therefore a hot-spot stress, not an edge stress. The lower/upper integration-point plane labels of the script are not verified (equal values appear at points the assumed ordering places at non-mirror-image positions); only the element-level location is used. The measured stress is below the fully constrained biaxial thermal stress E alpha dT/(1 - nu) in every case (ratios are computed by scripts/make_figures.py).

Stress layer finding (2026-11-03): scripts/make_field_figures.py (centroid height splits the 2 element layers) gives the largest element von Mises of the top / bottom layer at 30 s: 1 mm restrained 56.51 / 62.89 MPa; 2 mm restrained 61.71 / 89.60 MPa; 3 mm restrained 41.98 / 91.04 MPa; 3 mm clamped 66.31 / 110.52 MPa. The study maximum of every one of these cases lies in the bottom layer (the cold, insulated back face), at the plate centre. In the clamped 3 mm case the top-layer maximum (66.31 MPa) equals the edge-strip maximum. A first version of the field figure coloured each visible top face with its own element value and therefore showed only the top layer, which looked inconsistent with the stored peak stress (the clamped centre looked lighter than its edges); the figure now shows the maximum over both layers at each plan position, and every row is checked against the stored history (centre deflection, peak von Mises, peak temperature). Hypothesis, not tested: bowing toward the hot face adds bending tension to the hot face and bending compression to the cold face; bending and membrane stress were not separated.

### Amendment 16 (2026-11-04): consolidated record of decisions and verification stages not written into this specification earlier

Purpose. The decisions and stages below were taken or run between 2026-11-01 and 2026-11-03. They are recorded in the result files and the run scripts but were not written into this document. This amendment records them after the fact, from the stored result files. Nothing here changes a tolerance, a criterion or a result. The status line and the list of open decisions at the top of this file are history; where they differ from this amendment, this amendment and Amendments 4, 12 and 13 apply.

Decision index. P05-D1 to P05-D5: table at the top. P05-D6 and P05-D7: below. P05-D8: Amendment 4. P05-D9 and P05-D10: Amendment 12. The open decisions listed at the top (heat pulse, back-face condition, baseline thickness, panel size) were closed by P05-D8: 2 mm baseline, 500 x 500 mm panel, insulated back face.

P05-D6 (decided 2026-11-01, recorded late): the verification reference for V7 (nonlinear post-buckling) is the analytic single-mode imperfect-plate curve only, (W/h)^3 - (a^2 + c (tau - 1)) (W/h) - c a = 0 with c = 4/(3 (1 + nu)), a = A/h, tau = dT/dT_cr, tolerance 10 % on the centre deflection for 1 < tau <= 2. A shell-element cross-check and a digitised published curve were declined. The reference and the tolerance were later changed by Amendment 3 (corrected reference, 5 %); the Amendment 2 failure stays on record.

P05-D7 (decided 2026-11-01, recorded late): the imperfection has the shape of the first buckling mode with amplitude 0.1 h for the V7 reference, plus an imperfection-sensitivity sweep at 0.01 h and 0.5 h (three runs). The sign of the baseline imperfection is fixed in Amendment 10.

Frame support of the plate model (implementation of P05-D4, chosen 2026-11-01 as an implementation detail, recorded late). P05-D4 requires in-plane restrained edges. Three implementations of the in-plane restraint were tested in V3: uniform rise of 2 K from 20 C on the 0.5 x 0.5 x 0.002 m plate, mesh v_n40_z2, central-region stress against -E alpha dT/(1 - nu) = -2.921433 MPa, tolerance 1 % fixed before the run.
- V3a, pinning only the mid-thickness nodes of the four edges (results/v3a_midline_support.txt, as run first): mean deviation 0.951 %, maximum deviation 1.412 %, edge reaction 99.04 % of the exact value. FAIL (the criterion is the maximum deviation); kept on record.
- V3b, whole edge face restrained (results/v3b_full_edge.txt): maximum deviation 5.6e-6 % (round-off), edge reaction 100.00 %. PASS. This restraint also fixes the edge rotation, so it is not used for the plate model.
- V3c, membrane restraint (results/v3c_membrane.txt): the thickness-mean of the in-plane displacement normal to each edge face is set to zero on every through-thickness line with *EQUATION (scripts/frame_constraints.py); the edge can rotate. Maximum deviation 0.0343 %, mean 0.0169 %, edge reaction 99.98 %. PASS.
Slice study (scripts/slice_study.py, results/slice_study.txt; strip 0.5 x 0.002 x 0.002 m, refinement levels 2, 4, 8, 16, dT = 2 K): the mid-line pin error grows with refinement, 0.807, 1.055, 1.393, 1.735 % (increments of 0.25, 0.34, 0.34 points per halving of the element size); the membrane restraint error falls, 0.0641, 0.0329, 0.0160, 0.00787 % (about first order); the full-face restraint stays at 5.6e-6 %. The mid-line pin is therefore a mesh-dependent support, not a converged one, and was rejected. Cause: hypothesis, consistent with the slice but not isolated further: the pin acts as a line load at one height whose local effect does not vanish with refinement. Consequence: the membrane restraint is the frame support of V6, V7 and every later run.

V1 to V4 (run 2026-11-01, recorded late; results/v1_conduction.txt, v2_free_expansion.txt, v4_gradient.txt): V1, V2, V3b, V3c and V4 pass (V3a is the recorded failure above). V2: displacement errors 1.3e-14 %, maximum von Mises stress 7.0e-8 % of E alpha dT. V4: centre rise error 2.4e-4 % (tolerance 2 %).

V5 result (results/v5_transient.txt, scripts/run_v5.py): slab with constant top flux 50 kW/m2, T0 = 20 C, h = 2 mm (h^2/alpha = 1.391 s), 8 x 8 in-plane mesh, nz = 2, 4, 8, dt = 0.05, 0.025, 0.0125 s. Error measure: maximum nodal temperature error through the thickness relative to the surface rise at that time. Criterion cell: nz = 2, dt = 0.0125 s, tolerance 2 %. Result: 0.108 % at t = 0.5 s and 0.0039 % at t = 1 s. PASS. Time-step convergence at t = 0.5 s, nz = 8: 0.568, 0.277, 0.136 % for the three time steps (observed order 1.04 and 1.02: backward Euler, first order in time). The through-thickness effect is +0.027 points at dt = 0.05 s (nz = 8 against nz = 2). Differences at t = 2 s and 4 s are below 7e-5 % of the surface rise. Practical note: CalculiX allows 100 increments per step by default; the 4 s runs with dt = 0.0125 s need 320, so the decks set INC=2000 on *STEP.

V6a (record of a failed attempt, results/v6a_thermal_preload.txt): linear buckling with a thermal preload (*TEMPERATURE 20 -> 21 C in the *STATIC step, then *BUCKLE 6) returned factors of 4.2e7 to 4.7e8 against a closed form of 2.30101 K. The preload step itself is correct (mean sxx = syy = -1.46047e6 Pa against the expected -1.46072e6 Pa per K, edge reaction 1460.471 N). The thermal preload did not enter the buckling analysis in CalculiX 2.23; this agrees with forum reports, and the internal cause was not investigated. Consequence: V6 uses the mechanical equivalent of the preload (prescribed mean in-plane edge displacement alpha dT L/2 through dummy nodes in the membrane equations, scripts/frame_prestrain.py). The thermal load path is exercised in the nonlinear static steps (P05-D5).

V6b result (results/v6_buckling.txt, scripts/check_v6.py; mesh v_n40_z2): mode 1 at 2.29383 K against 2.30101 K (-0.312 %, tolerance 5 %); modes 2 to 6 within -0.154 to -0.347 %; mode-shape fits to sin(m pi x/L) sin(n pi y/L) with R^2 >= 0.99999; preload stress sxx 0.017 % from -E alpha/(1 - nu) per K. PASS.

V6c result (results/v6_convergence.txt, scripts/run_v6_conv.py): mode 1 at n = 10, 20, 40: 2.28916, 2.29368, 2.29383 K (nz = 2) and 2.28890, 2.29316, 2.29322 K (nz = 4); observed orders 4.96 and 6.03 and GCI(40) printed as 0.000 % (the two finest meshes differ by only 0.00015 K, so the order is indicative); through-thickness effect at n = 40: -0.026 %; edge reaction within 0.017 % of the exact value. The remaining -0.31 % gap to the closed form is not a mesh effect. PASS.

V6d result (results/v6_thickness.txt, scripts/run_v6_thick.py): reference preload 0.1 K, 20 x 20 x 2 C3D20R. dT_cr (FE / closed form) at h = 1, 2, 3, 5 mm: 0.57424 / 0.57525, 2.29368 / 2.30101, 5.15246 / 5.17727, 14.26429 / 14.38131 K; gaps -0.177, -0.319, -0.479, -0.814 %; mode 2 / mode 1 = 2.50185 to 2.50455 against 2.5. All within 5 %; log-log slope of the gap against h is 0.95, so the gap is about first order in h/L. Its cause (3D or shear effects, or the edge restraint) was not isolated. PASS.
V6d attempt 1 (results/v6d_attempt1_unit_preload_FAIL.txt, kept on record): with a 1 K reference preload the 1 mm plate returned a first factor of 1.43665 K (+149.7 %) and a degenerate pair (mode 2 / mode 1 = 1.00000). This equals the second mode group of attempt 2 (0.57424 K x 2.50185 = 1.43666 K), so the fundamental mode (0.575 K, below the shift value of 1) was not returned. The scripts/frame_prestrain.py header records that CalculiX returned no factor below its shift value of 1 in this case; the internal cause is a hypothesis, supported by the fix working. Change made: reference preload 0.1 K (scripts/frame_prestrain.py --dT), so that every factor of interest is above 1. No tolerance was changed.

Element-type limitation of the mesh generator: scripts/make_panel_mesh.py --etype C3D8 changes only the element label; it still writes the 20-node connectivity (mesh/v8x_n40_z4_c3d8_info.txt: "C3D8 (Gmsh type 17, 20 nodes)"). It therefore does not produce a first-order mesh, and no result obtained with that label is to be read as an 8-node result. The meshes named v8x_* are diagnostic scratch and are not reproduced by scripts/build_meshes.sh.

Reproducibility record (2026-11-04). The nine hand-written decks of V1, V2, V3a, V3b, V3c, V4, v6_buckling, v6_diag_prestress and v6b_prescribed are tracked in verification/decks/ (identical to the decks that produced the stored results, checked with cmp). scripts/build_meshes.sh builds the meshes that no run script builds (v_n40_z2 with its frame and prestrain files, the v6c_n20 and v6c_n40 bases, the four v7_* imperfect meshes, v8a_n20 and v8a_n40); its verify mode rebuilt them under scratch names and compared 19 files, all identical to the meshes of the stored results. scripts/run_early_verification.sh rebuilds missing meshes, solves the nine decks, runs the check scripts and compares the eight result files with the stored copies: all eight are byte-identical (V3a exits with status 1, as recorded), and the regenerated V4 temperature file is identical. On 2026-11-04 the V2 solver output files (.dat, .frd) found in simulation/ were empty although the V2 result file exists; the cause was not investigated; the rerun reproduced the stored V2 result exactly. Meshes and solver outputs are not committed (simulation/ and mesh/*.inp are ignored); mesh parameters are kept in mesh/*_info.txt. The later stages (V5, V6c, V6d, V7, V8a, V8b, baseline, study) are produced by their own run scripts; a clean rerun of those was not part of this audit, apart from the regression gate of Amendment 13.
