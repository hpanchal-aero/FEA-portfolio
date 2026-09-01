// Project 01 - Aerospace Mounting Bracket
// Stage: Baseline Verification Geometry
//
// Purpose: plain flat cantilever plate for Euler-Bernoulli analytical
// verification (stress + tip deflection) and mesh convergence study.
// This geometry deliberately excludes holes, fillet, and gusset -
// those are introduced in later geometry files after this baseline
// is verified against analytical theory.
//
// Idealization:
//   - Root (x=0 plane): fixed (encastre) boundary condition, applied
//     at the solver stage, not represented here.
//   - Free tip (x=L plane): statically equivalent transverse point
//     force of 235.44 N applied at the free-edge centroid
//     (x=L, y=0, z=h/2) during FEA - also applied at solver stage.
//
// Axes:
//   x - cantilever / load-path direction (length L)
//   y - flange width direction (width b)
//   z - thickness direction (height h)

// ---- Verification geometry parameters (mm) ----
L = 60;   // cantilever length, load-path direction
b = 40;   // flange width
h = 4;    // thickness

// ---- Geometry ----
cube([L, b, h]);
