// Project 01 - Aerospace Mounting Bracket
// STANDALONE TEST: structured transfinite hexahedral mesh
// SCOPED TO THE SIMPLE VERIFICATION PLATE ONLY - see chat log.
// This is NOT adopted as the general bracket meshing strategy; the
// final ribbed/filleted bracket will require a separate mesh
// methodology decision.
//
// Rationale: unstructured tetrahedral meshing (Delaunay and HXT,
// both with SecondOrderLinear straight-edge midsides, Netgen
// optimization) repeatedly produced isolated C3D10 elements with
// nonpositive Jacobian determinants that could not be resolved
// through mesh-quality/size controls (see chat log investigation).
// A structured transfinite hex mesh is generated analytically from
// a regular grid, not from Delaunay/HXT tetrahedralization, so
// element regularity is guaranteed by construction rather than
// discovered after the fact.
//
// Element type: C3D20 (20-node quadratic serendipity hexahedron)

SetFactory("OpenCASCADE");
L = 60.0; B = 40.0; H = 4.0;
Box(1) = {0, 0, 0, L, B, H};

eps = 1e-3;

// Divisions in POINTS (n_elements = n_points - 1), giving uniform
// 0.8mm elements in all three directions at this refinement level.
nx = 76; // 75 elements along L
ny = 51; // 50 elements along B
nz = 6;  // 5 elements along H

// --- Identify the 12 box edges by direction via bounding-box query
//     (same deterministic approach used for surfaces earlier -
//     never assume Gmsh's internal numbering) ---
xline1() = Curve In BoundingBox{-eps,-eps,-eps, L+eps,eps,eps};
xline2() = Curve In BoundingBox{-eps,-eps,H-eps, L+eps,eps,H+eps};
xline3() = Curve In BoundingBox{-eps,B-eps,-eps, L+eps,B+eps,eps};
xline4() = Curve In BoundingBox{-eps,B-eps,H-eps, L+eps,B+eps,H+eps};

yline1() = Curve In BoundingBox{-eps,-eps,-eps, eps,B+eps,eps};
yline2() = Curve In BoundingBox{-eps,-eps,H-eps, eps,B+eps,H+eps};
yline3() = Curve In BoundingBox{L-eps,-eps,-eps, L+eps,B+eps,eps};
yline4() = Curve In BoundingBox{L-eps,-eps,H-eps, L+eps,B+eps,H+eps};

zline1() = Curve In BoundingBox{-eps,-eps,-eps, eps,eps,H+eps};
zline2() = Curve In BoundingBox{-eps,B-eps,-eps, eps,B+eps,H+eps};
zline3() = Curve In BoundingBox{L-eps,-eps,-eps, L+eps,eps,H+eps};
zline4() = Curve In BoundingBox{L-eps,B-eps,-eps, L+eps,B+eps,H+eps};

Transfinite Curve {xline1(), xline2(), xline3(), xline4()} = nx;
Transfinite Curve {yline1(), yline2(), yline3(), yline4()} = ny;
Transfinite Curve {zline1(), zline2(), zline3(), zline4()} = nz;

// --- Surfaces: transfinite + recombine into quads ---
all_surfaces() = Surface In BoundingBox{-eps,-eps,-eps, L+eps,B+eps,H+eps};
Transfinite Surface {all_surfaces()};
Recombine Surface {all_surfaces()};

// --- Volume: transfinite structured meshing produces hexahedra ---
Transfinite Volume {1};

// --- Face identification for boundary conditions ---
root_surf()      = Surface In BoundingBox{-eps,-eps,-eps, eps,B+eps,H+eps};
tip_surf()       = Surface In BoundingBox{L-eps,-eps,-eps, L+eps,B+eps,H+eps};
bottom_surf()    = Surface In BoundingBox{-eps,-eps,-eps, L+eps,B+eps,eps};
top_surf()       = Surface In BoundingBox{-eps,-eps,H-eps, L+eps,B+eps,H+eps};
side_y0_surf()   = Surface In BoundingBox{-eps,-eps,-eps, L+eps,eps,H+eps};
side_y40_surf()  = Surface In BoundingBox{-eps,B-eps,-eps, L+eps,B+eps,H+eps};

Physical Surface("root", 101)     = root_surf();
Physical Surface("tip", 102)      = tip_surf();
Physical Surface("bottom", 103)   = bottom_surf();
Physical Surface("top", 104)      = top_surf();
Physical Surface("side_y0", 105)  = side_y0_surf();
Physical Surface("side_y40", 106) = side_y40_surf();
Physical Volume("plate", 201)     = {1};

// --- Element type: C3D20 (20-node serendipity hex, no interior nodes) ---
Mesh.ElementOrder = 2;
Mesh.SecondOrderIncomplete = 1;  // 20-node, not 27-node hex
Mesh.SecondOrderLinear = 1;      // straight edges - geometry is exactly flat

Mesh 3;
