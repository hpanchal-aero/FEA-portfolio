// Project 01 - Aerospace Mounting Bracket
// Stage 2: Plate with hole - stress concentration verification
// Single test mesh (far_size=0.5, hole_size=0.15) before convergence study.
//
// Geometry: exact OCC box with a through-thickness circular hole cut
// via boolean difference - no STL faceting, exact circular boundary
// (critical for a stress-concentration study).
//
// Hole: d=5mm, centered at (x=30, y=20), through full thickness.
// Located at mid-span to avoid the root Saint-Venant boundary layer
// characterized in the verification-plate stage.
//
// Mesh: unstructured tetrahedra (C3D10) with local refinement at the
// hole boundary via a Distance+Threshold field - structured hex
// cannot wrap a circular cutout.

SetFactory("OpenCASCADE");
L = 60.0; B = 40.0; H = 4.0;
hole_x = 30.0; hole_y = 20.0; hole_d = 5.0; hole_r = hole_d/2.0;

Box(1) = {0, 0, 0, L, B, H};
Cylinder(2) = {hole_x, hole_y, -1.0, 0, 0, H+2.0, hole_r};
BooleanDifference(3) = { Volume{1}; Delete; }{ Volume{2}; Delete; };

eps = 1e-3;

root_surf()      = Surface In BoundingBox{-eps,-eps,-eps, eps,B+eps,H+eps};
tip_surf()       = Surface In BoundingBox{L-eps,-eps,-eps, L+eps,B+eps,H+eps};
bottom_surf()    = Surface In BoundingBox{-eps,-eps,-eps, L+eps,B+eps,eps};
top_surf()       = Surface In BoundingBox{-eps,-eps,H-eps, L+eps,B+eps,H+eps};
side_y0_surf()   = Surface In BoundingBox{-eps,-eps,-eps, L+eps,eps,H+eps};
side_y40_surf()  = Surface In BoundingBox{-eps,B-eps,-eps, L+eps,B+eps,H+eps};
hole_surf()      = Surface In BoundingBox{hole_x-hole_r-eps, hole_y-hole_r-eps, -eps,
                                            hole_x+hole_r+eps, hole_y+hole_r+eps, H+eps};

Physical Surface("root", 101)     = root_surf();
Physical Surface("tip", 102)      = tip_surf();
Physical Surface("bottom", 103)   = bottom_surf();
Physical Surface("top", 104)      = top_surf();
Physical Surface("side_y0", 105)  = side_y0_surf();
Physical Surface("side_y40", 106) = side_y40_surf();
Physical Surface("hole", 107)     = hole_surf();
Physical Volume("plate", 201)     = {3};

// ---- Local mesh refinement at the hole boundary ----
far_size = 3.0;
hole_size = 0.25;
transition_width = 1.5;   // mm - WIDTH of the transition band, beyond
                            // the hole edge (not absolute distance
                            // from hole center - see chat log, this
                            // was the source of a field-definition
                            // bug that caused a 2.5M-element mesh)

Field[1] = Distance;
Field[1].SurfacesList = {hole_surf()};

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = hole_size;
Field[2].SizeMax = far_size;
Field[2].DistMin = hole_r;
Field[2].DistMax = hole_r + transition_width;

Background Field = 2;

Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.CharacteristicLengthMax = far_size;

Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 1;
Mesh.Algorithm3D = 10;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;

Mesh 3;
