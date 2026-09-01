// Project 01 - Aerospace Mounting Bracket
// Stage: Baseline verification mesh (coarse - convergence study level 1)
//
// Converts OpenSCAD STL export into a Gmsh volumetric tetrahedral mesh
// via surface classification (STL has no native CAD topology).
//
// NOTE: ClassifySurfaces splits each side face of the box into two
// triangulated patches (STL export artifact from OpenSCAD, confirmed
// via bounding-box inspection - see scripts/inspect_mesh_surfaces.py).
// This does not affect the volume mesh; Physical Surface groups below
// recombine each split pair into the correct complete physical face
// for boundary-condition tagging.

Merge "../geometry/exported/verification_plate.stl";

angle = 40;
ClassifySurfaces{angle * Pi/180, 1, 1, angle * Pi/180};
CreateGeometry;

Surface Loop(1) = Surface{:};
Volume(1) = {1};

// ---- Physical groups for boundary conditions ----
// Confirmed via inspect_mesh_surfaces.py bounding boxes:
Physical Surface("root", 101) = {10, 11};   // x=0 face - fixed BC
Physical Surface("tip", 102)  = {6, 7};     // x=60 face - load application
Physical Surface("bottom", 103) = {3};      // z=0 face
Physical Surface("top", 104)    = {2};      // z=4 face
Physical Surface("side_y0", 105)  = {4, 5}; // y=0 face
Physical Surface("side_y40", 106) = {8, 9}; // y=40 face
Physical Volume("plate", 201) = {1};

// ---- Mesh settings ----
Mesh.CharacteristicLengthMax = 1.3;
Mesh.CharacteristicLengthMin = 0.5;
Mesh.ElementOrder = 2;
Mesh.Algorithm3D = 1;

Mesh 3;
