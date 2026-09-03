// Project 01 - Aerospace Mounting Bracket
// Convergence study case: element size = 1.3 mm

Merge "/home/harsh/FEA-portfolio/01-aerospace-mounting-bracket/geometry/exported/verification_plate.stl";

angle = 40;
ClassifySurfaces{angle * Pi/180, 1, 1, angle * Pi/180};
CreateGeometry;

Surface Loop(1) = Surface{:};
Volume(1) = {1};

Physical Surface("root", 101) = {10, 11};
Physical Surface("tip", 102)  = {6, 7};
Physical Surface("bottom", 103) = {3};
Physical Surface("top", 104)    = {2};
Physical Surface("side_y0", 105)  = {4, 5};
Physical Surface("side_y40", 106) = {8, 9};
Physical Volume("plate", 201) = {1};

Mesh.CharacteristicLengthMax = 1.3;
Mesh.CharacteristicLengthMin = 0.5;
Mesh.ElementOrder = 2;
Mesh.Algorithm3D = 1;

Mesh 3;
