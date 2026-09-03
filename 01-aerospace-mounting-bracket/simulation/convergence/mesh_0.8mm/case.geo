// Project 01 - Aerospace Mounting Bracket
// Convergence study case: element size = 0.8 mm
// Geometry source: Gmsh OpenCASCADE exact box (NOT the OpenSCAD STL -
// see chat log for rationale: STL faceting caused nonpositive-Jacobian
// tets near box corners at fine mesh sizes).
// Reference design record: geometry/cad/01_verification_plate.scad
// (identical dimensions: L=60, B=40, H=4 mm)

SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 60.0, 40.0, 4.0};

eps = 0.001;
root_surf()      = Surface In BoundingBox{-eps,-eps,-eps, eps,40.0+eps,4.0+eps};
tip_surf()       = Surface In BoundingBox{60.0-eps,-eps,-eps, 60.0+eps,40.0+eps,4.0+eps};
bottom_surf()    = Surface In BoundingBox{-eps,-eps,-eps, 60.0+eps,40.0+eps,eps};
top_surf()       = Surface In BoundingBox{-eps,-eps,4.0-eps, 60.0+eps,40.0+eps,4.0+eps};
side_y0_surf()   = Surface In BoundingBox{-eps,-eps,-eps, 60.0+eps,eps,4.0+eps};
side_y40_surf()  = Surface In BoundingBox{-eps,40.0-eps,-eps, 60.0+eps,40.0+eps,4.0+eps};

Physical Surface("root", 101)     = root_surf();
Physical Surface("tip", 102)      = tip_surf();
Physical Surface("bottom", 103)   = bottom_surf();
Physical Surface("top", 104)      = top_surf();
Physical Surface("side_y0", 105)  = side_y0_surf();
Physical Surface("side_y40", 106) = side_y40_surf();
Physical Volume("plate", 201)     = {1};

Mesh.CharacteristicLengthMax = 0.8;
Mesh.CharacteristicLengthMin = 0.56;
Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 1;
Mesh.Algorithm3D = 10;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;

Mesh 3;
