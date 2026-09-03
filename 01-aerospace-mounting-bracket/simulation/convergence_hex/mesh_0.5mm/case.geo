// Project 01 - Aerospace Mounting Bracket
// Structured hex convergence case: nx=120 ny=80 nz=8
// Scoped to verification plate ONLY - see chat log.

SetFactory("OpenCASCADE");
L = 60.0; B = 40.0; H = 4.0;
Box(1) = {0, 0, 0, L, B, H};

eps = 0.001;

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

Transfinite Curve {xline1(), xline2(), xline3(), xline4()} = 120;
Transfinite Curve {yline1(), yline2(), yline3(), yline4()} = 80;
Transfinite Curve {zline1(), zline2(), zline3(), zline4()} = 8;

all_surfaces() = Surface In BoundingBox{-eps,-eps,-eps, L+eps,B+eps,H+eps};
Transfinite Surface {all_surfaces()};
Recombine Surface {all_surfaces()};
Transfinite Volume {1};

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

Mesh.ElementOrder = 2;
Mesh.SecondOrderIncomplete = 1;
Mesh.SecondOrderLinear = 1;

Mesh 3;
