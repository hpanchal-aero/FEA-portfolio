"""
Stage 3 (fillet) — surface identification check ONLY.

Rebuilds the exact same geometry as build_fillet_case.py (box union +
fillet), then queries every resulting surface's OCC type and bounding
box. Does NOT create a mesh size field, does NOT call mesh.generate(),
does NOT write any mesh file. Diagnostic only.
"""

import gmsh

TOL = 1e-6

LENGTH = 60.0
WIDTH  = 40.0
D_ROOT = 6.0
D_TIP  = 4.0
X_TRANS = 30.0
FILLET_R = 0.99

HALF_ROOT = D_ROOT / 2.0
HALF_TIP  = D_TIP / 2.0

gmsh.initialize()
gmsh.model.add("fillet_surface_check")
occ = gmsh.model.occ

root_box = occ.addBox(0.0, 0.0, -HALF_ROOT, X_TRANS, WIDTH, D_ROOT)
tip_box  = occ.addBox(X_TRANS, 0.0, -HALF_TIP, LENGTH - X_TRANS, WIDTH, D_TIP)
occ.synchronize()

out, out_map = occ.fuse([(3, root_box)], [(3, tip_box)])
occ.synchronize()
fused_vol = out[0][1]

all_curves = gmsh.model.getEntities(1)
top_edges = []
bottom_edges = []

for dim, tag in all_curves:
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
    x_at_trans = abs(xmin - X_TRANS) < TOL and abs(xmax - X_TRANS) < TOL
    spans_width = (ymax - ymin) > (WIDTH - 1.0)
    if not (x_at_trans and spans_width):
        continue
    z_mid = 0.5 * (zmin + zmax)
    z_flat = (zmax - zmin) < TOL
    if not z_flat:
        continue
    if abs(z_mid - HALF_TIP) < TOL:
        top_edges.append(tag)
    elif abs(z_mid + HALF_TIP) < TOL:
        bottom_edges.append(tag)

fillet_edges = top_edges + bottom_edges
radii = [FILLET_R] * len(fillet_edges)

filleted = occ.fillet([fused_vol], fillet_edges, radii, True)
occ.synchronize()
final_vol = filleted[0][1]

# ---- Surface identification, no meshing ----
print("=== ALL SURFACES: type + bounding box ===")
surfaces = gmsh.model.getEntities(2)
for dim, tag in sorted(surfaces, key=lambda s: s[1]):
    stype = gmsh.model.getType(dim, tag)
    bbox = gmsh.model.getBoundingBox(dim, tag)
    print(f"Surface {tag}: type={stype}  bbox={bbox}")

print("\n=== TARGETED CHECK: surfaces 3 and 11 ===")
for target_tag in [3, 11]:
    found = any(tag == target_tag for _, tag in surfaces)
    if not found:
        print(f"Surface {target_tag}: DOES NOT EXIST in this model")
        continue
    stype = gmsh.model.getType(2, target_tag)
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, target_tag)
    print(f"\nSurface {target_tag}:")
    print(f"  OCC type: {stype}")
    print(f"  Bounding box: [{xmin:.4f}, {xmax:.4f}] x [{ymin:.4f}, {ymax:.4f}] x [{zmin:.4f}, {zmax:.4f}]")

    # Expected fillet surface geometry:
    #  - spans x approx [X_TRANS - FILLET_R, X_TRANS + FILLET_R] = [29.01, 30.99]
    #  - spans full width y in [0, 40]
    #  - top fillet:    z approx [HALF_TIP, HALF_TIP + FILLET_R]     = [2.00, 2.99]
    #  - bottom fillet: z approx [-(HALF_TIP+FILLET_R), -HALF_TIP]  = [-2.99, -2.00]
    x_expected = abs(xmin - (X_TRANS - FILLET_R)) < 0.05 and abs(xmax - (X_TRANS + FILLET_R)) < 0.05
    y_expected = abs(ymin - 0.0) < 0.05 and abs(ymax - WIDTH) < 0.05
    z_top_expected = abs(zmin - HALF_TIP) < 0.05 and abs(zmax - (HALF_TIP + FILLET_R)) < 0.05
    z_bottom_expected = abs(zmax - (-HALF_TIP)) < 0.05 and abs(zmin - (-(HALF_TIP + FILLET_R))) < 0.05

    is_top_fillet = stype == "Cylinder" and x_expected and y_expected and z_top_expected
    is_bottom_fillet = stype == "Cylinder" and x_expected and y_expected and z_bottom_expected

    if is_top_fillet:
        print("  CONFIRMED: matches expected TOP fillet surface geometry")
    elif is_bottom_fillet:
        print("  CONFIRMED: matches expected BOTTOM fillet surface geometry")
    else:
        print("  NOT CONFIRMED as a fillet surface — geometry does not match expected fillet location/extent")

gmsh.finalize()
