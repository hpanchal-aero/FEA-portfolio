"""
Diagnostic: compare pre- and post-toe-fillet geometry directly,
including exact bounding box of every surface, to locate what
actually changed and where -- rather than relying on area-only
matching, which produced a false-positive identification.
"""

import gmsh
import os

PRE_BREP = "mesh/full_bracket_study/frame_plus_gusset_check.brep"
POST_BREP = "mesh/full_bracket_study/frame_gusset_toefillet_check.brep"

def dump_surfaces(label, brep_path):
    gmsh.initialize()
    gmsh.model.add(label)
    gmsh.model.occ.importShapes(brep_path)
    gmsh.model.occ.synchronize()

    vols = gmsh.model.getEntities(dim=3)
    print(f"\n=== {label} ===")
    print(f"File: {brep_path}")
    print(f"Volumes: {len(vols)}")
    for (dim, tag) in vols:
        mass = gmsh.model.occ.getMass(3, tag)
        print(f"  volume {tag}: mass = {mass:.6f} mm^3")

    surfs = gmsh.model.getEntities(dim=2)
    print(f"Surfaces: {len(surfs)}")
    data = []
    for (dim, tag) in surfs:
        area = gmsh.model.occ.getMass(2, tag)
        bb = gmsh.model.occ.getBoundingBox(dim, tag)
        data.append((tag, area, bb))
    data.sort(key=lambda d: d[1])
    for tag, area, bb in data:
        print(f"  surf {tag}: area={area:.4f}  bbox=[{bb[0]:.4f},{bb[1]:.4f},{bb[2]:.4f}] "
              f"-> [{bb[3]:.4f},{bb[4]:.4f},{bb[5]:.4f}]")

    gmsh.finalize()
    return data

pre_data = dump_surfaces("PRE-fillet (frame+gusset)", PRE_BREP)
post_data = dump_surfaces("POST-fillet", POST_BREP)

print("\n=== Focused check: any surface near the toe region (x~35-40, z~0-6)? ===")
for tag, area, bb in post_data:
    xmid = (bb[0] + bb[3]) / 2
    zmid = (bb[2] + bb[5]) / 2
    if 34.0 <= xmid <= 41.0 and -1.0 <= zmid <= 7.0:
        print(f"  surf {tag}: area={area:.4f}  bbox=[{bb[0]:.4f},{bb[1]:.4f},{bb[2]:.4f}] "
              f"-> [{bb[3]:.4f},{bb[4]:.4f},{bb[5]:.4f}]  (near-toe candidate)")
