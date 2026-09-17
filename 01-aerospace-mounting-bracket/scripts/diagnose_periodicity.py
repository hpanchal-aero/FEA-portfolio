"""
Stage 5 -- investigate the suspected periodicity mechanism behind the
462 persistently-curved order-2 elements. The extrude() operations
used to build the gusset and base frame in OCC can create implicit
periodic (master/slave) entity links between source and target faces
of an extrusion, even without the user requesting periodic meshing.
If Gmsh's order-2 node placement snaps a "slave" midside node to
match its "master" counterpart's position rather than computing it
locally from its own edge, that would explain deterministic,
setting-independent offsets exactly like what we observed.

This script queries gmsh.model.getPeriodic() for every entity in the
geometry to check whether any such links exist.
"""

import gmsh

INPUT_BREP = "mesh/full_bracket_study/full_bracket_geometry_final.brep"

gmsh.initialize()
gmsh.model.add("periodicity_check")
gmsh.model.occ.importShapes(INPUT_BREP)
gmsh.model.occ.synchronize()

print("Checking for periodic entity links at each dimension...\n")

found_any = False
for dim in [0, 1, 2, 3]:
    entities = gmsh.model.getEntities(dim)
    for (d, tag) in entities:
        try:
            master_tag, node_map, affine = gmsh.model.getPeriodic(d, [tag])
        except Exception as e:
            continue
        # getPeriodic returns the tag itself if no periodicity; check if
        # it actually returned something different / meaningful
        if master_tag and master_tag != [tag]:
            found_any = True
            print(f"dim={d} tag={tag}: PERIODIC, master={master_tag}")

if not found_any:
    print("No periodic entity links found via getPeriodic() at any dimension.")
    print("(This would rule out the periodicity hypothesis entirely.)")

gmsh.finalize()
