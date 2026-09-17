"""
Direct geometric test on the 462 bad elements: compute the signed
volume of each element's LINEAR (corner-node-only) tetrahedron.

If corner-only volume is positive and reasonably sized -> the base
linear element is fine; the problem is specifically introduced by
curved midside-node placement (order-2 curving artifact).

If corner-only volume is near-zero or negative -> the underlying
linear tet itself is degenerate/inverted -- a real meshing defect,
unrelated to curving.
"""

import gmsh
import numpy as np

MSH_FILE = "mesh/full_bracket_study/stage5_level1_v3.msh"
BAD_START, BAD_END = 170259, 170720

gmsh.initialize()
gmsh.open(MSH_FILE)

etypes, etags, enodes = gmsh.model.mesh.getElements(3)
tet_idx = [i for i, et in enumerate(etypes) if et == 11][0]
etags = np.array(etags[tet_idx])
enodes = np.array(enodes[tet_idx]).reshape(len(etags), 10)

mask = (etags >= BAD_START) & (etags <= BAD_END)
bad_etags = etags[mask]
bad_enodes = enodes[mask]

def signed_volume(p0, p1, p2, p3):
    return np.dot(np.cross(p1 - p0, p2 - p0), p3 - p0) / 6.0

vols = []
for tag, nn in zip(bad_etags, bad_enodes):
    corners = nn[:4]
    coords = np.array([gmsh.model.mesh.getNode(int(n))[0] for n in corners])
    v = signed_volume(coords[0], coords[1], coords[2], coords[3])
    vols.append(v)

vols = np.array(vols)
print(f"Linear (corner-only) tet volumes for {len(vols)} bad elements:")
print(f"  min={vols.min():.6f}  mean={vols.mean():.6f}  max={vols.max():.6f}")
print(f"  negative or near-zero (<1e-6): {(vols < 1e-6).sum()} / {len(vols)} "
      f"({100*(vols<1e-6).mean():.1f}%)")
print(f"  negative: {(vols < 0).sum()} / {len(vols)}")

# For reference: typical GOOD element volume in this mesh, for scale comparison
good_sample = etags[~mask][:2000]
good_vols = []
for tag in good_sample:
    idx = np.where(etags == tag)[0][0]
    corners = enodes[idx][:4]
    coords = np.array([gmsh.model.mesh.getNode(int(n))[0] for n in corners])
    good_vols.append(signed_volume(coords[0], coords[1], coords[2], coords[3]))
good_vols = np.array(good_vols)
print(f"\nFor scale -- 2000 sampled GOOD elements' linear volumes:")
print(f"  min={good_vols.min():.6f}  mean={good_vols.mean():.6f}  max={good_vols.max():.6f}")

gmsh.finalize()
