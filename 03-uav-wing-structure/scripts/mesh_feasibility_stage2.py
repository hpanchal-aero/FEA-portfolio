"""
Stage 2 mesh feasibility test -- generates a linear tet mesh at several
maxh sizes via Netgen and reports node/element counts and timing.
Mesh-only, no C3D10 promotion, no solve. Mirrors the Stage 1 pattern.
"""
import time
from netgen.occ import OCCGeometry
import netgen.meshing as ngm

STEP_FILE = "mesh/stage2_ribs/stage2_ribs.step"
MAXH_VALUES = [3.0, 2.0, 1.5]

# Empirical WSL host-crash ceilings established in Stage 1 (this machine,
# 11GB RAM, no .wslconfig cap): ~230,000 elements / ~1.18M equations for
# a direct-solver CalculiX job. These are mesh-only node/tet counts here
# (pre-C3D10-promotion); C3D10 promotion roughly doubles node count and
# each tet becomes one C3D10 element, so tet count here is directly
# comparable to the eventual element count.
ELEMENT_CEILING = 230_000

print("=" * 70)
print("STAGE 2 MESH FEASIBILITY (spar + ribs)")
print("=" * 70)

geo = OCCGeometry(STEP_FILE)

for maxh in MAXH_VALUES:
    t0 = time.time()
    ngmesh = geo.GenerateMesh(maxh=maxh)
    elapsed = time.time() - t0
    n_tets = ngmesh.ne
    n_nodes_linear = len(ngmesh.Points())
    pct_ceiling = 100 * n_tets / ELEMENT_CEILING
    print(f"maxh={maxh:.1f}mm: {n_nodes_linear} linear nodes, {n_tets} tets, "
          f"{elapsed:.1f}s ({pct_ceiling:.1f}% of {ELEMENT_CEILING} tet ceiling)")
