"""
Project 02 -- Configuration A (solid-panel frame) baseline geometry.

Locked spec: 100x100x100mm envelope, four 8.5x8.5x100mm corner rails
(outer faces flush with envelope), four 2mm-thick side panels spanning
the gaps between rails (outer faces flush with envelope), fused into
a single solid. Top/bottom panels and bolted-joint fastener detail
are explicitly OMITTED -- documented scope simplifications, not
oversights.

Geometry construction and validation ONLY. Does not mesh or solve.

REVISION NOTE: the original version left the z=50 section face from the
sectional-area check in the model, so the saved BRep was a compound of the
solid plus one stray 953 mm^2 face. This version removes the section
residue, asserts the model is exactly 1 volume / 18 surfaces before saving,
and also writes a STEP file.
"""

import gmsh
import sys
import os

TOL = 1e-6

ENV = 50.0          # half-envelope, mm (full envelope 100x100)
RAIL = 8.5           # rail cross-section side, mm
HEIGHT = 100.0        # z-height, mm
PANEL_T = 2.0         # panel thickness, mm
RHO = 2700e-9         # Al 6061-T6 density, kg/mm^3

RAIL_INNER = ENV - RAIL   # 41.5, inner face of each rail

EXPECTED_RAIL_VOL_EACH = RAIL * RAIL * HEIGHT              # 7225
EXPECTED_PANEL_LEN = 2 * RAIL_INNER                          # 83.0
EXPECTED_PANEL_VOL_EACH = EXPECTED_PANEL_LEN * PANEL_T * HEIGHT  # 16600
EXPECTED_TOTAL_VOL = 4 * EXPECTED_RAIL_VOL_EACH + 4 * EXPECTED_PANEL_VOL_EACH  # 95300
EXPECTED_MASS = EXPECTED_TOTAL_VOL * RHO
EXPECTED_XSEC_AREA = 4 * (RAIL * RAIL) + 4 * (EXPECTED_PANEL_LEN * PANEL_T)  # 953
EXPECTED_N_VOLUMES = 1
EXPECTED_N_SURFACES = 18

OUT_DIR = "mesh/config_a_study"
os.makedirs(OUT_DIR, exist_ok=True)

gmsh.initialize()
gmsh.model.add("config_a_base_frame")
occ = gmsh.model.occ

print("=" * 70)
print("STEP 1: RAILS -- pre-fuse cross-section check")
print("=" * 70)

rail_positions = [
    (-ENV, -ENV),
    (RAIL_INNER, -ENV),
    (-ENV, RAIL_INNER),
    (RAIL_INNER, RAIL_INNER),
]
rail_tags = []
for i, (x0, y0) in enumerate(rail_positions, start=1):
    tag = occ.addBox(x0, y0, 0.0, RAIL, RAIL, HEIGHT)
    occ.synchronize()
    bbox = gmsh.model.getBoundingBox(3, tag)
    xext = bbox[3] - bbox[0]
    yext = bbox[4] - bbox[1]
    zext = bbox[5] - bbox[2]
    ok = (abs(xext - RAIL) < TOL and abs(yext - RAIL) < TOL and abs(zext - HEIGHT) < TOL)
    print(f"  Rail {i}: bbox extent x={xext:.4f} y={yext:.4f} z={zext:.4f} "
          f"[{'OK' if ok else 'MISMATCH -- STOP'}]")
    if not ok:
        print("FATAL: rail pre-fuse dimension check failed.")
        gmsh.finalize()
        sys.exit(1)
    rail_tags.append(tag)

print("\n" + "=" * 70)
print("STEP 2: PANELS -- pre-fuse dimension check")
print("=" * 70)

# Panels span the gap between rails on each of the 4 sides, full height,
# PANEL_T thick measured inward from the envelope face.
panel_defs = [
    ("y=-ENV", -RAIL_INNER, -ENV, EXPECTED_PANEL_LEN, PANEL_T),           # x:[-41.5,41.5], y:[-50,-48]
    ("y=+ENV", -RAIL_INNER, ENV - PANEL_T, EXPECTED_PANEL_LEN, PANEL_T),  # y:[48,50]
    ("x=-ENV", -ENV, -RAIL_INNER, PANEL_T, EXPECTED_PANEL_LEN),           # x:[-50,-48], y:[-41.5,41.5]
    ("x=+ENV", ENV - PANEL_T, -RAIL_INNER, PANEL_T, EXPECTED_PANEL_LEN),  # x:[48,50]
]
panel_tags = []
for label, x0, y0, dx, dy in panel_defs:
    tag = occ.addBox(x0, y0, 0.0, dx, dy, HEIGHT)
    occ.synchronize()
    bbox = gmsh.model.getBoundingBox(3, tag)
    xext = bbox[3] - bbox[0]
    yext = bbox[4] - bbox[1]
    zext = bbox[5] - bbox[2]
    ok = (abs(xext - dx) < TOL and abs(yext - dy) < TOL and abs(zext - HEIGHT) < TOL)
    print(f"  Panel [{label}]: bbox extent x={xext:.4f} y={yext:.4f} z={zext:.4f} "
          f"[{'OK' if ok else 'MISMATCH -- STOP'}]")
    if not ok:
        print("FATAL: panel pre-fuse dimension check failed.")
        gmsh.finalize()
        sys.exit(1)
    panel_tags.append(tag)

print("\n" + "=" * 70)
print("STEP 3: VOLUME OVERLAP CHECK (rails vs panels, before fusing)")
print("=" * 70)
# Copy rail/panel solids so the intersection test doesn't consume the
# originals we still need for the real fuse below.
rail_copies = [occ.copy([(3, t)])[0][1] for t in rail_tags]
panel_copies = [occ.copy([(3, t)])[0][1] for t in panel_tags]
try:
    intersect_result, _ = occ.intersect(
        [(3, t) for t in rail_copies], [(3, t) for t in panel_copies],
        removeObject=True, removeTool=True
    )
    occ.synchronize()
    overlap_vol = sum(occ.getMass(3, t) for d, t in intersect_result) if intersect_result else 0.0
except Exception:
    overlap_vol = 0.0
print(f"  Rail/panel intersection volume: {overlap_vol:.6f} mm^3 "
      f"[{'OK, negligible' if overlap_vol < 1e-6 else 'WARNING -- nonzero overlap, expected volume calc invalid'}]")

print("\n" + "=" * 70)
print("STEP 4: FUSE all rails + panels into a single solid")
print("=" * 70)
all_tags = rail_tags + panel_tags
fused, _ = occ.fuse([(3, all_tags[0])], [(3, t) for t in all_tags[1:]])
occ.synchronize()
if len(fused) != 1:
    print(f"FATAL: fuse produced {len(fused)} volumes, expected 1.")
    gmsh.finalize()
    sys.exit(1)
vol_tag = fused[0][1]
print(f"  Fuse OK: 1 volume, tag={vol_tag}")

print("\n" + "=" * 70)
print("STEP 5: POST-FUSE VALIDATION")
print("=" * 70)

# 5a. Mass check
actual_mass_vol = occ.getMass(3, vol_tag)  # this returns volume (mm^3) for a 3D entity
pct_diff = 100 * (actual_mass_vol - EXPECTED_TOTAL_VOL) / EXPECTED_TOTAL_VOL
print(f"  Volume: actual={actual_mass_vol:.4f} mm^3, expected={EXPECTED_TOTAL_VOL:.4f} mm^3, "
      f"%diff={pct_diff:+.4f}%")
vol_ok = abs(pct_diff) < 0.1
print(f"  Volume check: [{'OK' if vol_ok else 'MISMATCH -- STOP'}]")

# 5b. Bounding box check
bbox = gmsh.model.getBoundingBox(3, vol_tag)
print(f"  Bbox: x=[{bbox[0]:.4f},{bbox[3]:.4f}] y=[{bbox[1]:.4f},{bbox[4]:.4f}] "
      f"z=[{bbox[2]:.4f},{bbox[5]:.4f}]  (expect x,y=[-50,50], z=[0,100])")
bbox_ok = (abs(bbox[0] + ENV) < TOL and abs(bbox[3] - ENV) < TOL and
           abs(bbox[1] + ENV) < TOL and abs(bbox[4] - ENV) < TOL and
           abs(bbox[2] - 0.0) < TOL and abs(bbox[5] - HEIGHT) < TOL)
print(f"  Bbox check: [{'OK' if bbox_ok else 'MISMATCH -- STOP'}]")

# 5c. Sliver/duplicate surface check
surfaces = gmsh.model.getEntities(2)
areas = [(occ.getMass(2, t), t) for d, t in surfaces]
areas.sort()
slivers = [(a, t) for a, t in areas if a < 1.0]
print(f"  Surfaces: {len(surfaces)}, slivers (<1.0mm^2): {len(slivers)}")
if slivers:
    for a, t in slivers[:10]:
        print(f"    sliver surf {t}: area={a:.6f}")
sliver_ok = len(slivers) == 0

# 5d. Sectional check (per your adjustment -- replaces post-fuse rail
# cross-section indexing, which boolean ops can invalidate by splitting
# surfaces). Cut a plane at z=50 and measure the resulting cross-section
# area against the analytically expected 953 mm^2.
print(f"\n  Sectional check at z=50 (expected area {EXPECTED_XSEC_AREA:.4f} mm^2):")
section_result = []
try:
    plane_tag = occ.addDisk(0, 0, 50.0, 200, 200)  # oversized disk, will be clipped by intersection
    occ.synchronize()
    solid_copy = occ.copy([(3, vol_tag)])[0][1]
    section_result, _ = occ.intersect([(2, plane_tag)], [(3, solid_copy)],
                                       removeObject=True, removeTool=True)
    occ.synchronize()
    section_area = sum(occ.getMass(2, t) for d, t in section_result if d == 2)
    xsec_pct_diff = 100 * (section_area - EXPECTED_XSEC_AREA) / EXPECTED_XSEC_AREA
    print(f"    Measured section area: {section_area:.4f} mm^2, %diff={xsec_pct_diff:+.4f}%")
    xsec_ok = abs(xsec_pct_diff) < 0.5
except Exception as e:
    print(f"    Sectional check FAILED to compute: {e}")
    xsec_ok = False
print(f"  Sectional check: [{'OK' if xsec_ok else 'MISMATCH -- STOP'}]")

# 5e. CLEANUP: the section check leaves its 2D result face in the model.
# Remove it so it is not written into the geometry file.
if section_result:
    occ.remove(section_result, recursive=True)
    occ.synchronize()

# 5f. Model-cleanliness check: exactly 1 volume and 18 surfaces must remain.
n_vols = len(gmsh.model.getEntities(3))
n_surfs = len(gmsh.model.getEntities(2))
print(f"\n  Model after cleanup: {n_vols} volume(s), {n_surfs} surface(s) "
      f"(expect {EXPECTED_N_VOLUMES} and {EXPECTED_N_SURFACES})")
clean_ok = (n_vols == EXPECTED_N_VOLUMES and n_surfs == EXPECTED_N_SURFACES)
print(f"  Cleanliness check: [{'OK' if clean_ok else 'MISMATCH -- STOP'}]")

print("\n" + "=" * 70)
all_ok = vol_ok and bbox_ok and sliver_ok and xsec_ok and clean_ok
print(f"OVERALL VALIDATION: {'PASS' if all_ok else 'FAIL'}")
print("=" * 70)

if all_ok:
    brep_path = os.path.join(OUT_DIR, "config_a_base_frame.brep")
    step_path = os.path.join(OUT_DIR, "config_a_clean.step")
    gmsh.write(brep_path)
    gmsh.write(step_path)
    print(f"\nSaved: {brep_path}")
    print(f"Saved: {step_path}")
    print(f"Mass at Al 6061-T6 density: {actual_mass_vol * RHO * 1000:.2f} g "
          f"(expected ~{EXPECTED_MASS*1000:.2f} g)")
else:
    print("\nGeometry NOT saved -- validation failed. Fix before proceeding.")

gmsh.finalize()
sys.exit(0 if all_ok else 1)
