# Project 02 — 1U CubeSat Primary Structure — Development Log

## Specification (locked)

**Research question**: How does CubeSat structural configuration
affect mass, stiffness, and stress margin under representative
launch loading?

**Material**: Al 6061-T6 (E=68,900 MPa, ν=0.33, ρ=2700 kg/m³,
yield≈276 MPa — standard published values; exact source to be
confirmed and cited when material properties are locked into the
first .inp file).

**Loading**: 12g quasi-static axial launch acceleration. Explicitly a
**preliminary design assumption, not a sourced launch-vehicle
requirement** — same status as Project 01's 12g load assumption,
carried forward for methodological consistency across the portfolio.

**Structural configurations to compare**:
1. Configuration A — solid-panel frame (flat panels, corner rails, no
   lightening features)
2. Configuration B — ribbed/lightened frame (same envelope, panels
   with lightening pockets or rib pattern)

**Constraints**:
- Max von Mises stress ≤ yield/FoS, FoS=1.5 → allowable ≈184 MPa
- Mass checked against the 1U CDS budget (1.33 kg) as a design
  target, not necessarily a hard pass/fail constraint — exact current
  CDS revision figure to be verified before being treated as binding
- No displacement limit set — CubeSat panel structures are typically
  stiffness/frequency-governed rather than displacement-governed;
  open to revision if the FEA results suggest otherwise

**Explicitly out of scope**: modal/frequency response, random
vibration loading — deferred to Stage 7 per the original portfolio
sequencing (Project 07, CubeSat modal analysis).

**Verification strategy**: no single closed-form reference exists for
the full assembled structure, unlike Project 01's progressive
bracket stages. Plan: verify individual features analytically where
a clean reference exists (single panel as plate bending, corner rail
as simple column/beam) before combining into the full assembly. The
full assembly itself will be treated as an internal consistency
study (equilibrium checks, mesh convergence) rather than claimed as
closed-form-verified — the same honest framing Project 01 used for
its own Stage 5 full-bracket integration.

**Toolchain**: OpenSCAD/Gmsh/Netgen (geometry/meshing), CalculiX
2.23 (solver, ccx223 conda environment), Python (NumPy, PyVista,
Matplotlib) — same locked toolchain as Project 01's Stage 5 onward.

**Directory structure**:
02-cubesat-structure/
├── README.md
├── development-log.md
├── geometry/
├── mesh/
├── simulation/
├── scripts/
├── results/
└── figures/

## Status Log

- [x] Specification defined and approved
- [ ] Baseline geometry (Configuration A)
- [ ] Baseline mesh
- [ ] Baseline solve + sanity check
- [ ] Mesh convergence study
- [ ] Feature-level analytical verification (panel, corner rail)
- [ ] Configuration B (ribbed/lightened) geometry
- [ ] Comparative analysis (A vs. B: mass, stress, stiffness)
- [ ] Documentation and figures
- [ ] Git commit and push
