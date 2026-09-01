# Project 01 — Ribbed L-Bracket: Mass-Minimized Aerospace Mounting Bracket

**Status:** In progress — specification and skeleton stage.

## Research Question

How can the mass of a ribbed L-bracket mounting interface be minimized
while satisfying stress, displacement, and factor-of-safety constraints
under a representative static equipment load?

## Engineering Motivation

Mounting brackets are ubiquitous in aerospace structures (avionics,
instruments, subsystem hardware) and are a common target for structural
mass optimization, since bracket mass is rarely mission-critical
functionally but directly impacts overall vehicle mass budget.

## Geometry (Baseline)

Angle bracket with two mounting flanges (perpendicular planes) joined by a
central triangular stiffening gusset:
- Flange A — fixed to primary structure (bolted, encastre boundary)
- Flange B — bolted to mounted equipment (load application surface)

## Material

Al 7075-T6
- E = 71.7 GPa, ν = 0.33, ρ = 2810 kg/m³
- σ_yield ≈ 503 MPa (source to be cited explicitly in material definition file)

## Loading

- Equipment mass: 2 kg
- Design acceleration: 12g quasi-static — **preliminary assumed value, not
  a sourced launch/vibration requirement**
- Resultant static load: F = 2 kg × 9.81 m/s² × 12 = 235.44 N, applied at
  Flange B
- Flange A: fixed (encastre) at bolt-hole locations

## Methodology (planned)

OpenSCAD (geometry) → Gmsh (mesh) → CalculiX (linear static FEA) → Python
(post-processing, mesh convergence, parametric study) → ParaView (visualization)

## Verification Strategy

- Mesh convergence study (h-refinement) monitoring max von Mises stress and
  Flange B displacement
- Analytical cross-check against simplified cantilever plate/beam bending
  solution

## Validation Strategy

No physical test article or published dataset exists for this specific
custom geometry. Validation in the strict sense (does the model represent
reality) is **not performed** for this project — this is documented as an
explicit limitation. Only verification (is the FE model solving the
mathematical model correctly) is claimed.

## Directory Structure

01-aerospace-mounting-bracket/
├── README.md
├── geometry/
│ ├── cad/ — OpenSCAD source files
│ └── exported/ — STEP/STL exports for meshing
├── mesh/ — Gmsh .geo and .msh files
├── simulation/ — CalculiX .inp files
├── scripts/ — Python automation/post-processing
├── results/ — CalculiX outputs, CSV summaries
└── figures/ — plots, ParaView renders


## Reproducibility

Environment and exact commands will be documented here as each stage is
completed.

## Status Log

- [x] Specification defined and approved
- [ ] Baseline geometry (OpenSCAD)
- [ ] Baseline mesh (Gmsh)
- [ ] Baseline CalculiX model
- [ ] Baseline solve + sanity check
- [ ] Mesh convergence study
- [ ] Analytical verification
- [ ] Parametric study
- [ ] Final documentation and figures
