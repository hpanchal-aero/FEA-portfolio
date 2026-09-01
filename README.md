# FEA Portfolio — Computational Aerospace Structures & Mechanics

## Who I Am

I am an aerospace engineering student/researcher building a computational
engineering portfolio focused on structural analysis, finite element methods,
and (eventually) coupled thermal/fluid-structure problems relevant to
aerospace systems. This repository is aimed at demonstrating graduate-level
computational engineering capability, not software tutorials.

## What This Portfolio Is

A sequence of FEA and computational mechanics projects, each treated as a
self-contained engineering study following:

Problem → Theory → Modeling → Numerical Implementation → Verification →
Validation → Parametric Investigation → Engineering Interpretation →
Engineering Decision → Limitations

## Why I Am Building This

To demonstrate, with reproducible open-source tooling, that I can carry out
credible computational structural/aerospace engineering analysis suitable
for graduate research and technically serious engineering work.

## Toolchain

- **Geometry:** OpenSCAD
- **Meshing:** Gmsh
- **FEA Solver:** CalculiX
- **CFD:** OpenFOAM (later projects)
- **Post-processing:** ParaView, PyVista
- **Analysis/Automation:** Python (NumPy, SciPy, Pandas, Matplotlib)
- **Version Control:** Git / GitHub

## Engineering Philosophy

- A simple, well-verified problem is more valuable than a complex,
  unverified one.
- Mesh convergence and verification are mandatory, not optional.
- Verification (solving the math correctly) and validation (representing
  reality correctly) are always distinguished explicitly.
- Automation is introduced only after a single baseline case is confirmed
  correct.
- No result is reported without a stated basis for trusting it.

## Project Status

### In Progress
- **01 — Aerospace Mounting Bracket** — parametric ribbed L-bracket mass
  minimization under static equipment loading. OpenSCAD → Gmsh → CalculiX →
  Python → ParaView workflow.

### Planned
- **02 — 1U CubeSat Primary Structure** — structural architecture trade
  study under representative launch loading.
- **03 — UAV Wing Structural Analysis** — spar/rib/skin configuration study
  under aerodynamic loading.
- **04 — NACA 0012 CFD → Structural FEA** — one-way aerodynamic pressure
  coupling to structural response.
- **05 — Thermo-Structural Aerospace Panel** — thermal gradient → thermal
  stress workflow.
- **06 — UAV Wing Spar Mass Optimization** — constrained mass minimization
  with Pareto-style trade study.
- **07 — CubeSat Modal Analysis** — natural frequency / mode shape study.
- **08 — Thermal Protection / High-Temperature Structure** — aerothermal
  loading → thermal-structural response.

## Repository Structure

Each project directory follows:

NN-project-name/
├── README.md
├── geometry/
│ ├── cad/
│ └── exported/
├── mesh/
├── simulation/
├── scripts/
├── results/
└── figures/


## Reproducibility

Each project README documents the full environment, exact commands, and
verification/validation basis needed to reproduce results independently.
