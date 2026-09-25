# FEA Portfolio — Computational Aerospace Structures & Mechanics

A structured computational mechanics portfolio focused on aerospace structures, finite element analysis, and the progression toward coupled thermal and fluid–structural problems.

The projects move from fundamental structural analysis toward aerospace-specific loading, optimization, modal analysis, and eventually coupled aerodynamic and aerothermal problems.

The objective is to develop the ability to **formulate, model, verify, validate, and interpret computational structural problems** rather than simply demonstrate FEA software.

---

## What This Portfolio Builds

The projects are designed to develop experience across:

* Static structural analysis
* Aerospace structural modelling
* Load-path and structural configuration studies
* Mesh convergence and numerical verification
* Modal analysis
* Structural optimization
* Aerodynamic pressure → structural response
* Thermal → structural coupling
* High-temperature aerospace structures
* Fluid–structure and thermo-structural analysis

The sequence is intentional: establish reliable structural mechanics fundamentals before introducing increasingly coupled multiphysics problems.

---

## Project Roadmap

| #  | Project                                                                               | Status         | Focus                                                             |
| -- | ------------------------------------------------------------------------------------- | -------------- | ----------------------------------------------------------------- |
| 01 | [**Aerospace Mounting Bracket**](01-aerospace-mounting-bracket/)                      | 🟢 Completed | Parametric ribbed L-bracket, static loading, mass minimization    |
| 02 | [**1U CubeSat Primary Structure**](02-cubesat-structure/)                             | 🟢 Completed | Structural architecture trade under representative launch loading |
| 03 | [**UAV Wing Structural Analysis**](03-uav-wing-structural-analysis/)                  | 🟡 In Progress | Spar/rib/skin configuration under aerodynamic loading             |
| 04 | [**NACA 0012 CFD → Structural FEA**](04-naca0012-cfd-structural/)                     | ⚪ Planned      | One-way aerodynamic pressure coupling to structural response      |
| 05 | [**Thermo-Structural Aerospace Panel**](05-thermo-structural-panel/)                  | ⚪ Planned      | Thermal gradient → thermal stress                                 |
| 06 | [**UAV Wing Spar Mass Optimization**](06-uav-wing-spar-optimization/)                 | ⚪ Planned      | Constrained mass minimization and structural trade study          |
| 07 | [**CubeSat Modal Analysis**](07-cubesat-modal-analysis/)                              | ⚪ Planned      | Natural frequencies and mode shapes                               |
| 08 | [**Thermal Protection / High-Temperature Structure**](08-high-temperature-structure/) | ⚪ Planned      | Aerothermal loading → thermal-structural response                 |

---

## Project Standards

Each project is expected to establish, where applicable:

* A clearly defined engineering problem
* Relevant structural theory and governing equations
* Explicit material and loading assumptions
* A justified structural model
* Appropriate boundary conditions
* Mesh sensitivity or convergence assessment
* Numerical verification
* Validation against analytical, experimental, or published data where credible references exist
* Parametric investigation where relevant
* Physical interpretation of structural behaviour
* Engineering conclusions and limitations

**A numerically solved model is not automatically a physically credible structural analysis.**

Verification and validation are treated separately. Where meaningful validation data are unavailable, that limitation is documented explicitly.

---

## Computational Stack

`CalculiX` · `Gmsh` · `OpenSCAD` · `Python` · `PyVista` · `ParaView` · `Git`

`OpenFOAM` is introduced in later projects where aerodynamic or aerothermal loading is coupled to the structural model.

---

## Reproducibility

Each project is maintained as a self-contained computational study.

Project-level documentation covers:

* Geometry and idealization
* Material properties
* Loading and boundary conditions
* Mesh generation
* Solver configuration
* Numerical settings
* Verification and validation approach
* Post-processing
* Reproduction procedure

The objective is for each study to be reproducible from a clean environment using the instructions provided in its repository.

---

## Repository Structure

Each project follows a structure similar to:

```text
NN-project-name/
├── README.md
├── geometry/
│   ├── cad/
│   └── exported/
├── mesh/
├── simulation/
├── scripts/
├── results/
└── figures/
```

The exact structure may vary depending on the requirements of the individual analysis.

---

## Status

* 🟢 **Completed** — project meets its documented completion standard
* 🟡 **In progress** — active development or investigation is underway
* ⚪ **Planned** — project has not yet begun

Being **In progress** does not imply that the model has been fully verified or validated.

A project is marked **Complete** only when the documented evidence supports its stated engineering conclusions and the reproducibility requirements have been satisfied.

---

## Direction

This portfolio provides the computational structures foundation for the broader aerospace research program.

The progression is:

**Structural mechanics → Aerospace structures → Structural optimization → Aerodynamic loading → Thermal loading → Fluid–structure interaction → Aerothermoelasticity**

The long-term goal is to develop the structural and computational mechanics capability required to study **coupled high-speed aerospace systems**, where aerodynamic, thermal, and structural physics interact.

---

## Author

**Harsh Panchal**
Aerospace Engineering · Computational Mechanics · Aerospace Structures
