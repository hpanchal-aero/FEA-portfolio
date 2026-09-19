# Project 02 — 1U CubeSat Primary Structure

A staged FEA study of a 1U CubeSat primary structure under
quasi-static launch loading, built with CalculiX, Netgen, and Python.
Full development history lives in
[`development-log.md`](development-log.md); this file is the summary.

## Research Question

How does CubeSat structural configuration affect mass, stiffness,
and stress margin under representative launch loading?

## Methodology

**Status: specification locked, geometry not yet started.**

Two structural configurations of a 1U CubeSat frame will be compared
under an assumed quasi-static axial launch load:

- **Configuration A** — solid-panel frame (flat panels, corner rails,
  no lightening features)
- **Configuration B** — ribbed/lightened frame (same envelope, panels
  with lightening features), testing the mass-vs-stiffness trade-off

Material: Al 6061-T6. Load: 12g quasi-static axial (explicitly a
preliminary design assumption, not a sourced launch-vehicle
requirement — consistent with Project 01's treatment of its own load
assumption). Constraint: von Mises stress ≤ yield/1.5 FoS (~184 MPa).
Mass checked against the standard 1U CDS budget (1.33 kg) as a
design target.

Modal/frequency response is explicitly out of scope for this project
and deferred to Stage 7 (CubeSat modal analysis) per the portfolio's
original sequencing.

**Verification strategy**: unlike Project 01's bracket, no single
closed-form reference exists for the full assembled structure.
Individual features (a single panel as plate bending, a corner rail
as a simple column) will be verified analytically where possible
before combining; the full assembly is treated as an internal
consistency study (equilibrium, mesh convergence) rather than a
closed-form-verified result.

## Findings

Not yet available — project just started.

## Conclusion

Not yet available — project just started.
