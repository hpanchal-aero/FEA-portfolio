# Project 02 — 1U CubeSat Primary Structure

A staged FEA study of a 1U CubeSat primary structure under
quasi-static launch loading, built with CalculiX, Netgen, and Python.
Full development history lives in
[`development-log.md`](development-log.md); this file is the summary.

## Research Question

How does CubeSat structural configuration affect mass, stiffness,
and stress margin under representative launch loading?

## Methodology

**Status: Configuration A baseline complete and verified within the
limits stated below. Configuration B (ribbed/lightened) has not begun
as a study; it starts after Configuration A is finished.**

Two structural configurations of a 1U CubeSat frame will be compared
under an assumed quasi-static axial launch load:

- **Configuration A** — solid-panel frame: four 8.5 x 8.5 x 100 mm
  corner rails and four 2 mm x 83 mm side panels inside a
  100 x 100 x 100 mm envelope, fused into a single solid. Top and
  bottom panels are omitted, and bolted joints are not modeled (both
  documented simplifications).
- **Configuration B** — ribbed/lightened frame with the same
  envelope. Not started.

Material: Al 6061-T6 (E = 68,900 MPa, nu = 0.33, rho = 2700 kg/m^3).
Load: 12 g quasi-static axial, a preliminary design assumption and not
a sourced launch-vehicle requirement. A 1.33 kg total mass assumption
gives 156.6 N, applied as a uniform 0.1643 MPa pressure on the
953 mm^2 top face; equal axial strain in all members (common rigid
base, common length) gives equal axial stress. Base face (z = 0) fixed.
Constraint: von Mises <= yield / 1.5 (~184 MPa). Mass is checked
against the 1U CDS budget (1.33 kg) as a design target; that figure is
not verified against the current CDS text, and the yield value
(276 MPa) is a commonly quoted typical value with no source cited yet.

Modal and random-vibration response are out of scope and deferred to
Project 07 (CubeSat modal analysis).

**Toolchain**: geometry built with the Gmsh OCC kernel (geometry only,
no Gmsh meshing), meshed with Netgen (straight-sided C3D10), solved
with CalculiX 2.23, post-processed with Python (NumPy, PyVista).

**Verification strategy**: no single closed-form reference exists for
the full assembly. Individual members (an isolated panel and an
isolated rail) were verified against uniaxial (axial-bar) theory,
sigma = P/A and delta = PL/(AE), with closed-form buckling checks. The
assembly is an internal-consistency study: force equilibrium, two mesh
levels, and a roller-base diagnostic. It is not claimed as
closed-form-verified.

## Findings

Isolated members (1.5 mm mesh, fixed base, uniform top pressure):

| | Panel (2a) | Rail (2b) |
|---|---|---|
| Base reaction vs applied | +0.0000% | +0.0000% |
| Top displacement magnitude vs theory | 1.29% below | 0.33% below |
| Mean von Mises, 25 < z < 75 mm | 0.1638 MPa (-0.29%) | 0.1643 MPa (-0.001%) |
| Peak von Mises | 0.345 MPa, at base | 0.259 MPa, at base corners |
| Applied / critical load (closed form) | 0.0011 | 0.0016 |

Configuration A assembly (fixed base, C3D10):

| | 3.0 mm mesh | 2.0 mm mesh |
|---|---|---|
| Tets / equations | 31,013 / 186,228 | 83,854 / 474,207 |
| Base reaction vs 156.578 N | +0.0000% | -0.0000% |
| Top displacement magnitude vs theory | 1.71% below | 1.65% below |
| Mean von Mises, panel region | 0.16424 MPa | 0.16424 MPa |
| Mean von Mises, rail region | 0.16047 MPa | 0.16061 MPa |
| Peak von Mises (base corners) | 0.368 MPa | 0.422 MPa |

- Mean stress and displacement are **stable across the two mesh
  levels**. Two levels do not give a convergence order.
- The peak von Mises stress sits at the fixed-base corners. It rose
  15% on refinement and is **not converged and sensitive to the base
  boundary condition**. A 3.0 mm diagnostic with a roller base (axial
  motion of the base constrained, lateral motion free) gave a uniform
  0.1643 MPa stress everywhere and displacement matching theory to
  0.0002%, so the peak and the ~1.7% displacement offset come from the
  lateral restraint of the fixed base. The peak is not reported as a
  design stress.
- The mean stress (0.164 MPa) is about 1,100 times below the ~184 MPa
  allowable, so stress margin does not discriminate between
  configurations under this load case.

## Limitations

- **Axial-only loading (locked specification).** In a prismatic
  configuration such as Configuration A, stress is set by total
  cross-section area and stiffness by EA/L, so this load case cannot
  distinguish arrangements of equal area; for such variants the
  comparison is primarily a mass-versus-area trade. Lateral or
  combined quasi-static loading is a proposed extension, not part of
  the specification.
- Two mesh levels only; the roller-base run is a diagnostic at 3.0 mm.
- One element spans each 2 mm wall at 3.0 mm (checked in one panel
  region); local stress gradients at the rail-panel junction and the
  base are not resolved.
- Fully fused rails and panels: no bolted-joint contact or local joint
  stress.
- The rail-region mean stress is about 2% below the panel-region mean
  in the assembly; this is not explained.

## Conclusion

Configuration A is a verified baseline: equilibrium closes on both
meshes, mean stress and displacement match uniaxial theory to within
about 2% and are stable across two mesh levels, and the base-corner
peak is a boundary-condition effect that is not converged. Conclusions
about configuration trade-offs await the Configuration B study.
