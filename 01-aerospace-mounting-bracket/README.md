# Project 01 — Ribbed L-Bracket: Mass-Minimized Aerospace Mounting Bracket

A staged, verification-driven FEA study of an aerospace mounting
bracket, built with CalculiX, Netgen, and Python. Full development
history, every intermediate stage, and all diagnostic detail live in
[`development-log.md`](development-log.md); this file is the summary.

## Research Question

How can the mass of a ribbed L-bracket mounting interface be
minimized while satisfying stress, displacement, and factor-of-safety
constraints under a representative static equipment load?

## Methodology

The bracket geometry — a right-angle frame with a flange-junction
fillet, a stiffening gusset, a toe fillet, and a mounting hole — was
built up incrementally rather than analyzed all at once. Each
geometric feature (hole, fillet, gusset, bend) was introduced and
verified in isolation, against a closed-form or published reference
where one exists, before being combined into the final integrated
bracket. Once the full geometry was locked and verified (Stage 5),
three design variables (flange thickness, gusset leg length, toe
fillet radius) were swept parametrically via Latin Hypercube sampling
to search for the lowest-mass design satisfying a stress limit
(335 MPa, from an assumed Al 7075-T6 yield with a 1.5 factor of
safety) and a displacement limit (1.5mm) — both explicitly treated as
preliminary design assumptions, not sourced certification
requirements.

Toolchain: OpenSCAD/Gmsh/Netgen for geometry and meshing, CalculiX for
solving, Python (NumPy, PyVista, Matplotlib) for automation,
extraction, and visualization — all open-source, run from the
terminal, version-controlled with Git.

## Findings

- The verification pipeline (geometry → mesh → solve → analytical
  comparison) was validated stage by stage: a plain cantilever plate
  against Euler-Bernoulli beam theory, a hole against classical
  stress-concentration references, a shoulder fillet against
  Peterson's charts, and a two-flange frame against closed-form
  statics — each within a few percent of theory.
- An un-tapered gusset termination was found to introduce a genuine,
  non-convergent stress singularity; a fillet at that location was
  shown to resolve it, confirmed by a dedicated mesh-convergence
  study.
- In the fully integrated bracket, the governing stress location was
  the mounting hole edge, not the more heavily-instrumented gusset
  toe — a result only surfaced by inspecting the full stress field
  rather than a narrower, assumption-driven search.
- In the mass-minimization sweep, **zero of the 14 successfully
  solved fixed-resolution designs satisfied the displacement
  constraint.** Investigating why revealed that the two
  highest-thickness designs in the original sample — the ones most
  likely to be stiff enough to pass — had been excluded from the
  study for computational-resource reasons, not engineering ones. A
  targeted re-check confirmed both would have passed, and a follow-up
  sample in that region resolved a genuine feasible design space.

## Conclusion

The selected design — flange thickness 4.98mm, gusset leg 19.58mm,
toe fillet radius 2.83mm — is the lowest-mass configuration found
that satisfies both constraints (mass 31,137 mm³, 7% margin on
displacement, 44% margin on stress), verified mesh-convergent to
within 0.5% across two refinement levels. It is heavier than the
project's original, unconstrained baseline geometry: the real result
of this study is not that mass could be reduced from that baseline,
but that meeting the stated constraints requires more material than
an arbitrary starting design — and that finding this out relied on
questioning the study's own computational exclusions rather than
accepting them at face value. The search covered 34 sampled points
across the design space, not an exhaustive search, and the result is
reported as the best design found within that budget, not a proven
global optimum.
