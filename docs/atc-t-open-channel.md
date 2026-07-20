# Open-channel thermal adjoint validation status

This note records the public status of the open-channel co-optimisation example
and the separate fixed-point beta sensitivity path.

## Summary

`examples/coOptimiseChannelAndBody` is a generic open-channel
co-optimisation example. It is useful for exercising design-domain handling,
thermal objectives, pressure-drop monitoring and body-fitted verification
workflows. It is not used as quantitative evidence for an exact
variable-property SST adjoint.

The production exact path is separate. In `v0.4.2`,
`useAnalyticFixedPointAdjointBeta` enables the guarded matrix-free fixed-point
beta sensitivity for the constant-property, one-way-coupled serial Stokes
branch. That path solves

```text
(I - M_x)^T psi = J_x
```

and assembles

```text
dJ/dbeta = J_beta + psi^T M_beta.
```

The implementation replaces the stock Brinkman sensitivity with the validated
fixed-point flow derivative and retains the direct conductivity contribution.
It does not double-count the two terms.

## Validated Scope

The exact fixed-point beta mode is runtime-guarded to the branch that was
validated:

- OpenFOAM v2512;
- serial execution;
- non-coupled patches;
- Stokes flow;
- consistent SIMPLE;
- zero pressure non-orthogonal correctors;
- constant fluid and solid diffusivities;
- constant volumetric heat capacity;
- one-way flow-to-thermal coupling;
- one active `topOSource` momentum option;
- no active temperature `fvOptions`;
- zero explicit thermal non-orthogonal diffusion correction.

The guards are implemented in
`thermalAdjointSimple::validateProductionFixedPointBetaScope()` and
`thermalAdjointSimple::validateThermalDiffusionCorrectionScope()`, with
property-mode introspection supplied by `thermalPropertyTables`.

Structural lower-bound relaxation ties at `beta=0` are handled as feasible
right derivatives with respect to increasing Brinkman/source coefficient. They
are not described as general two-sided classical derivatives at the bound.

## Open-channel Interpretation

The open-channel example helped expose why a continuous thermal-to-momentum
source is not, by itself, a complete discrete transpose of the segregated SIMPLE
map in an open-flow topology optimisation setting. The validated production
path therefore does not rely on the continuous source alone. It uses the
fixed-point map transpose and the corresponding beta reverse.

For broader SST or variable-property searches, the current public claim remains
more limited: the search path uses frozen-coefficient thermal adjoints and
body-fitted re-verification. Exact nonlinear thermal-property derivatives and
full SST variable-property fixed-point sensitivities are future work.

## Reproducibility

The repository keeps the case files and diagnostic drivers used to reproduce the
verification campaigns. Headline public documentation states the scope and
result classes rather than embedding case-specific numerical tables. Exact
numeric rows belong in reproducible run logs, scripts and reviewer-requested
evidence packs.
