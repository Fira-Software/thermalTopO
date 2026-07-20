# coOptimiseChannelAndBody

**Co-optimisation of a coolant path and its surrounding conductive body, in one
design field.**

> **Experimental example.** The design-domain machinery below works and is
> reproducible. This open-channel SST case is not used as quantitative evidence
> for an exact variable-property SST topology sensitivity. The verified
> production capability lives in [`cases/`](../../cases), and the exact
> fixed-point beta mode is a separately guarded constant-property Stokes path.

## What it demonstrates

Most topology-optimisation examples designate the fluid region as the design
space and freeze the solid. This one puts the coolant path **and** the
surrounding conductive body in the **same** design field, so both co-evolve,
while a solid crust under the heated surface is pinned and the coolant can never
be sealed off from its inlet. That is the design-domain pattern a "fluid region
plus structural mass, minus a protected surface layer" problem needs.

```
 top      +-----------------------------------------+  heatedSurface (fixed flux in)
          |             crust  (FIXED SOLID)         |
          +-----------------------------------------+
          |          upper body  (DESIGN)            |
          +======+---------------------------+======+
   inlet ->      |      channel (DESIGN)      |      -> outlet
          +======+---------------------------+======+
          |          lower body  (DESIGN)            |
          +-----------------------------------------+
          |              base  (FIXED SOLID)         |
 bottom   +-----------------------------------------+  bottomWall (adiabatic)
               |<-- plug -->|         |<-- plug -->|      (FIXED FLUID)
```

Generic 2-D cooled slab. Nothing in it is taken from any specific component
geometry.

| Capability | How |
|---|---|
| designable fluid **and** solid body in one field | `adjointPorousZones (designSpace)`, spanning channel *and* body |
| protected solid crust under the heated surface | `fixedPorousZones (crustZone baseZone)` + `fixedPorousValues (1 1)` |
| inlet/outlet always open | `fixedZeroPorousZones (inletPlug outletPlug)`, pinned fluid |
| peak hot-wall temperature objective | `patchTemperaturePNorm` (P=8) on `heatedSurface` |
| pressure loss | `PtLosses` on (inlet outlet), monitored each cycle |
| no checkerboard / small-scale porosity | Helmholtz filter + tanh projection |
| **design history every cycle** | β iso-surface written **every** cycle to `optimisation/topOIsoSurfaces/` |
| **geometry export** | `.stl` *and* `.vtp` per cycle, straight from the framework |
| objective + constraint histories | `optimisation/objective/` (CSV-ready) |

Run: `./Allrun`.

## Confirmed working

- Primal and adjoint both converge to `residualControl`, not to an iteration cap.
- Sensitivities are confined to the design space; pinned zones stay pinned, and
  the thermal term is explicitly masked either side of the chain.
- The regularisation/projection chain is the exact transpose of the true β(α)
  map, verified numerically by `fd_chain.py`.
- The objective descends: the optimiser grows conductive structure into the
  coolant, as the physics predicts.
- Design history and geometry export work out of the box, every cycle.

## Verification status

The **sign** of the gradient is correct here and the design-domain controls
above are sound. The quantitative production claim is intentionally narrower:
`useAnalyticFixedPointAdjointBeta` is validated for a constant-property,
one-way-coupled, serial Stokes branch with runtime scope guards. This example
uses the broader open-channel SST search path, where temperature-dependent
coefficient fields and turbulent thermal diffusivity are treated as frozen in
the thermal adjoint.

Full working note, with every hypothesis tested and rejected:
[`docs/atc-t-open-channel.md`](../../docs/atc-t-open-channel.md).

**None of this affects the verified capability in `cases/`**, which is unchanged
and independently reproducible.

## Gradient-verification harnesses (reusable)

| Script | What it does |
|---|---|
| `fd_check.py` | central FD vs adjoint on single design cells |
| `fd_sweep.py` | FD step-size sweep, finds the plateau where FD is a valid gate |
| `fd_chain.py` | chain-only FD: tests the filter/projection transpose with **no PDE physics** |
| `fd_localise.py` | ranks which cells actually drive a sensitivity |
| `fd_flowonly.py` | the same gate on stock OpenFOAM, with thermalTopO out of the path |
| `fd_ablation_{B,C}.py` | filter / projection ablations |

Two numerical rules learned here, applicable to any porous thermal topology
optimisation:

1. **Confirm `residualControl` was met, not just that the run finished.** With
   `limitedLinear` on `div(phi,T)`, T's outer residual can sit on a flat plateau
   while the linear solve converges happily: the run hits the `nIters` cap, the
   objective is still drifting, and across optimisation cycles that looks exactly
   like the optimiser going uphill. `bounded Gauss upwind` avoids that failure
   mode in this example.
2. **The adjoint flow solve wants the reverse relaxation pattern from the
   primal**: lower relaxation on velocity and higher relaxation on pressure.
   Otherwise the stock upstream `PtLosses` adjoint can stall in stiff porous
   geometry.
