# coOptimiseChannelAndBody

**Co-optimisation of a coolant path and its surrounding conductive body, in one
design field.**

> **Experimental example.** The design-domain machinery below works and is
> reproducible. Quantitative gradient verification for this regime
> (open-channel flow) is still in progress — see [Verification status](#verification-status).
> Do not use its optimisation results as evidence. The verified capability lives
> in [`cases/`](../../cases).

## What it demonstrates

Most topology-optimisation examples designate the fluid region as the design
space and freeze the solid. This one puts the coolant path **and** the
surrounding conductive body in the **same** design field, so both co-evolve,
while a solid crust under the heated surface is pinned and the coolant can never
be sealed off from its inlet. That is the design-domain pattern a "fluid region
plus structural mass, minus a protected surface layer" problem needs.

```
 y=0.024  +-----------------------------------------+  heatedSurface (fixed flux in)
          |             crust  (FIXED SOLID)         |
 y=0.020  +-----------------------------------------+
          |          upper body  (DESIGN)            |
 y=0.016  +======+---------------------------+======+
   inlet ->      |      channel (DESIGN)      |      -> outlet
 y=0.008  +======+---------------------------+======+
          |          lower body  (DESIGN)            |
 y=0.004  +-----------------------------------------+
          |              base  (FIXED SOLID)         |
 y=0      +-----------------------------------------+  bottomWall (adiabatic)
        x=0    x=0.012                    x=0.108   x=0.12
               |<-- plug -->|         |<-- plug -->|      (FIXED FLUID)
```

Generic 2-D cooled slab. Nothing in it is taken from any specific component
geometry.

| Capability | How |
|---|---|
| designable fluid **and** solid body in one field | `adjointPorousZones (designSpace)` — 1536 cells spanning channel *and* body |
| protected solid crust under the heated surface | `fixedPorousZones (crustZone baseZone)` + `fixedPorousValues (1 1)` |
| inlet/outlet always open | `fixedZeroPorousZones (inletPlug outletPlug)` — pinned fluid |
| peak hot-wall temperature objective | `patchTemperaturePNorm` (P=8) on `heatedSurface` |
| pressure loss | `PtLosses` on (inlet outlet), monitored each cycle |
| no checkerboard / small-scale porosity | Helmholtz filter + tanh projection |
| **design history every cycle** | β iso-surface written **every** cycle to `optimisation/topOIsoSurfaces/` |
| **geometry export** | `.stl` *and* `.vtp` per cycle, straight from the framework |
| objective + constraint histories | `optimisation/objective/` (CSV-ready) |

Run: `./Allrun`.

## Confirmed working

- Primal and adjoint both converge to `residualControl`, not to an iteration cap.
- Sensitivities are confined to the design space; pinned zones stay pinned (β
  movement in fixed zones contributes 0.01–0.02 % of the sensitivity, and the
  thermal term is explicitly masked either side of the chain).
- The regularisation/projection chain is the exact transpose of the true β(α)
  map — verified numerically to 5 significant figures (`fd_chain.py`).
- The objective descends: the optimiser grows conductive structure into the
  coolant, as the physics predicts.
- Design history and geometry export work out of the box, every cycle.

## Verification status

The **sign** of the gradient is correct here and the design-domain controls above
are sound. What is not yet quantitatively verified is the **thermal-to-momentum
adjoint coupling in open-channel flow**.

The coupling term is the continuous form `Ta·∇T`. It is correct in the smooth,
low-velocity limit, which is why it verifies to **0.2–1.6 %** against finite
differences in every case in `cases/` — all of which place the design region in a
low-velocity Brinkman sponge, where the term contributes almost nothing. This
example is the first with an **open channel** carrying real velocity through the
design space, and there the continuous form is not the discrete transpose of the
primal convection operator.

The cause is understood and the exact replacement is derived and **proven**. With
`bounded Gauss <scheme>`, OpenFOAM assembles `div(phi,T) − Sp(div(phi),T)`, so at
convergence `div(phi) ≈ 0` leaves the primal *residual* unchanged but **not** its
derivative with respect to the face flux. The exact flux sensitivity

    gPhi_f = T_f (A_P − A_N) − χ (A_P T_P − A_N T_N),   A_P = C_P Ta_P,  χ = 1

matches a finite difference of the face Lagrangian to **10–11 significant
figures**. `utilities/testAdjointTranspose` is a no-solve gate that checks it, and
`couplingForm exactFluxTranspose` implements it. Remaining work: inject `gPhi_f`
through the framework's Rhie–Chow flux linearisation rather than as a cell
momentum source.

Full working note, with every hypothesis tested and rejected:
[`docs/atc-t-open-channel.md`](../../docs/atc-t-open-channel.md).

**None of this affects the verified capability in `cases/`**, which is unchanged
and independently reproducible.

## Gradient-verification harnesses (reusable)

| Script | What it does |
|---|---|
| `fd_check.py` | central FD vs adjoint on single design cells |
| `fd_sweep.py` | FD step-size sweep — finds the plateau where FD is a valid gate |
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
   like the optimiser going uphill. `bounded Gauss upwind` converges in 171
   iterations here.
2. **The adjoint flow solve wants the reverse relaxation pattern from the
   primal**: low on velocity (Ua 0.3, not 0.6), high on pressure. Otherwise the
   stock upstream `PtLosses` adjoint stalls at residual ~0.1 in stiff porous
   geometry.
